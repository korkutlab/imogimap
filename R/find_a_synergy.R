#' Evalues synergies for a data frame containing values of an immune feature and stratified expression levels of two genes.

#' @importFrom stats median mad
#' @param fdata A numeric matrix with 3 columns: A numeric indicating value of an immune feature and two integers, each interpreted as a coded label for the expression level of a gene. Each row of fdata is a different sample or experiment.
#' @param method a character string indicating which synergy score method to be used. one of "max" or "independence". Default is "max".
#' @param oncogene1 An optional factor indicating assumed expression of gene1 in relation to IAP. One of "Expressed","Inhibited",or NA. Default is NA
#' @param oncogene2 An optional factor indicating assumed expression of  gene2. One of 1,-1,or NA. Default is NA
#' @keywords synergy scoring
#' @return A data frame containing: A numeric synergy score measuring synergistic impact of two genes on a feature, and two characters indicating whether expression or inhibition of each gene positively impacts feature. Returns NA if no synergistic interaction is found.
#' @details
#'
#'  fdata is a numeric matrix with 3 columns: The first column contains the numeric value of a single immune feature. The second and third columns are any of the coded labels 1 or 4, representing low or high expression levels of two genes as outputted by get_quantile_rank. Column-names should be names of feature and two genes accordingly. For details of quantile stratification see get_quantile_rank.
#'
#' For synergy score calculations, input data is stratified based on  its last two columns. Groups are labeled as low-low (LL), low-high (LH), high-low(HL), and high-high(HH) representing the expression levels of the corresponding two genes in each group. Groups are then reordered four times: (LL LH HL HH) (HH HL LH LL) (LH LL HH HL) and (HL HH LL LH). The right-most group is recognized as the base group and the left-most group is recognized as the impacted group for synergy score calculation. The order of groups in each permutation makes a unique assumption about whether each of the genes is an oncogene or a tumor suppressor (See publication for more details). If oncogene1 or oncogene2 values are provided, one of the four permutations will be chosen accordingly.
#'
#' Oncogene1 and oncogene2 force our prior knowledge on expression of single genes in in relation to IAP, and determine which groups should be taken as base-groups. If a gene is marked as "Expressed" it means it is expressed at higher IAP levels. If a gene is marked as "Inhibited" it indicates its lack of expression at higher IAP levels. Hence  groups with low(L) expression of gene are chosen as base groups (one of LL or LH) and only the two permutations with these base groups will be chosen.  Hence  groups with high(H) expression of gene are chosen as base groups (one of HH or LH) and only the two permutations with these base groups will be chosen.If Oncogene1 and oncogene2 are missing, all four permutations will be assessed as described above.
#'
#'
#' A synergy score is calculated for each permutation. The score with the maximum absolute value is chosen as the synergy score and the corresponding permutation order determines the oncogenic properties of two genes.  For details of synergy score calculations see calculate_syng_score.
#'
#'Three pvalues are calculated using Wilcoxon.test by comparing means of the target group with the other three groups. The maximum of the three pvalue is taken as the significance Pvalue of the score
#'
#' A data frame is returned that contains: Feature and genes name, a numeric synergy score measuring synergistic association of two genes on a feature, a numeric representing significance pvalue, and two characters indicating whether expression or inhibition of each gene positively associates with feature. If no synergistic association is found,score will be set to NA.
#'
#'
#' @examples
#' dft <- TCGA_Leukocyte_fraction$Leukocyte_fraction
#' mydata<- as.matrix(cbind(feature=dft,
#' gene1=sample(c(1,4),length(dft),replace = TRUE),
#' gene2=sample(c(1,4),length(dft),replace = TRUE)))
#' find_a_synergy(fdata=mydata,method="max")
#'
#' @seealso [get_quantile_rank()]
#' @export
#'
find_a_synergy=function(fdata,method,oncogene1,oncogene2){

  if(nrow(fdata)<2){
    sscoreij <- data.frame(agent1=colnames(fdata)[2],
                           agent2=colnames(fdata)[3],
                           Immune_feature=colnames(fdata)[1],
                           Synergy_score= NA,
                           agent1_expression = NA,
                           agent2_expression = NA,
                           wilcox_Pvalue=NA)
  }else{
    fdata <- cbind(fdata,as.integer((fdata[,2]*2+fdata[,3])/3))

    dft_LL <- fdata[fdata[,4]==1,1]
    dft_LH <- fdata[fdata[,4]==2,1]
    dft_HL <- fdata[fdata[,4]==3,1]
    dft_HH <- fdata[fdata[,4]==4,1]

    n1 <- length(dft_LL)
    n2 <- length(dft_LH)
    n3 <- length(dft_HL)
    n4 <- length(dft_HH)

    if(n1<1 || n2<1 || n3<1 || n4<1){
      CS <- NA
      agent1_expression <- NA
      agent2_expression <- NA
      pval <- NA
    }else{

      m_LL <- median(dft_LL,na.rm=TRUE)
      m_LH <- median(dft_LH,na.rm=TRUE)
      m_HL <- median(dft_HL,na.rm=TRUE)
      m_HH <- median(dft_HH,na.rm=TRUE)

      SEM2_LL <-  ( stats::mad(dft_LL,na.rm = T)^2 ) /sum( !is.na(dft_LL) )
      SEM2_LH <-  ( stats::mad(dft_LH,na.rm = T)^2 ) /sum( !is.na(dft_LH) )
      SEM2_HL <-  ( stats::mad(dft_HL,na.rm = T)^2 ) /sum( !is.na(dft_HL) )
      SEM2_HH <-  ( stats::mad(dft_HH,na.rm = T)^2 ) /sum( !is.na(dft_HH) )

      CS1 <- NA
      CS2 <- NA
      CS3 <- NA
      CS4 <- NA
      if(missing(oncogene1)) oncogene1<- NA
      if(missing(oncogene2)) oncogene2<- NA

      if(all(is.na(c(oncogene1,oncogene2)))){
        CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
        CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
        CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
        CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
      }else{
        if(!is.na(oncogene1) && is.na(oncogene2)){
          if(oncogene1=="Expressed"){
            CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
            CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
          }else{
            CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
            CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
          }
        }else{
          if(!is.na(oncogene2) && is.na(oncogene1)){
            if(oncogene2=="Expressed"){
              CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
              CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
            }else{
              CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
              CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
            }
          }else{
            if(oncogene1=="Expressed"){
              if(oncogene2=="Expressed"){
                CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
              }else{
                CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
              }
            }else{
              if(oncogene2=="Expressed"){
                CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
              }else{
                CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
              }
            }
          }
        }
      }


      fdata2 <- as.data.frame(fdata[,c(1,4)])
      colnames(fdata2) <- c("IAP","group")

      CS0 <- c(CS1,CS2,CS3,CS4)

      if(all(is.na(CS0))){
        CS <- NA
        agent1_expression <- NA
        agent2_expression <- NA
        pval <- NA
      }else{
        CS<-max(abs(CS0),na.rm = T)
        i <- which(abs(CS0)==CS)
        if(length(i)>1) i <- min(i)
        CS <- CS0[i]

        if(i==1){
          agent1_expression <- "Expressed"
          agent2_expression <- "Expressed"
          pval <-  max(ggpubr::compare_means(IAP~group,fdata2,paired=F,ref.group = 4)$p,na.rm=TRUE)
        }else{
          if(i==2){
            agent1_expression <- "Expressed"
            agent2_expression <- "Inhibited"
            pval <-  max(ggpubr::compare_means(IAP~group,fdata2,paired=F,ref.group = 3)$p,na.rm=TRUE)
          }else{
            if(i==3){
              agent1_expression <- "Inhibited"
              agent2_expression <- "Expressed"
              pval <-  max(ggpubr::compare_means(IAP~group,fdata2,paired=F,ref.group = 2)$p,na.rm=TRUE)
            }else{
              agent1_expression <- "Inhibited"
              agent2_expression <- "Inhibited"
              pval <- max(ggpubr::compare_means(IAP~group,fdata2,paired=F,ref.group = 1)$p,na.rm=TRUE)
            }
          }
        }
      }
    }

    sscoreij <- data.frame(agent1=colnames(fdata)[2],
                           agent2=colnames(fdata)[3],
                           Immune_feature=colnames(fdata)[1],
                           Synergy_score=CS,
                           agent1_expression=agent1_expression,
                           agent2_expression=agent2_expression,
                           wilcox_pvalue= pval )
  }
  return(sscoreij)
}
