#' Finds a statistical synergistic score for a data frame containing values of an immune feature and stratified expression levels of two genes.
#'
#' @param fdata A numeric matrix with 3 columns: A numeric indicating value of an immune feature and two integers, each interpreted as a coded label for the expression level of a gene. Each row of fdata is a different sample or experiment.
#' @param method a character string indicating which synergy score method to be used. one of "max" or "independence". Default is "max".
#' @param oncogene1 A factor indicating oncogenic properties of gene1. One of 1,-1,or NA. Default is NA
#' @param oncogene2 A factor indicating oncogenic properties of gene2. One of 1,-1,or NA. Default is NA
#' @keywords synergy scoring
#' @return A data frame containing: A numeric synergy score measuring synergistic impact of two genes on a feature, and two characters indicating whether expression or inhibition of each gene positively impacts feature. Returns NA if no synergistic interaction is found.
#' @details
#'
#'  fdata is a numeric matrix with 3 columns: The first column contains the numeric value of a single immune feature. The second and third columns are any of the coded labels 1 or 4, representing low or high expression levels of two genes as outputted by get_quantile_rank. Column-names should be names of feature and two genes accordingly. For details of quantile stratification see get_quantile_rank.
#'
#' For synergy score calculations, input data is stratified based on  its last two columns. Groups are labeled as low-low (LL), low-high (LH), high-low(HL), and high-high(HH) representing the expression levels of the corresponding two genes in each group. Groups are then reordered four times: (LL LH HL HH) (HH HL LH LL) (LH LL HH HL) and (HL HH LL LH) and the right-most group is recognized as the base group for synergy score calculation. The order of groups in each permutation makes a unique assumption about whether each of the genes is an oncogene or a tumor suppressor (See publication for more details). If oncogene1 or oncogene2 values are provided, one of the four permutations will be chosen accordingly.
#'
#' Oncogene1 and oncogene2 force our prior knowledge on oncogenic properties of genes, and determine which groups should be taken as base-groups. For example, a value of 1 for oncogene1 means groups with low(L) expression of gene1 can be chosen as base groups (one of LL or LH) and only the two permutations with these base groups will be chosen. If Oncogene1 and oncogene2 are missing, all four permutations will be assessed as described above.
#'
#'
#' A synergy score is calculated for each permutation. The score with the maximum absolute value is chosen as the synergy score and the corresponding permutation order determines the oncogenic properties of two genes.  For details of synergy score calculations see calculate_syng_score.
#'
#' A data frame is returned that contains: Feature and genes name, a numeric synergy score measuring synergistic impact of two genes on a feature, and two characters indicating whether expression or inhibition of each gene positively impacts feature. If no synergistic impact is found,score will be set to NA.
#'
#'
#' @examples
#' dft <- TCGA_EMT$EMTscore
#' mydata<- as.matrix(cbind(feature=dft,
#' gene1=sample(c(1,4),length(dft),replace = T),
#' gene2=sample(c(1,4),length(dft),replace = T)))
#' find_a_synergy(mydata)
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
                           agent2_expression = NA)
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
    }else{

      m_LL <- median(dft_LL,na.rm=T)
      m_LH <- median(dft_LH,na.rm=T)
      m_HL <- median(dft_HL,na.rm=T)
      m_HH <- median(dft_HH,na.rm=T)

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
          if(oncogene1==1){
            CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
            CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
          }else{
            CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
            CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
          }
        }else{
          if(!is.na(oncogene2) && is.na(oncogene1)){
            if(oncogene2==1){
              CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
              CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
            }else{
              CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
              CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
            }
          }else{
            if(oncogene1==1){
              if(oncogene2==1){
                CS1 <- calculate_syng_score(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method)
              }else{
                CS2 <- calculate_syng_score(m_LH,m_LL,m_HH,m_HL,SEM2_LH,SEM2_LL,SEM2_HH,SEM2_HL,method)
              }
            }else{
              if(oncogene2==1){
                CS3 <- calculate_syng_score(m_HL,m_HH,m_LL,m_LH,SEM2_HL,SEM2_HH,SEM2_LL,SEM2_LH,method)
              }else{
                CS4 <- calculate_syng_score(m_HH,m_HL,m_LH,m_LL,SEM2_HH,SEM2_HL,SEM2_LH,SEM2_LL,method)
              }
            }
          }
        }
      }


      CS0 <- c(CS1,CS2,CS3,CS4)
      if(all(is.na(CS0))){
        CS <- NA
        agent1_expression <- NA
        agent2_expression <- NA
      }else{
        CS<-max(abs(CS0),na.rm = T)
        i <- which(abs(CS0)==CS)
        if(length(i)>1) i <- min(i)
        CS <- CS0[i]

        if(i==1){
          agent1_expression <- "Expressed"
          agent2_expression <- "Expressed"
        }else{
          if(i==2){
            agent1_expression <- "Expressed"
            agent2_expression <- "Inhibited"
          }else{
            if(i==3){
              agent1_expression <- "Inhibited"
              agent2_expression <- "Expressed"
            }else{
              agent1_expression <- "Inhibited"
              agent2_expression <- "Inhibited"
            }
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
                         agent2_expression=agent2_expression)

  return(sscoreij)
}
