#' finds a statistical synergistic score for the combinatorial association of gene pairs on a feature
#'
#' @param fdata A numeric data frame with 3 columns: A numeric indicating value of an immune feature and two integers, each interpreted as a coded label for the expression level of a gene. Each row of fdata is a different sample or experiment.
#' @param method a charachter string indicating which synergy score method to be used. one of "max" or "independence". Default is "max".
#' @keywords synergy scoring
#' @return A data frame containing synergy scores
#' @details
#'
#'  fdata is a numeric data frame with 3 columns: The first column contains the numeric value of a single immune feature. The second and third columns are any of the coded labels 1 or 4, representing low or high expression levels of two genes as outputed by get_quantile_rank. For details of quantile startification see get_quantile_rank.
#' For synergy score calculations, fdata is stratified based on coded labels of the last two columns. Groups are labled as low-low, low-high, high-low and high-high representing the expression levels of the corresponding two genes in each group. Median values of the feature column are then determined for each group.
#'
#' Synergy scores are defined as the deviation of immune feature's median value in high-high group from the its expected median as estimated using one of the "max" or "independence" models. Lets define ma and mb as the difference between the levels of immune feature in low-low group and low-high/high-low groups respectively. The expected median of high-high group is then defined as max(ma , mb) in "max" model and as ma+mb-ma*mb in "independence" model. Synergy score is then calculated as the difference between mc and the expected median. Scores are called synergistic if mc exceeds the expected median. Antagonistic scores are ignored. The value of the synergy score are multipled by 100 and a sign factor depending on the dirrection of the effect. By default synergy score will be zero if median values change in opposite directions. Differences in median values that are smaller that combined standard errors are considered as zero.
#'
#'
#'
#' @examples
#' dft <- TCGA_EMT$EMTscore
#' mydata<- cbind(feature=dft,gene1=sample(c(1,4),length(dft),replace = T),gene2=sample(c(1,4),length(dft),replace = T))
#' get_syng_score(mydata)
#' @export
#'
get_syng_score=function(fdata,method){

  if(nrow(fdata)<2){
    sscoreij <- data.frame(agent1=colnames(fdata)[2],
      agent2=colnames(fdata)[3],
      Immune_Feature=colnames(fdata)[1],
      Synergy_score= NA)
  }else{
    fdata$state <- as.integer((fdata[,2]*2+fdata[,3])/3)
    dft_LL <- fdata[fdata$state==1,1]
    dft_LH <- fdata[fdata$state==2,1]
    dft_HL <- fdata[fdata$state==3,1]
    dft_HH <- fdata[fdata$state==4,1]

    n1 <- length(dft_LL)
    n2 <- length(dft_LH)
    n3 <- length(dft_HL)
    n4 <- length(dft_HH)

    if(n1<1 || n2<1 || n3<1 || n4<1){
      CS <- NA
    }else{
      m_LL <- median(dft_LL,na.rm=T)
      m_LH <- median(dft_LH,na.rm=T)
      m_HL <- median(dft_HL,na.rm=T)
      m_HH <- median(dft_HH,na.rm=T)

      ma <- (m_LH - m_LL)
      mb <- (m_HL - m_LL)
      mc <- (m_HH - m_LL)

      SEM2_LL <-  ( stats::mad(dft_LL,na.rm = T)^2 ) /sum( !is.na(dft_LL) )
      SEM2_LH <-  ( stats::mad(dft_LH,na.rm = T)^2 ) /sum( !is.na(dft_LH) )
      SEM2_HL <-  ( stats::mad(dft_HL,na.rm = T)^2 ) /sum( !is.na(dft_HL) )
      SEM2_HH <-  ( stats::mad(dft_HH,na.rm = T)^2 ) /sum( !is.na(dft_HH) )

      SEM_a <- sqrt( SEM2_LH + SEM2_LL )
      SEM_b <- sqrt( SEM2_HL + SEM2_LL )
      SEM_c <- sqrt( SEM2_HH + SEM2_LL )

      sign_a <- ifelse(sign(ma+SEM_a)==sign(ma-SEM_a),sign(ma),NA)
      sign_b <- ifelse(sign(mb+SEM_b)==sign(mb-SEM_b),sign(mb),NA)
      sign_c <- ifelse(sign(mc+SEM_c)==sign(mc-SEM_c),sign(mc),NA)

      my_sign <- c(sign_a,sign_b,sign_c)
      CS_sign <- sign_c*floor( abs(mean(my_sign,na.rm=T)))

      if(CS_sign==0 || is.na(CS_sign)){
        CS <- NA
      }else{
        ma <- abs(ma)
        mb <- abs(mb)
        mc <- abs(mc)

        if( missing(method) ){
          method=="max"
        }

        if(method=="max" ){
          me <- mc - max(ma,mb)
        }else{
          if(method=="independence"){
            me <- mc-(ma + mb - ma*mb)
          }else{
            stop("ERROR: Method is not found. Please choose a method from:  max or independence.")
          }
        }
        if( me < 0 ){
          CS<- NA
        }else{
          CS <- CS_sign * me
        }
      }
    }
    sscoreij <- data.frame(agent1=colnames(fdata)[2],
      agent2=colnames(fdata)[3],
      Immune_feature=colnames(fdata)[1],
      Synergy_score=CS)
  }
  return(sscoreij)
}
