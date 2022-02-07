#'Calculate four synergy scores given median and standard error of means of four groups
#' @param m_LL median in Low_Low expression group
#' @param m_LH median in Low_High expression group
#' @param m_HL median in High_Low expression group
#' @param m_HH median in High_High expression group
#' @param SEM2_LL Standard error of the mean of the representing group
#' @param SEM2_LH Standard error of the mean of the representing group
#' @param SEM2_HL Standard error of the mean of the representing group
#' @param SEM2_HH Standard error of the mean of the representing group
#' @param method A character string indicating which synergy score to be used. one of "max" or "independence".
#' Default is "max".
#' @keywords synergy score
#' @return A numerical  scores
#' @details
#'
#' Calculates a synergy score given median and standard error of means of four groups each representing different expression levels of two genes.
#'
#' Synergy scores are defined as the deviation of median value in high-high group from the its expected median as estimated using one of the "max" or "independence" models. Lets define ma and mb as the difference between the levels of immune feature in low-low group and low-high/high-low groups respectively. The expected median of high-high group is then defined as max(ma , mb) in "max" model and as ma+mb-ma*mb in "independence" model. Synergy score is then calculated as the difference between mc and the expected median. Scores are called synergistic if mc exceeds the expected median. Antagonistic scores are ignored. The value of the synergy score are multiplied by 100 and a sign factor depending on the direction of the effect. By default synergy score will be zero if median values change in opposite directions. Differences in median values that are smaller that combined standard errors are considered as zero.
#'@export

calculate_syng_score=function(m_LL,m_LH,m_HL,m_HH,SEM2_LL,SEM2_LH,SEM2_HL,SEM2_HH,method){

  ma <- (m_LH - m_LL)
  mb <- (m_HL - m_LL)
  mc <- (m_HH - m_LL)

  SEM_a <- sqrt( SEM2_LH + SEM2_LL )
  SEM_b <- sqrt( SEM2_HL + SEM2_LL )
  SEM_c <- sqrt( SEM2_HH + SEM2_LL )

  if(sign(ma+SEM_a)==sign(ma-SEM_a)){
    sign_a <- sign(ma)
  }else{
    sign_a <- NA}
  if(sign(mb+SEM_b)==sign(mb-SEM_b)){
    sign_b <- sign(mb)
  }else{
    sign_b <- NA}
  if(sign(mc+SEM_c)==sign(mc-SEM_c)){
    sign_c <- sign(mc)
  }else{
    sign_c <- NA}

  my_sign <- c(sign_a,sign_b,sign_c)
  CS_sign <- sign_c*floor( abs(mean(my_sign,na.rm=T)))

  if(CS_sign==0 || is.na(CS_sign)){
    CS <- NA
  }else{
    ma <- abs(ma)
    mb <- abs(mb)
    mc <- abs(mc)
    if(method=="max" ){
      me <- mc - max(ma,mb)
    }else{
      me <- mc - (ma + mb - ma*mb)
    }
    if( me < 0 ){
      CS<- NA
    }else{
      CS <- CS_sign * me
    }
  }
  return(CS)
}
