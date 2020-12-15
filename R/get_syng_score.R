#' Calculate synergy score for the association effect of two genes/proteins on a feature or response
#'
#' @param fdata A numeric data frame with 3 columns: A numeric value of an immune feature and two group ID's based on quantile stratification of two genes/proteins.
#' @param method a charachter string indicating which synergy score to be used. one of "HSA" or "Bliss". Default is HSA.
#' @keywords synergy scoring, HSA, Bliss
#' @return A data frame containing synergy scores
#' @details
#'
#'  fdata is a numeric data frame with 3 columns: The first column contains the numeric value of a single feature. The second and third column are any of the numbers 1 or 4, representing low or high levels of expression of two genes/proteins respectely. For details of quantile startification see get_quantile_rank.
#'
#' Synergy score calculation are based on Bliss defenition of independence: First median values of the feature column are determined for four groups of low-low, low-high, high-low and high-high expressions. Lets define ma, mb, and mc as median of low-high, high-low and high-high groups minus median of low-low group. The synergy score is then defined as
#' floor(abs(sign(ma) + sign(mb) + sign(mc)) / 3)*(abs(mc) - max(abs(ma) , abs(mb))). The value of the score is positive for synergistic effects and negative for antagonistic effects. By default the synergy score will be zero if ma, mb , and mc have opposite signs.
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
    sscoreij <- data.frame(Gene=colnames(fdata)[2],
      ICP=colnames(fdata)[3],
      Immune_Feature=colnames(fdata)[1],
      CScore= NA)
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
      Cp1 <- NA
      Cp2 <- NA
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
        CS <- 0
      }else{
        ma <- abs(ma)
        mb <- abs(mb)
        mc <- abs(mc)
        method <-  toupper(method)

        if( (missing(method)) || (method=="HSA") ){
          CS <- CS_sign * (mc/max(ma , mb))
        }
        else{
          if(method=="BLISS"){
            CS <- CS_sign * (mc/(ma + mb - ma*mb))
          }else{
            stop("ERROR: Method is not  HSA or BLISS.")
          }
        }
      }
    }
    sscoreij <- data.frame(Co_target=colnames(fdata)[2],
      Immune_checkpoint=colnames(fdata)[3],
      Immune_feature=colnames(fdata)[1],
      Synergy_score=CS)
  }
  return(sscoreij)
}
