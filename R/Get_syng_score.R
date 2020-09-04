#' calculates a synergy score
#'
#' Generates coopretivity boxplots immune checkpoints
#' @param fdata a formatted dataframe
#' @keywords synergy
#' @return a dataframe with gene names and synergy scores
#' @export
#' @examples im_syng(fdata=myformatteddata)
#' Get_syng_score()

Get_syng_score=function(fdata){
  if(nrow(fdata)==0){
    sscoreij <- data.frame(Gene=colnames(fdata)[3],
      ICP=colnames(fdata)[4],
      Immune_Feature=colnames(fdata)[2],
      synergy_score= NA, pvalueA_AB=NA,pvalueB_AB=Na)
  }else{
    fdata$state <- as.integer((fdata[,3]*2+fdata[,4])/3)
    dft1 <- fdata[fdata$state==1,2]
    dft2 <- fdata[fdata$state==2,2]
    dft3 <- fdata[fdata$state==3,2]
    dft4 <- fdata[fdata$state==4,2]

    n1 <- length(dft1)
    n2 <- length(dft2)
    n3 <- length(dft3)
    n4 <- length(dft4)

    if(n1<2 || n2<2 || n3<2 || n4<2){
      ss <- NA
      sp1 <- NA
      sp2 <- NA
    }else{
      m1 <- median(dft1,na.rm=T)
      m2 <- median(dft2,na.rm=T)
      m3 <- median(dft3,na.rm=T)
      m4 <- median(dft4,na.rm=T)
      if(m1==0.0){
        ss <-  m2 * m3  - m4
      }else{
        ss <-  m2 * m3 / (m1 * m1) - (m4 / m1)
      }
      sp1 <- wilcox.test(dft3 , dft4 , paired = F , exact = F)$p.value
      sp2 <- wilcox.test(dft2 , dft4 , paired = F , exact = F)$p.value
    }
    sscoreij <- data.frame(Gene=colnames(fdata)[3],
      ICP=colnames(fdata)[4],
      Immune_Feature=colnames(fdata)[2],
      synergy_score=ss,pvalueA_AB=sp1,pvalueB_AB=sp2)
  }
  return(sscoreij)
}
