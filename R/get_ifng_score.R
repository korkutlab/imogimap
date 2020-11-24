#' Calculates Z score for IFNG
#' Calculates Z score, the mean-centered and normalized by standard deviation
#'     expression values for IFNG within the specified dataset
#' @param  pdf a formated mRNA data frame
#' @keywords immune checkpoints
#' @return a dataframe of correlation coefficient and p.values
#' @examples get_ifng_score(pdf=sample_mRNA_data)
#' @export
#'
get_ifng_score=function(pdf){
  IFNG_gene_list <- "IFNG"
  missing_IFNG <- IFNG_gene_list[-which(IFNG_gene_list %in% rownames(pdf))]
  if(length(missing_IFNG)>0){
    warning(length(missing_IFNG)," missing EMT signature genes:  \n  ",
      lapply(missing_IFNG, function(x)paste0(x,"  ")),"\nCheck EMT_gene_list for all signature genes.\n")
  }

  sub_list <- IFNG_gene_list[!(IFNG_gene_list %in% missing_IFNG)]
  pdf_sub <- t(as.matrix(pdf[rownames(pdf) %in% sub_list,]))
  if(nrow(pdf_sub)==0){
    warning("No signature genes found!")
    pscore <-tibble("score"=numeric())

  }else{
    pdf_sub <- log2(pdf_sub+1)
    pdf_sub <- as.data.frame( t(scale(t(pdf_sub), center = T, scale = T)))

    pscore <-as.data.frame(colMeans(pdf_sub,na.rm = T))
    colnames(pscore) <- "IFNGscore"
    pscore$IFNGscore <- 2^pscore$IFNGscore
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","IFNGscore")]
  }
  return(pscore)
}
