#' Calculates a score for IFNG expression
#' @param  pdf A numeric data frame or matrix with samples as columns and gene/proteins as rows.
#' @keywords IFNG, gene/protein expression
#' @return a dataframe of IFNG score for each sample
#' @details
#' To calculate IFNG score, IFNG expression is mean-centered and normalized by standard deviation based on distribution of  expression in all samples.
#'
#' @examples get_ifng_score(pdf=sample_mRNA_data)
#' @export
#'
get_ifng_score=function(pdf){
  pdf <- as.matrix(pdf)
  IFNG_gene_list <- "IFNG"
  missing_IFNG <- IFNG_gene_list[-which(IFNG_gene_list %in% rownames(pdf))]
  if(length(missing_IFNG)>0){
    warning(length(missing_IFNG)," missing IFNG genes:",
      lapply(missing_IFNG, function(x)paste0(x,"  ")),"\n")
  }

  sub_list <- IFNG_gene_list[!(IFNG_gene_list %in% missing_IFNG)]
  pdf_sub <- as.matrix(pdf[rownames(pdf) %in% sub_list,])
  if(nrow(pdf_sub)==0){
    warning("No signature genes found!")
    pscore <-tibble("score"=numeric())

  }else{
    pdf_sub <- log2(pdf_sub+1)
    pscore <- as.data.frame( scale((pdf_sub), center = T, scale = T) )

    colnames(pscore) <- "IFNGscore"
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","IFNGscore")]
  }
  return(pscore)
}
