#' Calculates Angiogenesis score
#' @param  pdf A numeric data frame or matrix with samples as columns and gene/proteins as rows.
#' @keywords  Angiogenesis, pathway score
#' @return a dataframe of Angiogenesis scores for each sample
#' @details
#'  The angiogenesis score is defined as mean value of gene expressions in AG_gene_list. Scores are mean-centered and normalized by standard deviation based on distribution of scores in all samples.
#'
#' @examples get_angio_score(pdf=sample_mRNA_data)
#' @export
#'
#'
get_angio_score=function(pdf){
  pdf <- as.matrix(pdf)
  missing_genes <- AG_gene_list[-which(AG_gene_list %in% rownames(pdf))]
  if(length(missing_genes)>0){
    warning(length(missing_genes)," missing Angiogenesis signature genes:  \n  ",
      lapply(missing_genes, function(x)paste0(x,"  ")),"\n Check AG_gene_list for
      all signature genes.\n")
  }

  sub_list <- AG_gene_list[!(AG_gene_list %in% missing_genes)]
  pdf_sub <- pdf[rownames(pdf) %in% sub_list,]
  if(nrow(pdf_sub)==0){
    warning("No signature genes found!")
    pscore <-tibble("score"=numeric())

  }else{
    pdf_sub <- log2(pdf_sub+1)
    pdf_sub <- as.data.frame( t(scale(t(pdf_sub), center = T, scale = T)))

    pscore <-as.data.frame(colMeans(pdf_sub,na.rm = T))
    colnames(pscore) <- "AGscore"
    #pscore$AGscore <- 2^pscore$AGscore
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","AGscore")]
  }
  return(pscore)
}
