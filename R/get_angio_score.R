#' Calculates Angiogenesis score
#' Calculates Angiogenesis score from a gene list
#' @param  pdf a formated mRNA data frame
#' @keywords immune checkpoints
#' @return a dataframe of correlation coefficient and p.values
#' @examples get_angio_score(pdf=sample_mRNA_data)
#' @export
#'
#'
get_angio_score=function(pdf){

  missing_genes <- AG_gene_list[-which(AG_gene_list %in% rownames(pdf))]
  if(length(missing_genes)>0){
    warning(length(missing_genes)," missing Angiogenesis signature genes:  \n  ",
      lapply(missing_EMT, function(x)paste0(x,"  ")),"\n Check AG_gene_list for
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
    pscore$AGscore <- 2^pscore$AGscore
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","AGscore")]
  }
  return(pscore)
}
