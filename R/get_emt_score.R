#' Calculates EMT score
#' Calculates EMT score from a gene list
#' @param  pdf a formated mRNA data frame
#' @keywords immune checkpoints
#' @return a dataframe of correlation coefficient and p.values
#' @examples get_emt_score(pdf=sample_mRNA_data)
#' @export
#'
#'
get_emt_score=function(pdf){

  missing_EMT <- EMT_gene_list[-which(EMT_gene_list$genes %in% rownames(pdf)),]$gene
  if(length(missing_EMT)>0){
    warning(length(missing_EMT)," missing EMT signature genes:  \n  ",
      lapply(missing_EMT, function(x)paste0(x,"  ")),"\nCheck EMT_gene_list for all signature genes.\n")
  }

  sub_list <- EMT_gene_list[!(EMT_gene_list$genes %in% missing_EMT),]
  pdf_sub <- pdf[rownames(pdf) %in% sub_list$genes,]
  if(nrow(pdf_sub)==0){
    warning("No signature genes found!")
    pscore <-tibble("score"=numeric())

  }else{
    pdf_sub <- log2(pdf_sub+1)
    pdf_sub <- as.data.frame( t(scale(t(pdf_sub), center = T, scale = T)))
    pdf_sub <- sweep(pdf_sub, 1, sub_list$sign, "*")

    pscore <-as.data.frame(colMeans(pdf_sub,na.rm = T))
    colnames(pscore) <- "EMTscore"
    pscore$EMTscore <- 2^pscore$EMTscore
    pscore$Tumor_Sample_ID <- rownames(pscore)
    pscore <- pscore[,c("Tumor_Sample_ID","EMTscore")]
  }
  return(pscore)
}
