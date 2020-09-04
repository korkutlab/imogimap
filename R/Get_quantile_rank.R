#' A correlation function
#'
#' calculates quantile ranking for each column of a datarfame
#' @param df a numeric dataframe
#' @keywords quantile ranking
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @export
#' @examples im_coop( df=sample_mRNA_data)
#' Get_qunatile_rank()
#'
Get_qunatile_rank<-function(df){

  dfr <-apply(df,2,function(x){
    df_breaks<- quantile(x, na.rm=T,probs=c(0,0.25,0.5,0.75,1))
    if(!any(duplicated(df_breaks))){
      as.integer(cut(x,df_breaks,include.lowest=TRUE))
    }
  })
  if(inherits(dfr, "list")){
    lgenes <- names(dfr)[sapply(dfr, is.null)]
    warning("Genes with duplications in quantile values are removed:  ", lapply(lgenes, function(x)paste0(x,"  ")))
    dfr <- as.data.frame(Filter(Negate(is.null), dfr))
  }
  rownames(dfr)<-rownames(df)
  dfr<- as.data.frame(dfr)
  dfr$Tumor_Sample_ID <- rownames(dfr)
  dfr <- dfr[,c("Tumor_Sample_ID",c(setdiff(colnames(dfr), "Tumor_Sample_ID")))]
  return(dfr)
}
