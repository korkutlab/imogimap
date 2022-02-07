#' Stratify each column of data using quantile values.
#' @importFrom  stats quantile
#' @param df S numeric dataframe
#' @keywords Quantile ranking, Quartiles, Median, dataframe, Matrix
#' @return  A numeric dataframes with quartile ranking group IDs
#' @details
#' The generic function get_quantile_rank stratifies each column of a numeric data frame or matrix based on quartile or median values of its distribution. By default, quartile values will be used to stratify data into groups 1,2,3 and 4, with 1 representing the group with values below the first quartile and 4 representing the group with values higher than the fourth quantile. When there are duplication in quartile values, median values will be automatically used for stratification. In this case data will be divided into two groups with group 1 containing values lower than or equal to the median and  group 4 containing values higher than median respectively.
#'
#' @examples get_quantile_rank(df=sample_mRNA_data)
#' @export

get_quantile_rank<-function(df){

  if(min(df)==max(df)) return(NA)

  dfr <-apply(df,2,function(x){
    df_breaks<- quantile(x, na.rm=TRUE,probs=c(0,0.25,0.5,0.75,1))
    if(!any(duplicated(df_breaks))){
      as.integer(cut(x,df_breaks,include.lowest=TRUE))
    }else{
      df_breaks<- quantile(x, na.rm=TRUE,probs=c(0,0.5,1))
      if(!any(duplicated(df_breaks))){
        3*as.integer(cut(x,df_breaks,include.lowest=TRUE))-2
      }else{
        df_breaks<- quantile(x, na.rm=TRUE,probs=c(0,1))
        x <- x-(df_breaks[1]-df_breaks[2])/2
        3*as.integer(cut(x,2,include.lowest=TRUE))-2
      }
    }
  })
  if(is.null(dfr)){
    lgenes <- names(dfr)[sapply(dfr, is.null)]
    warning(length(lgenes)," column(s) with duplications in quantile values are removed.")
    excluded_expressions<<-lapply(lgenes, function(x)paste0(x,"  "))
    return(dfr)
  }

  if(inherits(dfr, "list")){
    lgenes <- names(dfr)[sapply(dfr, is.null)]
    warning(length(lgenes)," column(s) with duplications in quantile values are removed.")
    excluded_expressions<<-lapply(lgenes, function(x)paste0(x,"  "))
    dfr <- as.data.frame(Filter(Negate(is.null), dfr))
  }
  rownames(dfr)<-rownames(df)
  dfr<- as.data.frame(dfr)
  dfr$Tumor_Sample_ID <- rownames(dfr)
  dfr <- dfr[,c("Tumor_Sample_ID",c(setdiff(colnames(dfr), "Tumor_Sample_ID")))]
  return(dfr)
}
