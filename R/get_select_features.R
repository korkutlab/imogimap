#' Calculates and transforms default immune features within a specified sample cohort .
#' @importFrom  stats sd
#' @param  pdf A numeric data frame or matrix with samples as columns and gene/proteins as rows.
#' @param  feature A list of character identifying one or more of the four expression based features: EMTscore, AGscore, TCIscore, or IFNGscore
#' @keywords Immune feature, Probability integral transform
#' @return A list with two dataframes containing transformed values for  immune features and immune cell type fractions as listed in \code{TCGA_immune_features_list}.
#' @details
#' Given a cohort of samples, calculates immune features within the cohort and transforms their values using probability integral transform. The resulting feature values will have a standard uniform distribution between \code{[0,1]}
#'
#' @examples get_selected_features(pdf=sample_mRNA_data,feature="EMTscore")
#' @export

get_selected_features=function(pdf,feature){
  
  dft <-data.frame("Tumor_Sample_ID"=colnames(pdf))

  feature_var <- feature[feature %like% "score"]
  
  if(length(feature_var)>0){
    
    if("EMTscore" %in% feature){
      cohort_EMT <- get_emt_score(pdf)
      if(nrow(cohort_EMT) > 0){
        sd_EMT <- sd(cohort_EMT$EMTscore,na.rm = T)
        cohort_EMT$EMTscore <- ( tanh( cohort_EMT$EMTscore/sd_EMT ) + 1 ) / 2
      }else{
        cohort_EMT<-data.frame("Tumor_Sample_ID"=colnames(pdf),"EMTscore"=NA)
      }
      dft <- as.data.frame(merge(dft,cohort_EMT , by="Tumor_Sample_ID",all=TRUE))
    }
    
    if("AGscore" %in% feature){
      cohort_AG <- get_angio_score(pdf)
      if(nrow(cohort_AG) > 0){
        sd_AG <- sd(cohort_AG$AGscore,na.rm = T)
        cohort_AG$AGscore <- ( tanh( cohort_AG$AGscore/sd_AG ) + 1 ) / 2
      }else{
        cohort_AG<-data.frame("Tumor_Sample_ID"=colnames(pdf),"AGscore"=NA)
      }
      dft <- as.data.frame(merge(dft,cohort_AG , by="Tumor_Sample_ID",all=TRUE))
    }
    if("TCIscore" %in% feature){
      cohort_TCI <- get_TCI_score(pdf)
      if(nrow(cohort_TCI) > 0){
        sd_TCI <- sd(cohort_TCI$TCIscore,na.rm = T)
        cohort_TCI$TCIscore <- ( tanh( cohort_TCI$TCIscore/sd_TCI ) + 1 ) / 2
      }else{
        cohort_TCI<-data.frame("Tumor_Sample_ID"=colnames(pdf),"TCIscore"=NA)
      }
      dft <- as.data.frame(merge(dft,cohort_TCI, by="Tumor_Sample_ID",all=TRUE))
    }
    
    if("IFNGscore" %in% feature){
      cohort_IFNG <- get_ifng_score(pdf)
      if(nrow(cohort_IFNG) > 0){
        sd_IFNG <- sd(cohort_IFNG$IFNGscore,na.rm = T)
        cohort_IFNG$IFNGscore <- ( tanh( cohort_IFNG$IFNGscore/sd_IFNG ) + 1 ) / 2
      }else{
        cohort_IFNG<-data.frame("Tumor_Sample_ID"=colnames(pdf),"IFNGscore"=NA)
      }
      dft <- as.data.frame(merge(dft,cohort_IFNG , by="Tumor_Sample_ID",all=TRUE))
    }
  }

  rownames(dft)<-dft$Tumor_Sample_ID
  dft$Tumor_Sample_ID <- NULL
  dft <- as.matrix(dft)
  
  return(dft)
}
