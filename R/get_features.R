#' Calculates and transforms default immune features within a specified sample cohort .
#' @param  pdf A numeric data frame or matrix with samples as columns and gene/proteins as rows.
#' @keywords Immune feature, Probability integral transform
#' @return A list with two dataframes containing transformed values for  immune features and immune cell type fractions as listed in \code{TCGA_immune_features_list}.
#' @details
#' Given a cohort of samples, calculates immune features within the cohort and transforms their values using probability integral transform. The resulting feature values will have a standard uniform distribution between [0,1]
#'
#' @examples get_features(pdf=sample_mRNA_data)
#' @export

get_features=function(pdf){

  cohort_EMT <- get_emt_score(pdf)
  cohort_IFNG <- get_ifng_score(pdf)
  cohort_AG <- get_angio_score(pdf)
  df_lf <-  TCGA_Leukocyte_fraction[
    TCGA_Leukocyte_fraction$Tumor_Sample_ID %in% colnames(pdf),]
  df_tmb <- TCGA_TMB[TCGA_TMB$Tumor_Sample_ID %in% colnames(pdf),]

  #Transform features-----------------------------
  #-----------------------------------------------
  if(nrow(cohort_EMT) > 0){
    sd_EMT <- sd(cohort_EMT$EMTscore,na.rm = T)
    cohort_EMT$EMTscore <- ( tanh( cohort_EMT$EMTscore/sd_EMT ) + 1 ) / 2
  }else{
    cohort_EMT<-data.frame("Tumor_Sample_ID"=colnames(pdf),"EMTscore"=NA)}

  if(nrow(cohort_IFNG) > 0){
    sd_IFNG <- sd(cohort_IFNG$IFNGscore,na.rm = T)
    cohort_IFNG$IFNGscore <- ( tanh( cohort_IFNG$IFNGscore/sd_IFNG ) + 1 ) / 2
  }else{
    cohort_IFNG<-data.frame("Tumor_Sample_ID"=colnames(pdf),"IFNGscore"=NA)}

  if(nrow(cohort_AG) > 0){
    sd_AG <- sd(cohort_AG$AGscore,na.rm = T)
    cohort_AG$AGscore <- ( tanh( cohort_AG$AGscore/sd_AG ) + 1 ) / 2
  }else{
    cohort_AG<-data.frame("Tumor_Sample_ID"=colnames(pdf),"AGscore"=NA)}

  if(nrow(df_tmb) > 0){
    df_tmb$TMB_Non.silent_per_Mb <- tanh(  df_tmb$TMB_Non.silent_per_Mb/10 )
    df_tmb$TMB_Silent_per_Mb <- tanh(  df_tmb$TMB_Silent_per_Mb/10 )
  }else{
    df_tmb<-data.frame("Tumor_Sample_ID"=colnames(pdf),"TMB_Non.silent_per_Mb"=NA,"TMB_Silent_per_Mb"=NA)}

  if(nrow(df_lf) == 0){df_lf <- data.frame("Tumor_Sample_ID"=colnames(pdf),"Leukocyte_fraction"=NA)}

  #Merge all features-----------------------------
  #-----------------------------------------------

  dft <- as.data.frame(merge(cohort_EMT , cohort_AG, by="Tumor_Sample_ID",all=T))
  dft <- as.data.frame(merge(dft , cohort_IFNG, by="Tumor_Sample_ID",all=T))
  dft<- as.data.frame(merge(dft , df_lf, by="Tumor_Sample_ID",all=T))
  dft <- as.data.frame(merge(dft , df_tmb, by="Tumor_Sample_ID",all=T))


  return(dft)
}
