#' A correlation function
#'
#' Generates coopretivity boxplots immune checkpoints
#' @param gene1 a single gene name
#' @param gene2 a single gene name
#' @param Immune_Feautue an immune feature name as listed in im_syng output.
#' @param df_mrna a formated mRNA data frame
#' @param df_ict an optional formated immune cell type fractions data frame
#' @param df_lf an optional formated Leukocyte fraction data frame
#' @keywords boxplots
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @export
#' @examples im_boot(gene1 = "TGFB1",gene2="TNFSF4",
#'                   Immune_Feature="EMTscore",
#'                   N_iteration=1000,
#'                   df_mrna =  sample_mRNA_data,
#'                   df_lf = sample_Leukocyte_fraction_data,
#'                   df_ict = sample_immune_cell_fraction_data)
#' im_boot()


im_boot<-function(gene1 , gene2 , Immune_Feature , N_iteration , df_mrna , df_ict , df_lf){


  #read data
  mydata <- as.data.frame(df_mrna)
  df_selected <- t(mydata[rownames(mydata) %in% c(gene1,gene2),])
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }
  if(Immune_Feature=="EMTscore"){
    df_EMT <-  Get_EMTscore(mydata)
    if(nrow(df_EMT)==0){
      stop("No EMT signature gene found.\n")
    }else{
      df_feature <- df_EMT
    }
  }else{
    if(Immune_Feature=="Leukocyte_fraction"){
      if(missing(df_lf)){
        stop("No Leukocyte fraction provided.\n")
      }else{
        df_lf2 <- as.data.frame(df_lf)
        colnames(df_lf2)<-"Leukocyte_fraction"
        df_lf2$Tumor_Sample_ID <- rownames(df_lf2)
        df_feature <- df_lf2
      }
    }else{
      if(missing(df_ict)){
        stop("No immune cell fractions provided.\n")
      }else{
        df_ict2 <- as.data.frame(df_ict)
        df_ict2$Tumor_Sample_ID <- rownames(df_ict2)
        tmp <- which(colnames(df_ict2)==Immune_Feature)
        if(length(tmp)==0){
          stop(Immune_Feature," Not found.\n")
        }else{
          tmpID <- which(colnames(df_ict2)=="Tumor_Sample_ID")
          df_feature <- df_ict2[,c(tmpID,tmp)]
        }
      }
    }
  }
  #--------------------------------------------

  #construct quantile ranking matrices
  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- Get_qunatile_rank(df_selected)
  if(is.null(df_select_qr)){
    stop('CScore calculation failed!' )
  }else{
    if(ncol(df_select_qr)<3){
      stop('CScore calculation failed!' )
    }
  }
  #--------------------------------------------


  #Get C Score

  dft <- merge(df_feature,df_select_qr,by="Tumor_Sample_ID")
  dft <- as.data.frame(dft)
  dft <- dft[dft[,3] %in% c(1,4),]
  dft <- dft[dft[,4] %in% c(1,4),]
  dft <- dft[complete.cases(dft),]
  myscore <- Get_CScore(dft)$CScore
  if(is.na(myscore)) stop("Error: Not enough data to calculate combination score!")
  mysign <- sign(myscore)
  #--------------------------------------------

  #Bootstrap
  P_Count <- 0.0
  df <- df_mrna

  for(i in 1:N_iteration){

    df_selected <- as.data.frame(t(df[sample(nrow(df),2,replace = F),]))
    df_selected <- scale(df_selected,center = T,scale = T)
    df_select_qr <- Get_qunatile_rank(df_selected)
    if(is.null(df_select_qr)){
      next
    }else{
      if(ncol(df_select_qr)<3){
        next
      }
    }
    dft <- merge(df_feature,df_select_qr,by="Tumor_Sample_ID")
    dft <- as.data.frame(dft)
    dft <- dft[dft[,3] %in% c(1,4),]
    dft <- dft[dft[,4] %in% c(1,4),]
    dft <- dft[complete.cases(dft),]
    cc <- sum(mysign*Get_CScore(dft)$CScore > mysign*myscore)
    if(!is.na(cc)){
      P_Count <- P_Count + cc
    }
  }

  P_Count <- P_Count/N_iteration

  return(P_Count)
}



