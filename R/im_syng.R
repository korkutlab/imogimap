#' A correlation function
#'
#' Generates coopretivity boxplots immune checkpoints
#' @param gene_list a list of charachters
#' @param df_mrna a formated mRNA data frame
#' @param df_ict a formated immune cell type fractions data frame
#' @param df_lf a formated Leukocyte fraction data frame
#' @keywords immune features
#' @return a dataframes of bliss scores
#' @export
#' @examples im_syng(gene_list =  c("TGBFR1","TGBFR2","ACVRL1","ENG"),
#'                   df_mrna =  sample_mRNA_data,
#'                   df_lf = sample_Leukocyte_fraction_data,
#'                   df_ict = sample_immune_cell_fraction_data)
#' im_syng()


im_syng<-function(gene_list,df_mrna,df_ict,df_lf){


  #read data
  mydata <- as.data.frame(df_mrna)
  if(length(gene_list)==1){
    df_selected <- as.data.frame(df[rownames(mydata) %in% gene_list,])
    colnames(df_selected) <- gene_list
  }else{
    df_selected <- t(mydata[rownames(mydata) %in% gene_list,])
  }
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }
  missing_icp <- icp_list[-which(icp_list %in% rownames(mydata))]
  if(length(missing_icp)>0){
    warning("Missing immune checkpoints:   ",lapply(missing_icp, function(x)paste0(x,"  ")))
  }
  df_icp <- t(mydata[rownames(mydata) %in% icp_list,])
  if(nrow(df_icp)==0){
    stop("No immune checkpoint found. Check icp_list for list of immune checkpoints")
  }

  df_EMT <-  Get_EMTscore(mydata)
  if(nrow(df_EMT)==0){
    warning("No EMT signature gene found.\n")
  }
  if(missing(df_lf)){
    warning("No Leukocyte fraction provided.\n")
  }else{
    df_lf2 <- as.data.frame(df_lf)
    colnames(df_lf2)<-"Leukocyte_fraction"
    df_lf2$Tumor_Sample_ID <- rownames(df_lf2)
  }
  if(missing(df_ict)){
    warning("No immune cell fractions provided.\n")
  }else{
    df_ict2 <- as.data.frame(df_ict)
    df_ict2$Tumor_Sample_ID <- rownames(df_ict2)
  }
  df_feature <- merge(df_EMT,df_lf2,by="Tumor_Sample_ID")
  df_feature <- merge(df_feature,df_ict2,by="Tumor_Sample_ID")
  #--------------------------------------------

  #construct quantile ranking matrices
  df_selected <- scale(log2(df_selected+1),center = T,scale = T)
  df_icp <- scale(log2(df_icp+1),center = T,scale = T)
  df_select_qr <- Get_qunatile_rank(df_selected)
  df_icp_qr <- Get_qunatile_rank(df_icp)
  #--------------------------------------------

  #Calculate synergy scores for EMT

  df_syng <- data.frame(Gene=character(),
    ICP=character(),
    Immune_Feature=character(),
    synergy_score=numeric(),
    pvalueA_AB=numeric(),
    pvalueB_AB=numeric())
  for(gene_ID in 2:ncol(df_select_qr)){
    for(icp_ID in 2:ncol(df_icp_qr)){
      dft <- merge(df_select_qr[,c(1,gene_ID)],df_icp_qr[,c(1,icp_ID)],by="Tumor_Sample_ID")
      for(if_ID in 2:ncol(df_feature)){
        dft2 <- merge(df_feature[,c(1,if_ID)],dft,by="Tumor_Sample_ID")
        dft2 <- dft2[dft2[,3] %in% c(1,4),]
        dft2 <- dft2[dft2[,4] %in% c(1,4),]
        df_syng <- rbind(df_syng,Get_syng_score(dft2))
      }
    }
  }

  return(df_syng)
}


