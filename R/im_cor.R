#' A correlation function
#'
#' Calculates correlation with immune checkpoints
#' @param gene_list a list of charachters
#' @param df_mrna a formated mRNA data frame
#' @param df_ict a formated immune cell type fractions data frame
#' @param df_lf a formated Leukocyte fraction data frame
#' @keywords immune features
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @export
#' @examples im_cor(gene_list =  c("MCL1","TP53"),
#'                  df_mrna =  sample_mRNA_data,
#'                  df_lf = sample_Leukocyte_fraction_data,
#'                  df_ict = sample_immune_cell_fraction_data)
#' im_cor()


im_cor<-function(gene_list,df_mrna,df_ict,df_lf){

  results <- list()
  mydata <- as.data.frame(df_mrna)
  missing_icp <- icp_list[-which(icp_list %in% rownames(mydata))]

  if(length(missing_icp)>0){
    warning("Missing immune checkpoints:   ",lapply(missing_icp, function(x)paste0(x,"  ")))
  }

  if(length(gene_list)==1){
    df_selected <- as.data.frame(df[rownames(mydata) %in% gene_list,])
    colnames(df_selected) <- gene_list
  }else{
    df_selected <- t(mydata[rownames(mydata) %in% gene_list,])
  }
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }

  #Calculate immune checkpoint correlations
  df_icp <- t(mydata[rownames(mydata) %in% icp_list,])
  if(nrow(df_icp)==0){
    warning("No immune checkpoint found.")
    icp_cor <- tibble("rho"=numeric(),"padj"=numeric())
  }else{
    dft_cor <- corr.test(df_selected,df_icp,method="spearman",adjust = "BH")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$test_gene <- rownames(df_rho)
    df_rho <- gather(df_rho, icp_gene,rho,-test_gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$test_gene <- rownames(df_padj)
    df_padj <- gather(df_padj, icp_gene,padj,-test_gene)
    icp_cor <- merge(df_rho,df_padj,by=c("test_gene","icp_gene"))
    colnames(icp_cor)[1:2]<- c("Genes","Immune checkpoints")
  }
  #--------------------------------------------


  df_selected <- as.data.frame(df_selected)
  df_selected$Tumor_Sample_ID <- rownames(df_selected)
  #--------------------------------------------


  #Calculate EMT correlation
  df_EMT <-  Get_EMTscore(mydata)
  if(nrow(df_EMT)==0){
    warning("No EMT signature gene found.")
    EMT_cor <- tibble("rho"=numeric(),"padj"=numeric())
  }else{
    df_EMT <- merge(df_selected,df_EMT,by="Tumor_Sample_ID")
    row.names(df_EMT) <- df_EMT$Tumor_Sample_ID
    dft1 <- as.matrix(df_EMT[,c(2:(ncol(df_EMT)-1))])
    dft2 <- as.matrix(df_EMT$score)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "BH")
    EMT_cor <- data.frame("Genes"=gene_list,"rho"=dft_cor$r,"padj"=dft_cor$p)
  }
  #--------------------------------------------


  #Calculate infiltration correlation
  if(missing(df_lf)){
    warning("No Leukocyte fraction provided 1")
    infiltration_cor <- tibble("rho"=numeric(),"padj"=numeric())
  }else{
    df_lf <- as.data.frame(df_lf)
    colnames(df_lf)<-"Infiltration"
    df_lf$Tumor_Sample_ID <- rownames(df_lf)
    df_inf <- merge(df_selected,df_lf,by="Tumor_Sample_ID")
    if(nrow(df_inf)==0){
      warning("No Leukocyte fraction provided 2")
      infiltration_cor <- tibble("rho"=numeric(),"padj"=numeric())
    }else{

      row.names(df_inf) <- df_inf$Tumor_Sample_ID
      dft1 <- as.matrix(df_inf[,c(2:(ncol(df_inf)-1))])
      dft2 <- as.matrix(df_inf$Infiltration)
      dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "BH")
      infiltration_cor <- data.frame("Genes"=gene_list,"rho"=dft_cor$r,"padj"=dft_cor$p)
    }
  }
  #--------------------------------------------


  #Calculate immune cell fraction correlation
  if(missing(df_ict)){
    warning("No immune cell fractions provided 1")
    ict_cor <- tibble("rho"=numeric(),"padj"=numeric())
  }else{
    df_ict <- as.data.frame(df_ict)
    df_ict$Tumor_Sample_ID <- rownames(df_ict)
    df_ict <- merge(df_selected,df_ict,by="Tumor_Sample_ID")
    if(nrow(df_ict)==0){
      warning("No immune cell fractions provided 2")
      ict_cor <- tibble("rho"=numeric(),"padj"=numeric())
    }else{
      row.names(df_ict) <- df_ict$Tumor_Sample_ID
      strid <- as.numeric(which(colnames(df_selected)=="Tumor_Sample_ID"))
      dft1 <- as.matrix(select(df_ict,colnames(df_selected)[-strid]))
      dft2 <- as.matrix(select(df_ict,-colnames(df_selected)))
      dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "BH")
      df_rho <- as.data.frame(dft_cor$r)
      df_rho$gene <- rownames(df_rho)
      df_rho <- gather(df_rho, ict,rho,-gene)
      df_padj <- as.data.frame(dft_cor$p)
      df_padj$gene <- rownames(df_padj)
      df_padj <- gather(df_padj, ict,padj,-gene)
      ict_cor <- merge(df_rho,df_padj,by=c("gene","ict"))
      colnames(ict_cor)[1:2]<- c("Genes","Immune cell types")
    }
  }

  results <- list("Immune checkpoint"=icp_cor,"Immune cell types"=ict_cor,"Leukocyte fraction"=infiltration_cor,"EMT score"=EMT_cor)

  return(results)
}
