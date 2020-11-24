#' Calculate Spearman correlation with immune checkpoints
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @param cotarget A charachter vector of gene/protein IDs.
#' @param checkpoint An optional charachter vector of immune checkpoint gene/protein IDs.
#' @param data_expression A numeric matrix or data frame containing gene/protein expressions.
#' @param data_feature An optional numeric matrix or data frame containing immune features.
#' @param add_features An optional logical indicating if EMT score, angiogenesis score and IFNG expression should be added to immune features. Default is TRUE.
#' @keywords Spearman correlation,immune features, immune checkpoints
#' @return A list of dataframes containing spearman correlations and non-FDR adjusted probability values.
#' @details
#'
#' By default (if no checkpoint is specified), icp_gene_list will be used.
#'
#' data_expression is formatted with genes/proteins as rows and samples/patients as columns.
#' For data_expression sample formats see sample_mRNA_data.
#'
#' data_feature is formated with samples/patients as rows and immune features as columns.
#' For data_feature sample format see sample_immune_cell_fraction_data.
#'
#' @examples im_cor(cotarget =  c("BRAF","CTLA4"),
#'                  checkpoint= c("CD274","CTLA4"),
#'                  data_expression =  sample_mRNA_data,
#'                  data_feature = sample_immune_cell_fraction_data,
#'                  add_features=T)
#' @export

im_cor<-function(cotarget,checkpoint, data_expression,data_feature,add_features){

  results <- list()
  data_expression <- as.data.frame(data_expression)

  #Check for co-target expressions----------------------------
  if(length(cotarget)==1){
    df_selected <- as.data.frame(t(data_expression[rownames(data_expression) %in% cotarget,]))
    colnames(df_selected) <- cotarget
  }else{
    df_selected <- t(data_expression[rownames(data_expression) %in% cotarget,])
  }
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }

  #Check for immune checkpoint expressions---------------------
  if(missing(checkpoint)){
    checkpoint <- icp_gene_list
  }
  missing_icp <- checkpoint[-which(checkpoint %in% rownames(data_expression))]
  if(length(missing_icp)>0){
    warning("Missing immune checkpoints:   ",lapply(missing_icp, function(x)paste0(x,"  ")))
  }

  #Calculate correlations with checkpoints---------------------
  df_icp <- t(data_expression[rownames(data_expression) %in% checkpoint,])
  if(nrow(df_icp)==0){
    warning("No immune checkpoint found.")
    icp_cor <- tibble("rho"=numeric(),"pvalue"=numeric())
  }else{
    dft_cor <- psych::corr.test(df_selected,df_icp,method="spearman",adjust = "none")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$cotarget <- rownames(df_rho)
    df_rho <- gather(df_rho, checkpoint,rho,-cotarget)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$cotarget <- rownames(df_padj)
    df_padj <- gather(df_padj, checkpoint,padj,-cotarget)
    icp_cor <- merge(df_rho,df_padj,by=c("cotarget","checkpoint"))
    colnames(icp_cor)[1:2]<- c("Co-targets","Immune checkpoints")
  }

  #Check for immune features------------------------------
  if(missing(data_feature)){
    data_feature <- data.frame(row.names = colnames(data_expression))
    if(missing(add_features)){
      add_features <- TRUE
    }else{
      if(add_features==FALSE){
        warning("No feature data is specified. Using optional features. \n")
        add_features <- TRUE
      }
    }
  }else{
    data_feature <- as.data.frame(data_feature)
    if(missing(add_features)){
      add_features <- TRUE
    }
  }
  data_feature$Tumor_Sample_ID <- rownames(data_feature)


  #Check for additional immune features---------------------
  if(add_features==T){
    df_EMT <-  get_emt_score(data_expression)
    if(nrow(df_EMT)==0){
      warning("No EMT signature marker found.\n")
    }else{
      data_feature <- merge(df_EMT,data_feature,by="Tumor_Sample_ID")
    }

    df_ang <-  get_angio_score(data_expression)
    if(nrow(df_ang)==0){
      warning("No angiogenesis signature marker found.\n")
    }else{
      data_feature <- merge(df_ang,data_feature,by="Tumor_Sample_ID")
    }

    df_ifng <-  as.data.frame(t(data_expression[rownames(data_expression)=="IFNG",]))
    if(nrow(df_ang)==0){
      warning("No IFNG expression found.\n")
    }else{
      df_ifng$Tumor_Sample_ID <- rownames(df_ifng)
      data_feature <- merge(df_ifng,data_feature,by="Tumor_Sample_ID")
    }
  }


  #Calculate correlation with immune features---------------

  df_selected <- as.data.frame(df_selected)
  df_selected$Tumor_Sample_ID <- rownames(df_selected)
  df_merged <- merge(df_selected,data_feature,by="Tumor_Sample_ID")
  if(nrow(df_merged)==0){
    warning("No data with cotarget expression and immune feature found")
    ifeature_cor <- tibble("rho"=numeric(),"pvalue"=numeric())
  }else{
    row.names(df_merged) <- df_merged$Tumor_Sample_ID
    strid <- as.numeric(which(colnames(df_selected)=="Tumor_Sample_ID"))
    dft1 <- as.matrix(select(df_merged,colnames(df_selected)[-strid]))
    dft2 <- as.matrix(select(df_merged,-colnames(df_selected)))
    dft_cor <- psych::corr.test(dft1,dft2,method="spearman",adjust = "none")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$gene <- rownames(df_rho)
    df_rho <- gather(df_rho, ict,rho,-gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$gene <- rownames(df_padj)
    df_padj <- gather(df_padj, ict,pvalue,-gene)
    ifeature_cor <- merge(df_rho,df_padj,by=c("gene","ict"))
    colnames(ifeature_cor)[1:2]<- c("Co-targets","Immune features")

  }


  results <- list("Immune checkpoint"=icp_cor,"Immune features"=ifeature_cor)

  return(results)
}
