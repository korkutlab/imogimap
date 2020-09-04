#' A correlation function
#'
#' This function allows you to calculate correlation with immune checkpoints
#' @param gene_list a list of genes Hugo symbol
#' @param cohort a list of TCGA diseases
#' @keywords immune checkpoints
#' @return a dataframe of correlation coefficient and p.values
#' @export
#' @examples im_cor_tcga(gene_list=c("TP53","HER2"),cohort=c("acc","gbm"))
#' im_cor_tcga()


im_cor_tcga<-function(gene_list,cohort){

  cohort <- tolower(cohort)
  results <- list()

  for(cohortID in 1:length(cohort)){

    # Read data -----------------------
    cohort_study <- paste0(cohort[cohortID],"_tcga_pan_can_atlas_2018")
    df <- cBioDataPack(cohort_study,ask = F)@ExperimentList@listData
    df <- df$RNA_Seq_v2_expression_median
    df2 <- df@elementMetadata@listData
    df <- df@assays@data@listData[[1]]
    rownames(df)<- mapvalues(rownames(df),df2$Entrez_Gene_Id,df2$Hugo_Symbol,
      warn_missing = F)
    if(length(gene_list)==1){
      df_selected <- as.data.frame(df[rownames(df) %in% gene_list,])
      colnames(df_selected) <- gene_list
    }else{
      df_selected <- t(df[rownames(df) %in% gene_list,])
    }
    if(nrow(df_selected)==0){
      stop("ERROR: No valid Hugo symbol found")
    }
    df_icp <- t(df[rownames(df) %in% icp_gene_list,])


    # Calculate correlations with immune checkpoints -----------------------
    dft_cor <- corr.test(df_selected,df_icp,method="spearman",adjust = "BH")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$test_gene <- rownames(df_rho)
    df_rho <- gather(df_rho, icp_gene,rho,-test_gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$test_gene <- rownames(df_padj)
    df_padj <- gather(df_padj, icp_gene,padj,-test_gene)
    icp_cor <- merge(df_rho,df_padj,by=c("test_gene","icp_gene"))
    colnames(icp_cor)[1:2]<- c("Genes","Immune checkpoints")

    df_selected <- as.data.frame(df_selected)
    df_selected$Tumor_Sample_ID <- rownames(df_selected)

    # Calculate correlation with EMT score-----------------------
    df_EMT <- merge(df_selected,TCGA_EMT,by="Tumor_Sample_ID")
    row.names(df_EMT) <- df_EMT$Tumor_Sample_ID
    dft1 <- as.matrix(df_EMT[,c(2:(ncol(df_EMT)-1))])
    dft2 <- as.matrix(df_EMT$EMTscore)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "BH")
    EMT_cor <- data.frame("Genes"=gene_list,"rho"=dft_cor$r,"padj"=dft_cor$p)


    # Calculate correlation with Leukocyte fraction -----------------------
    df_inf <- merge(df_selected,TCGA_Leukocyte_fraction,by="Tumor_Sample_ID")
    row.names(df_inf) <- df_inf$Tumor_Sample_ID
    dft1 <- as.matrix(df_inf[,c(2:(ncol(df_inf)-1))])
    dft2 <- as.matrix(df_inf$Leukocyte_fraction)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "BH")
    infiltration_cor <- data.frame("Genes"=gene_list,"rho"=dft_cor$r,"padj"=dft_cor$p)


    # Calculate correlations with immune cell types -----------------------
    df_selected$PATIENT_BARCODE <- base::substr(df_selected$Tumor_Sample_ID, 1, 12)
    df_selected$Tumor_Sample_ID<-NULL
    df_selected <- df_selected %>% group_by(PATIENT_BARCODE) %>% mutate(across(cols = everything(),.fns = ~mean(.x, na.rm = TRUE))) %>% distinct
    df_ict <- merge(df_selected,TCGA_IMCell_fraction,by="PATIENT_BARCODE")
    row.names(df_ict) <- df_ict$PATIENT_BARCODE
    df_ict$PATIENT_BARCODE <- NULL
    df_selected <- df_selected[,c("PATIENT_BARCODE",c(setdiff(colnames(df_selected), "PATIENT_BARCODE")))]
    dft1 <- as.matrix(select(df_ict,colnames(df_selected)[-1]))
    dft2 <- as.matrix(select(df_ict,colnames(ICT_fraction)[-1]))
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "BH")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$gene <- rownames(df_rho)
    df_rho <- gather(df_rho, ict,rho,-gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$gene <- rownames(df_padj)
    df_padj <- gather(df_padj, ict,padj,-gene)
    ict_cor <- merge(df_rho,df_padj,by=c("gene","ict"))
    colnames(ict_cor)[1:2]<- c("Genes","Immune cell types")

    # Output -----------------------
    results[[cohort[cohortID]]] <- c(list("checkpoint genes"=icp_cor,"Immune cell types"=ict_cor,"Leukocyte fraction"=infiltration_cor,"EMT score"=EMT_cor))
  }
  return(results)
}
