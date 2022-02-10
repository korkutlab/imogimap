#' Finds Spearman correlation between an oncogene, an immune checkpoint and immune associated phenotypes.
#' @importFrom dplyr bind_rows across mutate group_by everything distinct
#' @importFrom psych corr.test
#' @import curatedTCGAData
#' @importFrom tidyr gather
#' @param onco_gene  A character vector of gene Hugo symbols.
#' @param icp_gene An optional character vector of immune checkpoint gene/protein IDs.
#' @param cohort a character vector of TCGA diseases
#' @param sample_list An optional character vector of TCGA samples barcodes indicating a subset of samples within a cohort.
#' @keywords correlation, immuno-oncology features, gene expression
#' @return a list of dataframes containing Spearman correlations and non-FDR adjusted probability values.
#' @details
#'
#' im_cor_tcga uses NASeq2GeneNorm expression data, as provided by \code{\link[curatedTCGAData]{curatedTCGAData}}, to find correlation  between onco_genes and immune checkpoints and immuno-oncology features as listed in TCGA_immune_features_list.
#'
#' By default (if no icp_gene is specified), icp_gene_list will be used.
#'
#' For TCGA disease list see TCGA_disease_list
#'
#'All barcodes in sample_list must be 15 character long and belong to the same cohort. When sample_list is provided, cohort should be the disease cohort that they belong to, otherwise only the first element of the cohort list will be used.
#'
#' A non-FDR-adjusted p.value is reported for each correlation value to allow for easier adjustments by user.
#'
#'All barcodes in sample_list must be 15 character long and belong to the same cohort. When sample_list is provided, cohort should be the disease cohort that they belong to, otherwise only the first element of the cohort list will be used.
#'
#' @examples im_cor_tcga(onco_gene=c("BRAF"),icp_gene=c("CD274","CTLA4"),cohort=c("gbm"))
#' @export

im_cor_tcga<-function(onco_gene,icp_gene,cohort,sample_list){

  rho <- PATIENT_BARCODE <- padj <- TMB <- gene <- ict <- pvalue <- NULL

  if(!missing(sample_list)){
    cohort <- cohort[1]
  }
  cohort <- tolower(cohort)
  results <- list()

  for(cohortID in 1:length(cohort)){

    # Read data -----------------------
    df <-curatedTCGAData( diseaseCode = cohort[cohortID],version = "1.1.38",
      assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]
    df <- df@assays$data@listData[[1]]
    colnames(df)<-  substr(colnames(df), 1, 15)

    if(!missing(sample_list)){
      df<-df[,sample_list ]
      if(ncol(df==0)){
        stop("ERROR: barcodes not found.")
      }
    }
    if(length(onco_gene)==1){
      df_selected <- as.data.frame(df[rownames(df) %in% onco_gene,])
      colnames(df_selected) <- onco_gene
    }else{
      df_selected <- t(df[rownames(df) %in% onco_gene,,drop=F])
    }
    if(nrow(df_selected)==0){
      stop("ERROR: No valid Hugo symbol found")
    }

    # Calculate correlations with immune checkpoints -----------------------
    if(missing(icp_gene)){
      icp_gene <- icp_gene_list
    }
    df_icp <- t(df[rownames(df) %in% icp_gene,,drop=F])
    dft_cor <- corr.test(df_selected,df_icp,method="spearman",adjust = "none")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$onco_gene <- rownames(df_rho)
    df_rho <- gather(df_rho, icp_gene,rho,-onco_gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$onco_gene <- rownames(df_padj)
    df_padj <- gather(df_padj, icp_gene,padj,-onco_gene)
    icp_cor <- merge(df_rho,df_padj,by=c("onco_gene","icp_gene"))
    colnames(icp_cor)[1:2]<- c("onco_gene","Immune checkpoints")


    #MANDATORY: add ID column for the rest of the code-----------
    df_selected <- as.data.frame(df_selected)
    df_selected$Tumor_Sample_ID <- rownames(df_selected)


    # Calculate correlation with EMT score-----------------------
    cohort_EMT <- get_emt_score(df)
    df_EMT <- merge(df_selected,cohort_EMT,by="Tumor_Sample_ID")
    row.names(df_EMT) <- df_EMT$Tumor_Sample_ID
    dft1 <- as.matrix(df_EMT[,c(2:(ncol(df_EMT)-1))])
    dft2 <- as.matrix(df_EMT$EMTscore)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "none")
    EMT_cor <- data.frame("onco_gene"=onco_gene,"rho"=dft_cor$r,"padj"=dft_cor$p)

    # Calculate correlation with Angiogenesis score-----------------------
    cohort_AG <- get_angio_score(df)
    df_AG <- merge(df_selected,cohort_AG,by="Tumor_Sample_ID")
    row.names(df_AG) <- df_AG$Tumor_Sample_ID
    dft1 <- as.matrix(df_AG[,c(2:(ncol(df_AG)-1))])
    dft2 <- as.matrix(df_AG$AGscore)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "none")
    AG_cor <- data.frame("onco_gene"=onco_gene,"rho"=dft_cor$r,"padj"=dft_cor$p)

    # Calculate correlation with IFNG-----------------------
    cohort_IFNG <- get_ifng_score(df)
    df_IFNG <- merge(df_selected,cohort_IFNG,by="Tumor_Sample_ID")
    row.names(df_IFNG) <- df_IFNG$Tumor_Sample_ID
    dft1 <- as.matrix(df_IFNG[,c(2:(ncol(df_IFNG)-1))])
    dft2 <- as.matrix(df_IFNG$IFNGscore)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "none")
    IFNG_cor <- data.frame("onco_gene"=onco_gene,"rho"=dft_cor$r,"padj"=dft_cor$p)

    # Calculate correlation with Tumor mutation burden -----------------------
    df_TMB <- merge(df_selected,TCGA_TMB,by="Tumor_Sample_ID")
    row.names(df_TMB) <- df_TMB$Tumor_Sample_ID
    dft1 <- as.matrix(df_TMB[,c(2:(ncol(df_TMB)-2))])
    dft2 <- as.matrix(df_TMB[,c(ncol(df_TMB)-1,ncol(df_TMB))])
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "none")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$gene <- rownames(df_rho)
    df_rho <- gather(df_rho, TMB,rho,-gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$gene <- rownames(df_padj)
    df_padj <- gather(df_padj, TMB,padj,-gene)
    TMB_cor <- merge(df_rho,df_padj,by=c("gene","TMB"))
    TMB_cor <- data.frame("onco_gene"=onco_gene,"rho"=dft_cor$r,"padj"=dft_cor$p)

    # Calculate correlation with Leukocyte fraction -----------------------
    df_inf <- merge(df_selected,TCGA_Leukocyte_fraction,by="Tumor_Sample_ID")
    row.names(df_inf) <- df_inf$Tumor_Sample_ID
    dft1 <- as.matrix(df_inf[,c(2:(ncol(df_inf)-1))])
    dft2 <- as.matrix(df_inf$Leukocyte_fraction)
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "none")
    infiltration_cor <- data.frame("onco_gene"=onco_gene,"rho"=dft_cor$r,"padj"=dft_cor$p)


    # Calculate correlations with immune cell types -----------------------
    df_selected$PATIENT_BARCODE <- substr(df_selected$Tumor_Sample_ID, 1, 12)
    df_selected$Tumor_Sample_ID<-NULL
    df_selected <- df_selected %>% group_by(PATIENT_BARCODE) %>% mutate(across(cols = c(2),.fns = ~mean(.x, na.rm = TRUE))) %>% distinct
    df_ict <- merge(df_selected,TCGA_IMCell_fraction,by="PATIENT_BARCODE")
    row.names(df_ict) <- df_ict$PATIENT_BARCODE
    df_ict$PATIENT_BARCODE <- NULL
    df_selected <- df_selected[,c("PATIENT_BARCODE",c(setdiff(colnames(df_selected), "PATIENT_BARCODE")))]
    dft1 <- as.matrix(df_ict[,colnames(df_selected)[-1],drop=F])
    dft2 <- as.matrix(df_ict[,colnames(TCGA_IMCell_fraction)[-1],drop=F])
    dft_cor <- corr.test(dft1,dft2,method="spearman",adjust = "none")
    df_rho <- as.data.frame(dft_cor$r)
    df_rho$gene <- rownames(df_rho)
    df_rho <- gather(df_rho, ict,rho,-gene)
    df_padj <- as.data.frame(dft_cor$p)
    df_padj$gene <- rownames(df_padj)
    df_padj <- gather(df_padj, ict,pvalue,-gene)
    ict_cor <- merge(df_rho,df_padj,by=c("gene","ict"))
    colnames(ict_cor)[1:2]<- c("onco_gene","Immune cell types")

    # Output -----------------------
    results[[cohort[cohortID]]] <- c(list("Immune checkpoints"=icp_cor,"Immune cell types"=ict_cor,"Leukocyte fraction"=infiltration_cor,"EMT score"=EMT_cor, "Angiogenesis score"=AG_cor,"IFNG"=IFNG_cor))
  }
  return(results)
}
