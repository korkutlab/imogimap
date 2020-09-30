#' A correlation function
#'
#' This function allows you to calculate correlation with immune checkpoints
#' @param gene_list a list of genes Hugo symbol
#' @param cohort a list of TCGA diseases
#' @keywords immune checkpoints
#' @return a dataframe of correlation coefficient and p.values
#' @export
#' @examples im_syng_tcga(gene_list=c("TP53","TGFB1"),cohort=c("acc","gbm"))
#' im_syng_tcga()


im_syng_tcga<-function(gene_list,cohort){

  cohort <- tolower(cohort)

  df_syng <- data.frame(Cohort = character() ,
    Gene = character() , ICP = character() ,
    Immune_Feature = character() , synergy_score = numeric() ,
    pvalueA_AB = numeric() , pvalueB_AB = numeric())


  for(cohortID in 1:length(cohort)){

    # Read data -----------------------
    disease <- cohort[cohortID]
    cohort_study <- paste0(disease,"_tcga_pan_can_atlas_2018")
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

    # Construct quantile ranking matrices -----------------------
    df_selected <- scale(log2(df_selected+1),center = T,scale = T)
    df_icp <- scale(log2(df_icp+1),center = T,scale = T)
    df_select_qr <- Get_qunatile_rank(df_selected)
    df_icp_qr <- Get_qunatile_rank(df_icp)

    #Calculate synergy scores for EMT and Infiltration -----------------------
    df_feature <- as.data.frame(merge(TCGA_EMT , TCGA_Leukocyte_fraction, by="Tumor_Sample_ID"))


    for(gene_ID in 2:ncol(df_select_qr)){

      for(icp_ID in 2:ncol(df_icp_qr)){

        dft <- merge(df_select_qr[ , c(1 , gene_ID)] , df_icp_qr[ , c(1,icp_ID)] ,
          by="Tumor_Sample_ID")

        for(im_ID in 2:ncol(df_feature)){

          dft2 <- merge(df_feature[ , c(1 , im_ID)] , dft , by="Tumor_Sample_ID")
          dft2 <- dft2[dft2[ , 3] %in% c(1 , 4) , ]
          dft2 <- dft2[dft2[ , 4] %in% c(1 , 4) , ]
          dft2 <- dft2[complete.cases(dft2),]
          if(nrow(dft2>0)){
            dfts <- Get_CScore(dft2)
            dfts$Cohort <- disease
            dfts <- dfts[ , c("Cohort" , c(setdiff(colnames(dfts) , "Cohort")))]
            df_syng <- rbind(df_syng , dfts)
          }
        }
      }
    }

    #Calculate synergy scores for immune cell fractions -----------------------
    df_select_qr$PATIENT_BARCODE <- substr(df_select_qr$Tumor_Sample_ID, 1, 12)
    df_select_qr$Tumor_Sample_ID <- NULL
    df_select_qr <- df_select_qr %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct

    df_icp_qr$PATIENT_BARCODE <- substr(df_icp_qr$Tumor_Sample_ID, 1, 12)
    df_icp_qr$Tumor_Sample_ID <- NULL
    df_icp_qr <- df_icp_qr %>% group_by(PATIENT_BARCODE) %>%
      mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct

    df_select_qr <- df_select_qr[,c("PATIENT_BARCODE",
      c(setdiff(colnames(df_select_qr), "PATIENT_BARCODE")))]
    df_icp_qr <- df_icp_qr[,c("PATIENT_BARCODE",
      c(setdiff(colnames(df_icp_qr), "PATIENT_BARCODE")))]

    for(gene_ID in 2:ncol(df_select_qr)){
      for(icp_ID in 2:ncol(df_icp_qr)){
        dft <- merge(df_select_qr[,c(1,gene_ID)],df_icp_qr[,c(1,icp_ID)],by="PATIENT_BARCODE")
        for(if_ID in 2:ncol(TCGA_IMCell_fraction)){
          dft2 <- merge(TCGA_IMCell_fraction[,c(1,if_ID)],dft,by="PATIENT_BARCODE")
          dft2 <- dft2[dft2[,3] %in% c(1,4),]
          dft2 <- dft2[dft2[,4] %in% c(1,4),]
          dft2 <- dft2[complete.cases(dft2),]
          if(nrow(dft2>0)){
            dfts <- Get_CScore(dft2)
            dfts$Cohort <- disease
            dfts <- dfts[,c("Cohort",c(setdiff(colnames(dfts),"Cohort")))]
            df_syng <- rbind(df_syng,dfts)
          }
        }
      }
    }
  }
  return(df_syng)
}
