#' A bootstrapping function
#'
#' Generates a probabilistic measure for an observation useing bootstrapping
#' @param gene1 a single gene name
#' @param gene2 a single gene name
#' @param cohort a single TCGA disease
#' @param Immune_Feature an immune feature name as listed in im_syng_tcga output.
#' @keywords pvalue, bootstrap, random sampling
#' @return the probability of observation
#' @examples im_boot_tcga(gene1 = "ACVRL1",gene2="CD274",
#'                   cohort="brca",
#'                   Immune_Feature="Leukocyte_fraction",
#'                   N_iteration=100000)
#' @export

im_boot_tcga<-function(gene1,gene2,cohort,Immune_Feature, N_iteration){


  cohort <- tolower(cohort)


  # Read data -----------------------

  cohort_study <- paste0(cohort,"_tcga_pan_can_atlas_2018")
  df <- cBioDataPack(cohort_study,ask=F)@ExperimentList@listData
  df <- df$RNA_Seq_v2_expression_median
  df2 <- df@elementMetadata@listData
  df <- df@assays@data@listData[[1]]
  rownames(df)<- mapvalues(rownames(df),df2$Entrez_Gene_Id,df2$Hugo_Symbol,
    warn_missing = F)
  df <- df[-which(is.na(rownames(df))),]

  #Select genes
  df_selected <- as.data.frame(t(df[rownames(df) %in% c(gene1,gene2),]))
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }

  #Select Feature
  if(Immune_Feature=="EMTscore"){
    df_feature <- TCGA_EMT
    colnames(df_feature)[1] <- "PATIENT_BARCODE"
  }else{
    if(Immune_Feature=="Leukocyte_fraction"){
      df_feature <- TCGA_Leukocyte_fraction
      colnames(df_feature)[1] <- "PATIENT_BARCODE"
    }else{
      if(Immune_Feature=="AGscore"){
        df_feature <- get_angio_score(df)
        colnames(df_feature)[1] <- "PATIENT_BARCODE"
      }else{
        if(Immune_Feature=="IFNGscore"){
          df_feature <- get_ifng_score(df)
          colnames(df_feature)[1] <- "PATIENT_BARCODE"
        }else{
          if(Immune_Feature %like% "TMB"){
            df_feature <- TCGA_TMB[,c("Tumor_Sample_ID",Immune_Feature)]
            colnames(df_feature)[1] <- "PATIENT_BARCODE"
          }else{
            df_feature <- TCGA_IMCell_fraction
            tmp <- which(colnames(df_feature)==Immune_Feature)
            if(length(tmp)==0){
              stop(Immune_Feature," Not found. Check ... for list of immune features\n")
            }else{
              tmpID <- which(colnames(df_feature)=="PATIENT_BARCODE")
              df_feature <- df_feature[,c(tmpID,tmp)]
              df_selected$PATIENT_BARCODE <- substr(rownames(df_selected), 1, 12)
              df_selected <- df_selected %>% group_by(PATIENT_BARCODE) %>%
                mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
              df_selected <- as.data.frame(df_selected)
              rownames(df_selected) <- df_selected$PATIENT_BARCODE
              df_selected$PATIENT_BARCODE <- NULL
            }
          }
        }
      }
    }
  }
  #--------------------------------------------

  #construct quantile ranking matrices

  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- get_qunatile_rank(df_selected)
  if(is.null(df_select_qr)){
    stop('CScore calculation failed!' )
  }else{
    if(ncol(df_select_qr)<3){
      stop('CScore calculation failed!' )
    }
  }
  colnames(df_select_qr)[1]<- "PATIENT_BARCODE"
  #--------------------------------------------

  #Get C Score

  dft <- merge(df_feature,df_select_qr,by="PATIENT_BARCODE")
  dft <- as.data.frame(dft)
  dft <- dft[dft[,3] %in% c(1,4),]
  dft <- dft[dft[,4] %in% c(1,4),]
  dft <- dft[complete.cases(dft),]
  dft<- dft[,c(2,3,4)]
  myscore <- get_syng_score(dft)$CScore
  mysign <- sign(myscore)
  #--------------------------------------------

  #Bootstrap
  P_Count <- 0.0
  if(Immune_Feature %in% colnames(TCGA_IMCell_fraction)){

    for(i in 1:N_iteration){

      df_selected <- as.data.frame(t(df[sample(nrow(df),2,replace = F),]))
      df_selected$PATIENT_BARCODE <- substr(rownames(df_selected), 1, 12)
      df_selected <- df_selected %>% group_by(PATIENT_BARCODE) %>%
        mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>% distinct
      df_selected <- as.data.frame(df_selected)
      rownames(df_selected) <- df_selected$PATIENT_BARCODE
      df_selected$PATIENT_BARCODE <- NULL
      df_selected <- scale(df_selected,center = T,scale = T)
      df_select_qr <- get_qunatile_rank(df_selected)
      if(is.null(df_select_qr)){
        next
      }else{
        if(ncol(df_select_qr)<3){
          next
        }
      }
      colnames(df_select_qr)[1]<- "PATIENT_BARCODE"
      dft <- merge(df_feature,df_select_qr,by="PATIENT_BARCODE")
      dft <- as.data.frame(dft)
      dft <- dft[dft[,3] %in% c(1,4),]
      dft <- dft[dft[,4] %in% c(1,4),]
      dft <- dft[complete.cases(dft),]
      dft<- dft[,c(2,3,4)]
      cc <- sum(mysign*get_syng_score(dft)$CScore > mysign*myscore)
      if(!is.na(cc)){
        P_Count <- P_Count + cc
      }
    }
  }else{

    for(i in 1:N_iteration){

      df_selected <- as.data.frame(t(df[sample(nrow(df),2,replace = F),]))
      df_selected <- scale(df_selected,center = T,scale = T)
      df_select_qr <- get_qunatile_rank(df_selected)
      if(is.null(df_select_qr)){
        next
      }else{
        if(ncol(df_select_qr)<3){
          next
        }
      }
      colnames(df_select_qr)[1]<- "PATIENT_BARCODE"
      dft <- merge(df_feature,df_select_qr,by="PATIENT_BARCODE")
      dft <- as.data.frame(dft)
      dft <- dft[dft[,3] %in% c(1,4),]
      dft <- dft[dft[,4] %in% c(1,4),]
      dft <- dft[complete.cases(dft),]
      dft<- dft[,c(2,3,4)]
      cc <- sum(mysign*get_syng_score(dft)$CScore > mysign*myscore)
      if(!is.na(cc)){
        P_Count <- P_Count + cc
      }
    }
  }

  P_Count <- P_Count/N_iteration

  return(P_Count)
}


