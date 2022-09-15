#' Generates a stratified boxplot for immune feature values based of a two genes.
#'
#' @importFrom dplyr across mutate group_by everything
#' @import curatedTCGAData
#' @importFrom stats median
#'
#' @param onco_gene A character indicating a single onco_gene Hugo symbol.
#' @param icp_gene A character indicating a single immune checkpoint Hugo symbol.
#' @param cohort a single TCGA disease
#' @param sample_list An optional character vector of TCGA samples barcodes indicating a subset of samples within a cohort. All barcodes in sample_list must be 15 character long and belong to the same cohort.
#' @param Immune_phenotype an immune phenotype name as listed in TCGA_Immune_phenotypes_list.
#' @param logtrans An optional logical indicating if y axis should be displayed in logarithmic scale. Default is FALSE.
#' @keywords boxplots, immune features, immune checkpoints, TCGA
#' @details
#'
#' Feature data is stratified based on expression quartiles of onco_gene and icp_gene. High/Low categories include samples with expression values in lower/upper quartiles correspondingly. Samples with expression values in middle quartiles are discarded. For details of quartile calculation see \code{\link[imogimap]{get_quantile_rank}} function.
#'
#' Pvalues are calculated using Wilcoxon test.
#'
#' @return a list of multiple dataframes of correlation coefficients and p.values FIXME (EVEN BEFORE MY EDITS IT WAS A GGPLOT NOT THIS)
#' @examples im_boxplot_tcga(onco_gene = "BRAF", icp_gene="CD274",
#' cohort="acc", Immune_phenotype="Mast.Cells.Activated",logtrans=TRUE)
#' @seealso get_quantile_rank
#' @export
im_boxplot_tcga<-function(onco_gene,icp_gene,cohort,Immune_phenotype,sample_list,logtrans=FALSE){

  PATIENT_BARCODE <- NULL
  state <- NULL
  out <- NULL
  cohort <- tolower(cohort)
  results <- list()

  #Read data -----------------------
  df <-curatedTCGAData::curatedTCGAData( diseaseCode = cohort,version = "1.1.38",
                                         assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]
  df <- df@assays@data@listData[[1]]
  colnames(df)<-  substr(colnames(df), 1, 15)

  if(!missing(sample_list)){
    df<-df[,sample_list ]
    if(ncol(df)==0){
      stop("ERROR: barcodes not found.")
    }
  }
  icp_gene <- ifelse(icp_gene=="VSIR", "C10orf54",
                     ifelse(icp_gene=="NCR3LG1","DKFZp686O24166",icp_gene))
  onco_gene <- ifelse(onco_gene=="VSIR", "C10orf54",
                      ifelse(onco_gene=="NCR3LG1","DKFZp686O24166",onco_gene))

  df_selected <- as.data.frame(t(df[rownames(df) %in% c(onco_gene,icp_gene),]))
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }
  df_selected <- df_selected[,c(onco_gene,icp_gene)]
  if(Immune_phenotype=="EMTscore"){
    df_feature <- get_emt_score(df)
    colnames(df_feature)[1] <- "PATIENT_BARCODE"
  }else{
    if(Immune_phenotype=="Leukocyte_fraction"){
      df_feature <- TCGA_Leukocyte_fraction
      colnames(df_feature)[1] <- "PATIENT_BARCODE"
    }else{
      if(Immune_phenotype=="AGscore"){
        df_feature <- get_angio_score(df)
        colnames(df_feature)[1] <- "PATIENT_BARCODE"
      }else{
        if(Immune_phenotype=="IFNGscore"){
          df_feature <- get_ifng_score(df)
          colnames(df_feature)[1] <- "PATIENT_BARCODE"
        }else{
          if(Immune_phenotype=="TCIscore"){
            df_feature <- get_TCI_score(df)
            colnames(df_feature)[1] <- "PATIENT_BARCODE"
          }else{

            if(grepl("TMB",Immune_phenotype)){
              df_feature <- TCGA_TMB[,c("Tumor_Sample_ID",Immune_phenotype)]
              colnames(df_feature)[1] <- "PATIENT_BARCODE"
            }else{
              df_feature <- TCGA_IMCell_fraction
              tmp <- which(colnames(df_feature)==Immune_phenotype)
              if(length(tmp)==0){
                stop(Immune_phenotype," Not found.
                Choose a feature from TCGA_feature_list.\n")
              }else{
                tmpID <- which(colnames(df_feature)=="PATIENT_BARCODE")
                df_feature <- df_feature[,c(tmpID,tmp)]
                df_selected$PATIENT_BARCODE <- substr(rownames(df_selected), 1, 12)
                df_selected <- df_selected %>% group_by(PATIENT_BARCODE) %>%
                  mutate(across(.cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>%  distinct
                df_selected <- as.data.frame(df_selected)
                rownames(df_selected) <- df_selected$PATIENT_BARCODE
                df_selected$PATIENT_BARCODE <- NULL
              }
            }
          }
        }
      }
    }
  }

  #construct quantile ranking matrices------------
  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- get_quantile_rank(df_selected)
  colnames(df_select_qr)[1]<- "PATIENT_BARCODE"


  #Reformat---------------------------------------
  df_feature <- merge(df_feature,df_select_qr,by="PATIENT_BARCODE")
  df_feature <- as.data.frame(df_feature)
  df_feature <- df_feature[df_feature[,3] %in% c(1,4),]
  df_feature <- df_feature[df_feature[,4] %in% c(1,4),]

  #Find if expression or inhibition of genes positively impact feature
  my_score <- find_a_synergy(df_feature[,-1],method = "max",ndatamin=8)
  if(is.na(my_score$Synergy_score)){
    effect_onco <- "High"
    effect_icp <- "High"
  }else{
    effect_onco <- my_score$agent1_expression
    effect_icp <- my_score$agent2_expression
  }

  #Flip labels if inhibition of genes positively impact feature
  if(effect_onco=="Low"){
    df_feature[df_feature[,3]==4,3] <- 0
    df_feature[df_feature[,3]==1,3] <- 4
    df_feature[df_feature[,3]==0,3] <- 1

  }
  if(effect_icp=="Low"){
    df_feature[df_feature[,4]==4,4] <- 0
    df_feature[df_feature[,4]==1,4] <- 4
    df_feature[df_feature[,4]==0,4] <- 1
  }

  #Define labels
  df_feature$state <- as.factor(as.integer((df_feature[,3]*2+df_feature[,4])/3))

  if(effect_onco=="High" && effect_icp=="High"){
    mylabels <-   c("Both low",
                    paste0(icp_gene," high\n",onco_gene," low"),
                    paste0(onco_gene," high\n",icp_gene," low"),
                    "Both high")
  }else{
    if(effect_onco=="Low" && effect_icp=="Low"){
      mylabels <- c("Both high",
                    paste0(onco_gene," high\n",icp_gene," low"),
                    paste0(icp_gene," high\n",onco_gene," low"),
                    "Both low")
    }else{
      if(effect_onco=="Low" && effect_icp=="High"){
        mylabels <- c( paste0(icp_gene," low\n",onco_gene," high"),
                       "Both high",
                       "Both low",
                       paste0(onco_gene," low\n",icp_gene," high"))
      }else{
        mylabels <- c(paste0(onco_gene," low\n",icp_gene," high"),
                      "Both low",
                      "Both high",
                      paste0(icp_gene," low\n",onco_gene," high"))

      }
    }
  }

  results <- list(df_feature=df_feature, logtrans=logtrans, Immune_phenotype=Immune_phenotype, mylabels=mylabels)
  return(results)
}


