#' Generate a stratified boxplot for immune feature values
#' @import dplyr
#' @import cBioPortalData
#' @import ggplot2
#' @import ggpubr
#' @param cotarget A charachter indicating a single cotarget Hugo symbol.
#' @param checkpoint A charachter indicating a single checkpoint Hugo symbol.
#' @param cohort a single TCGA disease
#' @param Immune_Feauture an immune feature name as listed in TCGA_immune_features_list.
#' @param logtrans An optional logical indicating if y axis should be displayed in logarithmic scale. Default is FALSE.
#' @keywords boxplots, immune features, immune checkpoints, cbioportal data
#' @details
#'
#' Feature data is stratified base on expression quartiles of cotarget and checkpoint. High/Low categories include samples with expression values in lower/upper quartiles correspondingly. Samples with expression values in middle quartiles are discarded. For details of quartile calculation see get_quantile_rank function.
#'
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @examples im_boxplot_tcga(cotarget = "BRAF", checkpoint="CD274",
#' cohort="acc", Immune_Feature="Mast.Cells.Activated",logtrans=TRUE)
#' @export

im_boxplot_tcga<-function(cotarget,checkpoint,cohort,Immune_Feature,logtrans){


  cohort <- tolower(cohort)
  results <- list()

  #Read data -----------------------
  cohort_study <- paste0(cohort,"_tcga_pan_can_atlas_2018")
  df <- cBioDataPack(cohort_study,ask=F)@ExperimentList@listData
  df <- df$RNA_Seq_v2_expression_median
  df2 <- df@elementMetadata@listData
  df <- df@assays@data@listData[[1]]
  rownames(df)<- plyr::mapvalues(rownames(df),df2$Entrez_Gene_Id,df2$Hugo_Symbol,
    warn_missing = F)
  df_selected <- as.data.frame(t(df[rownames(df) %in% c(cotarget,checkpoint),]))
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }

  if(Immune_Feature=="EMTscore"){
    df_feature <- get_emt_score(df)
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
          if(grepl("TMB",Immune_Feature)){
            df_feature <- TCGA_TMB[,c("Tumor_Sample_ID",Immune_Feature)]
            colnames(df_feature)[1] <- "PATIENT_BARCODE"
          }else{
            df_feature <- TCGA_IMCell_fraction
            tmp <- which(colnames(df_feature)==Immune_Feature)
            if(length(tmp)==0){
              stop(Immune_Feature," Not found.
                Choose a feature from TCGA_feature_list.\n")
            }else{
              tmpID <- which(colnames(df_feature)=="PATIENT_BARCODE")
              df_feature <- df_feature[,c(tmpID,tmp)]
              df_selected$PATIENT_BARCODE <- substr(rownames(df_selected), 1, 12)
              df_selected <- df_selected %>% group_by(PATIENT_BARCODE) %>%
                mutate(across(cols = everything(),.fns = ~median(.x, na.rm = TRUE))) %>%  distinct
              df_selected <- as.data.frame(df_selected)
              rownames(df_selected) <- df_selected$PATIENT_BARCODE
              df_selected$PATIENT_BARCODE <- NULL
            }
          }
        }
      }
    }
  }


  #construct quantile ranking matrices------------
  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- get_qunatile_rank(df_selected)
  colnames(df_select_qr)[1]<- "PATIENT_BARCODE"


  #Reformat---------------------------------------
  df_feature <- merge(df_feature,df_select_qr,by="PATIENT_BARCODE")
  df_feature <- as.data.frame(df_feature)
  df_feature <- df_feature[df_feature[,3] %in% c(1,4),]
  df_feature <- df_feature[df_feature[,4] %in% c(1,4),]
  df_feature$state <- as.factor(as.integer((df_feature[,3]*2+df_feature[,4])/3))


  #Plot-------------------------------------------
  if(missing(logtrans)){
    logtrans<-FALSE
  }
  if(logtrans){
    transvalue <- "log2"
  }else{
    transvalue <- "identity"
  }

  p<-ggplot(df_feature, aes(x=state,y=get(Immune_Feature),col="grey")) +
    labs(title="", x="", y=Immune_Feature)+
    geom_boxplot(width=0.5, show.legend = T,
      outlier.size = 1.0, lwd=0.5, outlier.colour = "red",outlier.shape = 21)+
    geom_jitter(position = position_jitter(width=0.1), cex=1, pch=19, alpha=1)+
    theme( axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10,angle = 90),
      axis.title =element_text(size=10),legend.position ="top")+
    scale_x_discrete(labels =
        c(" Both low",
          paste0(checkpoint," high\n",cotarget," low"),
          paste0(cotarget," high\n",checkpoint," low"),
          "Both high"))+
    stat_compare_means(comparisons =
        list( c("1", "2"),c("1", "3"), c("3", "4"), c("2", "4")),
      size=5,method = "wilcox.test")+
    theme(legend.text=element_text(size=10))+
    theme_bw()+
    coord_trans(y=transvalue)

  gg<-ggplot_build(p)
  xx<-gg$data[[1]][c("group","outliers")]
  df_feature2<-merge(df_feature,xx,by.x="state",by.y="group")
  df_feature2$out<-apply(df_feature2,1,function(x) x[Immune_Feature] %in% x$outliers)

  p<-ggplot(df_feature2, aes(x=state,y=get(Immune_Feature))) +
    labs(title="", x="", y=Immune_Feature)+
    geom_boxplot(width=0.5, show.legend = T, lwd=0.5,
      outlier.shape = NA,border="grey")+
    geom_jitter(aes(col=out) ,position = position_jitter(width=0.1),
      cex=1, pch=19, alpha=1)+
    scale_color_manual(values=c("black","red"),guide="none")+
    scale_x_discrete(labels =  c(" Both low",paste0(checkpoint," high\n",cotarget," low"),
      paste0(cotarget," high\n",checkpoint," low"),"Both high"))+
    stat_compare_means(comparisons =
        list( c("1", "2"),c("1", "3"), c("3", "4"), c("2", "4")),
      size=3,method = "wilcox.test")+
    theme( axis.text.x=element_text(size=15),
      axis.text.y=element_text(size=10,angle = 90),
      axis.title =element_text(size=10),legend.position ="top")+
    theme(legend.text=element_text(size=10))+
    theme_bw()+
    coord_trans(y=transvalue)


  return(p)
}


