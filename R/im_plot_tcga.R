#' A correlation function
#'
#' Generates coopretivity boxplots immune checkpoints
#' @param gene1 a single gene name
#' @param gene2 a single gene name
#' @param cohort a single TCGA disease
#' @param Immune_Feautue an immune feature name as listed in im_syng_tcga output.
#' @keywords boxplots
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @export
#' @examples im_plot_tcga(gene1 = "TGFBR1",gene2="CSF1",
#'                   cohort="lihc", Immune_Feauture="Mast.Cells.Activated")
#' im_plot_tcga()


im_plot_tcga<-function(gene1,gene2,cohort,Immune_Feauture){


  cohort <- tolower(cohort)
  results <- list()

  # Read data -----------------------

  cohort_study <- paste0(cohort,"_tcga_pan_can_atlas_2018")
  df <- cBioDataPack(cohort_study,ask=F)@ExperimentList@listData
  df <- df$RNA_Seq_v2_expression_median
  df2 <- df@elementMetadata@listData
  df <- df@assays@data@listData[[1]]
  rownames(df)<- mapvalues(rownames(df),df2$Entrez_Gene_Id,df2$Hugo_Symbol,
    warn_missing = F)

  df_selected <- as.data.frame(t(df[rownames(df) %in% c(gene1,gene2),]))
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }
  if(Immune_Feauture=="EMTscore"){
    df_feature <- pathio::TCGA_EMT
    colnames(df_feature)[1] <- "PATIENT_BARCODE"
  }else{
    if(Immune_Feauture=="Leukocyte_fraction"){
      df_feature <- pathio::TCGA_Leukocyte_fraction
      colnames(df_feature)[1] <- "PATIENT_BARCODE"
    }else{
      df_feature <- pathio::TCGA_IMCell_fraction
      tmp <- which(colnames(df_feature)==Immune_Feauture)
      if(length(tmp)==0){
        stop(Immune_Feauture," Not found. Check ... for list of immune features\n")
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

  #--------------------------------------------

  #construct quantile ranking matrices
  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- Get_qunatile_rank(df_selected)
  colnames(df_select_qr)[1]<- "PATIENT_BARCODE"
  #--------------------------------------------

  #plot
  df_feature <- merge(df_feature,df_select_qr,by="PATIENT_BARCODE")
  df_feature <- as.data.frame(df_feature)
  df_feature <- df_feature[df_feature[,3] %in% c(1,4),]
  df_feature <- df_feature[df_feature[,4] %in% c(1,4),]
  df_feature$state <- as.factor(as.integer((df_feature[,3]*2+df_feature[,4])/3))



  p<-ggplot(df_feature, aes(x=state,y=get(Immune_Feauture),col="grey")) +
    labs(title="", x="", y=Immune_Feauture)+
    geom_boxplot(width=0.5, show.legend = T, outlier.size = 1.0, lwd=0.5, outlier.colour = "red",outlier.shape = 21)+
    geom_jitter(position = position_jitter(width=0.1), cex=1, pch=19, alpha=1)+
    theme( axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10,angle = 90),
      axis.title =element_text(size=10),legend.position ="top")+
    scale_x_discrete(labels =  c(" Both low",paste0(gene2," high\n",gene1," low"), paste0(gene1," high\n",gene2," low"),"Both high"))+
    stat_compare_means(comparisons =
        list( c("1", "2"),c("1", "3"), c("3", "4"), c("2", "4")),
      size=5,method = "wilcox.test")+
    theme(legend.text=element_text(size=10))+theme_bw()

  gg<-ggplot_build(p)
  xx<-gg$data[[1]][c("group","outliers")]
  df_feature2<-merge(df_feature,xx,by.x="state",by.y="group")
  df_feature2$out<-apply(df_feature2,1,function(x) x[Immune_Feauture] %in% x$outliers)

  p<-ggplot(df_feature2, aes(x=state,y=get(Immune_Feauture))) +
    labs(title="", x="", y=Immune_Feauture)+
    geom_boxplot(width=0.5, show.legend = T, lwd=0.5,
      outlier.shape = NA,border="grey")+
    geom_jitter(aes(col=out) ,position = position_jitter(width=0.1),
      cex=1, pch=19, alpha=1)+
    scale_color_manual(values=c("black","red"),guide="none")+
    scale_x_discrete(labels =  c(" Both low",paste0(gene2," high\n",gene1," low"),
      paste0(gene1," high\n",gene2," low"),"Both high"))+
    stat_compare_means(comparisons =
        list( c("1", "2"),c("1", "3"), c("3", "4"), c("2", "4")),
      size=3,method = "wilcox.test")+
    theme( axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10,angle = 90),
      axis.title =element_text(size=10),legend.position ="top")+
    theme(legend.text=element_text(size=10))+
    theme_bw()


  return(p)
}


