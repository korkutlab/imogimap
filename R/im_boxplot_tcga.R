#' Generates a stratified boxplot for immune feature values based of a two genes.
#' @import dplyr
#' @import curatedTCGAData
#' @import ggplot2
#' @import ggpubr
#' @param onco_gene A character indicating a single onco_gene Hugo symbol.
#' @param icp_gene A character indicating a single immune checkpoint Hugo symbol.
#' @param cohort a single TCGA disease
#' @param sample_list An optional character vector of TCGA samples barcodes indicating a subset of samples within a cohort. All barcodes in sample_list must be 15 character long and belong to the same cohort.
#' @param Immune_Feauture an immune feature name as listed in TCGA_immune_features_list.
#' @param logtrans An optional logical indicating if y axis should be displayed in logarithmic scale. Default is FALSE.
#' @keywords boxplots, immune features, immune checkpoints, TCGA
#' @details
#'
#' Feature data is stratified based on expression quartiles of onco_gene and icp_gene. High/Low categories include samples with expression values in lower/upper quartiles correspondingly. Samples with expression values in middle quartiles are discarded. For details of quartile calculation see \code{\link[imogene]{get_quantile_rank}} function.
#'
#' Pvalues are calculated using Wilcoxon test.
#'
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @examples im_boxplot_tcga(onco_gene = "BRAF", icp_gene="CD274",
#' cohort="acc", Immune_Feature="Mast.Cells.Activated",logtrans=TRUE)
#' @seealso get_quantile_rank
#' @export

im_boxplot_tcga<-function(onco_gene,icp_gene,cohort,Immune_Feature,sample_list,logtrans){

  cohort <- tolower(cohort)
  results <- list()

  #Read data -----------------------
  df <-curatedTCGAData::curatedTCGAData( diseaseCode = cohort,
                                         assays = c("RNASeq2GeneNorm"), dry.run = F)@ExperimentList@listData[[1]]
  df <- df@assays$data@listData[[1]]
  colnames(df)<-  substr(colnames(df), 1, 15)

  if(!missing(sample_list)){
    df<-df[,sample_list ]
    if(ncol(df==0)){
      stop("ERROR: barcodes not found.")
    }
  }
  df_selected <- as.data.frame(t(df[rownames(df) %in% c(onco_gene,icp_gene),]))
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
  my_score <- find_a_synergy(df_feature[,-1],method = "max")
  if(is.na(my_score$Synergy_score)){
    effect_onco <- "Expressed"
    effect_icp <- "Expressed"
  }else{
    effect_onco <- my_score$agent1_expression
    effect_icp <- my_score$agent2_expression
  }
  #Flip labels if inhibition of genes positively impact feature
  if(effect_onco=="Inhibited"){
    df_feature[df_feature[,3]==4,3] <- 0
    df_feature[df_feature[,3]==1,3] <- 4
    df_feature[df_feature[,3]==0,3] <- 1

  }
  if(effect_icp=="Inhibited"){
    df_feature[df_feature[,4]==4,4] <- 0
    df_feature[df_feature[,4]==1,4] <- 4
    df_feature[df_feature[,4]==0,4] <- 1
  }

  #Define labels
  df_feature$state <- as.factor(as.integer((df_feature[,3]*2+df_feature[,4])/3))

  if(effect_onco=="Expressed" && effect_icp=="Expressed"){
    mylabels <-   c("Both low",
                    paste0(icp_gene," high\n",onco_gene," low"),
                    paste0(onco_gene," high\n",icp_gene," low"),
                    "Both high")
  }else{
    if(effect_onco=="Inhibited" && effect_icp=="Inhibited"){
      mylabels <- c("Both high",
                    paste0(onco_gene," high\n",icp_gene," low"),
                    paste0(icp_gene," high\n",onco_gene," low"),
                    "Both low")
    }else{
      if(effect_onco=="Inhibited" && effect_icp=="Expressed"){
        mylabels <- c( paste0(icp_gene," low\n",onco_gene," high"),
                       "Both low",
                       "Both high",
                       paste0(onco_gene," low\n",icp_gene," high"))
      }else{
        mylabels <- c(paste0(onco_gene," low\n",icp_gene," high"),
                      "Both high",
                      "Both low",
                      paste0(icp_gene," low\n",onco_gene," high"))

      }
    }
  }
  #Plot-------------------------------------------
  if(missing(logtrans)){
    logtrans<-FALSE
  }
  if(logtrans){
    transvalue <- "log2"
    df_feature[df_feature[,2]==0,2]<-NA
    plot_min <- min(df_feature[,2],na.rm=T)
    df_feature[is.na(df_feature[,2]),2] <- plot_min
  }else{
    transvalue <- "identity"
  }

  p<-ggplot(df_feature, aes(x=state,y=get(Immune_Feature),col="grey")) +
    labs(title="", x="", y=Immune_Feature)+
    geom_boxplot(width=0.5, show.legend = T,
                 outlier.size = 1.0, lwd=0.5,
                 outlier.colour = "red",outlier.shape = 21)+
    geom_jitter(position = position_jitter(width=0.1), cex=1, pch=19, alpha=1)+
    theme( axis.text.x=element_text(size=10),
           axis.text.y=element_text(size=10,angle = 90),
           axis.title =element_text(size=10),legend.position ="top")+
    scale_x_discrete(labels = mylabels)+
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


  p <- ggplot(df_feature2, aes(x=state,y=get(Immune_Feature))) +
    labs(title="", x="", y=Immune_Feature)+
    geom_boxplot(width=0.5, show.legend = T, lwd=1,
                 outlier.shape = NA,border="grey")+
    geom_jitter(aes(col=out) ,position = position_jitter(width=0.1),
                cex=2, pch=19, alpha=1)+
    scale_color_manual(values=c("black","red"),guide="none")+
    scale_x_discrete(labels =  mylabels)+
    theme_bw()+
    stat_compare_means(comparisons =
                         list( c("1", "2"),c("1", "3"), c("3", "4"), c("2", "4")),
                       size=10,method = "wilcox.test",
                       bracket.size=1,vjust=1.5)+
    theme( axis.text.x=element_text(size=25),
           axis.text.y=element_text(size=25,angle = 90),
           axis.title =element_text(size=25),legend.position ="top")+
    theme(legend.text=element_text(size=25))+
    coord_trans(y=transvalue)
  return(p)
}


