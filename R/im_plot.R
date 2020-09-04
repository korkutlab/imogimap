#' A correlation function
#'
#' Generates coopretivity boxplots immune checkpoints
#' @param gene1 a single gene name
#' @param gene2 a single gene name
#' @param Immune_Feautue an immune feature name as listed in im_syng output.
#' @param df_mrna a formated mRNA data frame
#' @param df_ict an optional formated immune cell type fractions data frame
#' @param df_lf an optional formated Leukocyte fraction data frame
#' @keywords boxplots
#' @return a list of multiple dataframes of correlation coefficients and p.values
#' @export
#' @examples im_plot(gene1 = "TGFB1",gene2="TNFSF4",
#'                   Immune_Feautue="EMT_score",
#'                   df_mrna =  sample_mRNA_data,
#'                   df_lf = sample_Leukocyte_fraction_data,
#'                   df_ict = sample_immune_cell_fraction_data)
#' im_plot()


im_plot<-function(gene1,gene2,Immune_Feautue,df_mrna,df_ict,df_lf){


  #read data
  mydata <- as.data.frame(df_mrna)
  df_selected <- t(mydata[rownames(mydata) %in% c(gene1,gene2),])
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found. Select a gene name from your mRNA data")
  }
  if(Immune_Feautue=="EMTscore"){
    df_EMT <-  Get_EMTscore(mydata)
    if(nrow(df_EMT)==0){
      stop("No EMT signature gene found.\n")
    }else{
      df_feature <- df_EMT
    }
  }else{
    if(Immune_Feautue=="Leukocyte_fraction"){
      if(missing(df_lf)){
        stop("No Leukocyte fraction provided.\n")
      }else{
        df_lf2 <- as.data.frame(df_lf)
        colnames(df_lf2)<-"Leukocyte_fraction"
        df_lf2$Tumor_Sample_ID <- rownames(df_lf2)
        df_feature <- df_lf2
      }
    }else{
      if(missing(df_ict)){
        stop("No immune cell fractions provided.\n")
      }else{
        df_ict2 <- as.data.frame(df_ict)
        df_ict2$Tumor_Sample_ID <- rownames(df_ict2)
        tmp <- which(colnames(df_ict2)==Immune_Feautue)
        if(length(tmp)==0){
          stop(Immune_Feautue," Not found.\n")
        }else{
          tmpID <- which(colnames(df_ict2)=="Tumor_Sample_ID")
          df_feature <- df_ict2[,c(tmpID,tmp)]
        }
      }
    }
  }
  #--------------------------------------------

  #construct quantile ranking matrices
  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- Get_qunatile_rank(df_selected)
  #--------------------------------------------

  #plot
  df_feature <- merge(df_feature,df_select_qr,by="Tumor_Sample_ID")
  df_feature <- df_feature[df_feature[,3] %in% c(1,4),]
  df_feature <- df_feature[df_feature[,4] %in% c(1,4),]
  df_feature$state <- as.factor(as.integer((df_feature[,3]*2+df_feature[,4])/3))


  p=ggplot(df_feature, aes(x=state,y=get(Immune_Feautue))) +
    labs(title="", x="", y=Immune_Feautue)+
    geom_boxplot(width=0.5, show.legend = T, outlier.size = 1.0, lwd=0.5, outlier.colour = "red",outlier.shape = 21)+
    geom_jitter(position = position_jitter(width=0.1), cex=1, pch=19, alpha=1)+
    theme( axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10,angle = 90),
      axis.title =element_text(size=10),legend.position ="top")+
    scale_x_discrete(labels =  c(" Both low",paste0(gene2," high\n",gene1," low"), paste0(gene1," high\n",gene2," low"),"Both high"))+
    stat_compare_means(comparisons =
        list( c("1", "2"),c("1", "3"), c("3", "4"), c("2", "4")),
      size=5,method = "wilcox.test")+
    theme(legend.text=element_text(size=10))

  return(p)
}


