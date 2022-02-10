#' Generate a stratified boxplot for immune feature values
#' @importFrom ggplot2 aes coord_trans element_text geom_boxplot geom_jitter ggplot ggplot_build labs  position_jitter scale_color_manual scale_x_discrete theme theme_bw
#' @importFrom ggpubr stat_compare_means
#' @param onco_gene A charachter indicating a single onco_gene ID.
#' @param icp_gene A charachter indicating a single immune checkpoint ID.
#' @param data_expression A numeric matrix or data frame containing gene/protein expressions
#' @param data_feature A numeric matrix or data frame containing a single immune feature.
#' @keywords boxplot, immune features, immune checkpoints
#' @return a stratified boxplot
#' @details
#'
#' Feature data is stratified base on expression quartiles of onco_gene and icp_gene. High/Low categories include samples with expression values in lower/upper quartiles correspondingly. Samples with expression values in middle quartiles are discarded. For details of quartile calculation see get_quantile_rank function.
#'
#'Pvalues are calculated using Wilcoxon test.
#'
#' data_expression is formatted with genes/proteins as rows and samples/patients as columns.
#' For data_expression sample formats see \code{\link[imogimap]{sample_mRNA_data}}.
#'
#' data_feature is formatted with samples/patients as rows and immune feature as single column.
#' For data_feature sample format see \code{\link[imogimap]{sample_Leukocyte_fraction_data}}.
#'
#' @examples im_boxplot(onco_gene = "TGFB1",icp_gene="TNFSF4",
#'                   data_expression =  sample_mRNA_data,
#'                   data_feature = sample_Leukocyte_fraction_data)
#' @export

im_boxplot <- function(onco_gene,icp_gene,data_expression,data_feature){

  state <- NULL
  #Check for onco_gene and immune checkpoint expressions---------
  data_expression <- as.data.frame(data_expression)
  df_selected <- t(data_expression[rownames(data_expression) %in% c(onco_gene,icp_gene),])
  if(nrow(df_selected)==0){
    stop("ERROR: No gene found.Select genes from expression data")
  }

  #Construct quantile ranking matrices------------------
  df_selected <- scale(df_selected,center = T,scale = T)
  df_select_qr <- get_quantile_rank(df_selected)


  #Reformat------------------
  data_feature <- as.data.frame(data_feature)
  data_feature$Tumor_Sample_ID <- rownames(data_feature)
  data_feature <- merge(data_feature,df_select_qr,by="Tumor_Sample_ID")
  data_feature <- data_feature[data_feature[,3] %in% c(1,4),]
  data_feature <- data_feature[data_feature[,4] %in% c(1,4),]
  data_feature$state <- as.factor(as.integer((data_feature[,3]*2+data_feature[,4])/3))


  #--------------------------------------------

  #plot
  Immune_Feature <- colnames(data_feature)[2]
  p=ggplot(data_feature, aes(x=state,y=get(Immune_Feature))) +
    labs(title="", x="", y=Immune_Feature)+
    geom_boxplot(width=0.5, show.legend = T, outlier.size = 1.0, lwd=0.5, outlier.colour = "red",outlier.shape = 21)+
    geom_jitter(position = position_jitter(width=0.1), cex=1, pch=19, alpha=1)+
    theme( axis.text.x=element_text(size=10),
      axis.text.y=element_text(size=10,angle = 90),
      axis.title =element_text(size=10),legend.position ="top")+
    scale_x_discrete(labels =  c(" Both low",paste0(icp_gene," high\n",onco_gene," low"), paste0(onco_gene," high\n",icp_gene," low"),"Both high"))+
    stat_compare_means(comparisons =
        list( c("1", "4"), c("2", "4"), c("3", "4")),
      size=5,method = "wilcox.test")+
    theme(legend.text=element_text(size=10))

  return(p)
}


