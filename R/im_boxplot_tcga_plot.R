#' Plot data from im_bloxplot_tcga 
#' 
#' @importFrom ggplot2 aes coord_trans element_text geom_boxplot geom_jitter ggplot ggplot_build labs  position_jitter scale_color_manual scale_x_discrete theme theme_bw
#' @importFrom ggpubr stat_compare_means
#' 
#' @param obj a list obj produced by im_bloxplot_tcga
#' 
#' @return a ggplot object
#' 
#' @examples 
#' # FIXME
#' 
#' 
#' @export
im_boxplot_tcga_plot <- function(obj) {
  logtrans <- obj[["logtrans"]]
  df_feature <- obj[["df_feature"]]
  Immune_phenotype <- obj[["Immune_phenotype"]]
  mylabels <- obj[["mylabels"]]
  
  if(missing(logtrans)){
    logtrans<-FALSE
  }
  if(logtrans){
    transvalue <- "log2"
    df_feature[df_feature[,2]==0,2]<-NA
    plot_min <- min(df_feature[,2],na.rm=TRUE)
    df_feature[is.na(df_feature[,2]),2] <- plot_min
  }else{
    transvalue <- "identity"
  }
  
  p<-ggplot(df_feature, aes(x=state,y=get(Immune_phenotype),col="grey")) +
    labs(title="", x="", y=Immune_phenotype)+
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
  df_feature2$out<-apply(df_feature2,1,function(x) x[Immune_phenotype] %in% x$outliers)
  
  p <- ggplot(df_feature2, aes(x=state,y=get(Immune_phenotype))) +
    labs(title="", x="", y=Immune_phenotype)+
    geom_boxplot(width=0.5, show.legend = T, lwd=1,
                 outlier.shape = NA,border="grey")+
    geom_jitter(aes(col=out) ,position = position_jitter(width=0.1),
                cex=2, pch=19, alpha=1)+
    scale_color_manual(values=c("black","red"),guide="none")+
    scale_x_discrete(labels =  mylabels)+
    theme_bw()+
    stat_compare_means(comparisons =
                         list( c("1", "4"),c("2", "4"), c("3", "4")),
                       size=10,method = "wilcox.test",
                       bracket.size=1,vjust=1.5)+
    theme( axis.text.x=element_text(size=25),
           axis.text.y=element_text(size=25,angle = 90),
           axis.title =element_text(size=25),legend.position ="top")+
    theme(legend.text=element_text(size=25))+
    coord_trans(y=transvalue)
  
  return(p)
}
