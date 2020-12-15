#' Generates SIF formatted data-frame 
#'
#' @param df a data frame of synergy scores as outputed by im_syng
#' @param cohort a charachter string indicating a single cohort
#' @param Immune_Feature a charachter string indicating a single immune feature
#' @param cutoff a numeric value indicating the minimum cutoff value for absolute values of  synergy score. Default is 1 
#' @keywords SIF 
#' @details 
#' get_sif generates a data-frame in SIF format. (See http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats/)   
#' @return a binary SIF as a data.frame with three columns:
#'   "Cotarget", "Score", "Checkpoint"
#' @examples df <- im_syng_tcga(cotarget = c("BRAF","AKT1"),cohort = c("acc","gbm"))
#'           get_sif(df=df,cohort="acc", Immune_Feature="B.Cells.Naive"   )
#' @export

get_sif<- function(df,cohort, Immune_Feature,cutoff) {
  
  if(missing(cutoff)){
    cutoff <- 1
  }
  df <- df[complete.cases(df),]
  df <- df[df$Disease==cohort,]
  df <- df[df$Immune_feature== Immune_Feature,]
  df <- df[abs(df$Synergy_score) > cutoff,]
  df$INTERACTION_TYPE <- ifelse(sign(df$Synergy_score)==1, "Synergistic", "antagonistic")
  df <- df[,c("Co_target","INTERACTION_TYPE","Immune_checkpoint")]
  colnames(df)<- c("gene_A", "INTERACTION_TYPE", "gene_B")
  df <- df[complete.cases(df),]
  df <- df[!duplicated(t(apply(df, 1, sort))),]
  return(df)
}
