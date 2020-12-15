#' Convert SIF to to igraph network
#' @import igraph
#' @param sif a binary SIF as a data.frame with three columns:
#'   "gene_A", "INTERACTION_TYPE", "gene_B"
#' @param directed a boolean weather the returned graph should be directed (DEFAULT: TRUE)
#' @return a directed igraph network with interaction types
#' @examples
#' @details Users are likely to run into issues if the input SIF has factor levels
#' @export
#'

im_netplot<- function(df,cohort, Immune_Feature,cutoff) {

  if(missing(cutoff)){
    cutoff <- 1
  }
  df <- df[complete.cases(df),]
  df <- df[df$Disease==cohort,]
  df <- df[df$Immune_feature== Immune_Feature,]
  df <- df[abs(df$Synergy_score) > cutoff,]

  df$Synergy_type <- ifelse(sign(df$Synergy_score)==1, "Synergistic", "antagonistic")
  df$Synergy_score <- abs(df$Synergy_score)

  df <- df[,c("Co_target","Synergy_type","Synergy_score","Immune_checkpoint")]
  colnames(df)<- c("gene_A", "Synergy_type","Synergy_score", "gene_B")

  df <- df[complete.cases(df),]
  df <- df[!duplicated(t(apply(df, 1, sort))),]

  sif <- df[, c("gene_A", "gene_B")]
  g <- graph.edgelist(as.matrix(sif), directed=FALSE)

  g <- set_edge_attr(g, "Synergy_type", index=E(g), df[, "Synergy_type"])
  g <- set_edge_attr(g, "Synergy_score", index=E(g), df[, "Synergy_score"])

  ew <- E(g)$Synergy_score
  ecol <- ifelse(E(g)$Synergy_type=="Synergistic","red","blue")
  cg <- igraph::edge.betweenness.community(g)

  p <- plot(g, vertex.size = degree(g),edge.color=ecol,edge.width=ew,
    mark.groups = by(seq_along(cg$membership), cg$membership, invisible))

  return(p)
}
