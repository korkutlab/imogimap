#' Convert SIF to to igraph network
#' @importFrom igraph graph.edgelist E E<- set_edge_attr
#' @param df A dataframe as outputed by im_syng_tcga or im_syng
#' @param Immune_Feature a charachter string indicating name of animmune feature.
#' @param checkpoint An optional character vector of immune checkpoints for color coding vertex. default is icp_gene_list
#' @param cohort An optional charachter string indicating a single TCGA disease.
#' @param cutoff A numeric indicating cut-off for absolute values of synergy scores.
#' @param seed A single value, interpreted as an integer, or NULL (see set.seed for details).
#' @return A plotted igraph object showing network of synergistic (red) and antagonistic (blue) associations.
#' @examples df <- im_syng_tcga(cotarget = icp_gene_list,checkpoint = icp_gene_list,cohort = "brca")
#' im_netplot(df, Immune_Feature = "EMTscore",cutoff = 5, seed = 1234)
#' @details
#'
#' im_netplot egenrate and plots network from synergy-score data frame. vertexes are gene names in dataframe and edges are synergistic (red) and antagonistic (blue) associations. Edges with absolute values lower that the cut-off wil be removed from the graph.
#'
#' seed is an integer that is the starting point from which random numbers are generated in igraph. Seed value is used for the reproducibility of plot layout. different seed numbers will generate different
#' @export
#'

im_netplot<- function(df, checkpoint, cohort, Immune_Feature, cutoff,seed) {

  if(missing(cutoff)){
    cutoff <- 2
  }
  if(missing(checkpoint)){
    checkpoint <- icp_gene_list
  }
  if(!missing(cohort)){
   df <- df[df$Disease==cohort,]
  }

  df <- df[complete.cases(df),]

  df <- df[df$Immune_feature== Immune_Feature,]
  df <- df[abs(df$Synergy_score) > cutoff,]

  df$Synergy_type <- ifelse(sign(df$Synergy_score)==1, "Synergistic", "antagonistic")
  df$Synergy_score <- abs(df$Synergy_score)

  df <- df[,c("Co_target","Synergy_type","Synergy_score","Immune_checkpoint")]
  colnames(df)<- c("gene_A", "Synergy_type","Synergy_score", "gene_B")

  df <- df[complete.cases(df),]
  df <- df[!duplicated(t(apply(df, 1, sort))),]

  sif <- df[, c("gene_A", "gene_B")]
  g <- igraph::graph.edgelist(as.matrix(sif), directed=FALSE)

  g <- set_edge_attr(g, "Synergy_type", index=E(g), df[, "Synergy_type"])
  g <- set_edge_attr(g, "Synergy_score", index=E(g), df[, "Synergy_score"])

  cg <- igraph::edge.betweenness.community(g)
  ew <-  4*E(g)$Synergy_score/max( E(g)$Synergy_score)+1
  ecol <- ifelse(E(g)$Synergy_type=="Synergistic","red","blue")
  ecol <- alpha(ecol,0.7)
  V(g)$color <- ifelse(names(V(g)) %in% checkpoint,"tan2","plum2")
  vsize <- 4*degree(g)/max(degree(g))+1
  dq <- quantile(degree(g),probs = c(0,0.3,0.6,1))
  vcex <- ifelse(degree(g)<dq[2],0.7,ifelse(degree(g)<dq[3],1,1.3))
  set.seed(seed)
  par(mar=rep(0,4))
  p <- plot(g,
    vertex.size = vsize,
    vertex.label.family = "Helvetica",
    vertex.label.cex=vcex,
    vertex.label.dist=.1,
    edge.color = ecol,
    edge.width = ew,
    mark.groups = by(seq_along(cg$membership), cg$membership, invisible))

  return(p)
}
