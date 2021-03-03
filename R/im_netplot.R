#' Creates and plot an igraph network from dataframe
#' @importFrom igraph graph.edgelist E E<- set_edge_attr V V<- degree
#' @param df A dataframe as outputted by im_syng_tcga or im_syng
#' @param Immune_Feature a character string indicating name of an immune feature.
#' @param icp_gene An optional character vector of immune checkpoints for color coding vertex. default is icp_gene_list
#' @param cohort An optional character string indicating a single TCGA disease.
#' @param cutoff A numeric indicating cut-off for absolute values of synergy scores.
#' @param seed A single value, interpreted as an integer, or NULL (see set.seed for details).
#' @return A plotted igraph object showing network of synergistic (red) and antagonistic (blue) associations.
#' @examples
#' df <- im_syng_tcga(onco_gene = icp_gene_list,icp_gene = icp_gene_list,cohort = "acc")
#' im_netplot(df, Immune_Feature = "EMTscore",cutoff = 5, seed = 1234)
#' @details
#'
#' im_netplot constructs and plots network from synergy score data frame. Immune-checkpoints/onco-genes are depicted as black/white vertices, and positive/negative synergistic interactions are depicted as red/blue edges. Thickness of an edge is determined by the absolute value of the score, and the size of each vertex is determined by its degree. Edges with absolute values lower that the cut-off will be removed from the graph.
#'
#' seed is an integer that is the starting point from which random numbers are generated in igraph. Seed value is used for the reproducibility of plot layout. different seed numbers will generate different
#' @export
#'

im_netplot<- function(df, icp_gene, cohort, Immune_Feature, cutoff,seed) {

  if(missing(seed)){
    seed<- 1
  }
  if(missing(cutoff)){
    cutoff <- 0.05
  }
  if(missing(icp_gene)){
    icp_gene <- icp_gene_list
  }
  if(!missing(cohort)){
   df <- df[df$Disease==cohort,]
  }

  df <- df[complete.cases(df),]

  df <- df[df$Immune_feature== Immune_Feature,]
  df <- df[abs(df$Synergy_score) > cutoff,]
  if(nrow(df)==0){
    stop("ERROR: There is no data to show. Try decreasing cutoff.")
  }
  df$Synergy_type <- ifelse(sign(df$Synergy_score)==1, "Synergistic", "antagonistic")
  df$Synergy_score <- abs(df$Synergy_score)

  df <- df[,c("agent1","Synergy_type","Synergy_score","agent2")]
  colnames(df)<- c("gene_A", "Synergy_type","Synergy_score", "gene_B")

  df <- df[complete.cases(df),]
  df <- df[!duplicated(t(apply(df, 1, sort))),]

  sif <- df[, c("gene_A", "gene_B")]
  g <- igraph::graph.edgelist(as.matrix(sif), directed=FALSE)

  g <- set_edge_attr(g, "Synergy_type", index=E(g), df[, "Synergy_type"])
  g <- set_edge_attr(g, "Synergy_score", index=E(g), df[, "Synergy_score"])

  cg <- igraph::edge.betweenness.community(g)
  ew <-  20*E(g)$Synergy_score/max( E(g)$Synergy_score)+1
  ecol <- ifelse(E(g)$Synergy_type=="Synergistic","red","blue")
  ecol <- alpha(ecol,0.4)
  V(g)$color <- ifelse(names(V(g)) %in% icp_gene,"grey","white")
  vsize <- 4*degree(g)/max(degree(g))+1
  dq <- quantile(degree(g),probs = c(0,0.3,0.6,1))
  vcex <- ifelse(degree(g)<dq[2],2,ifelse(degree(g)<dq[3],3,4))
  set.seed(seed)
  par(mar=rep(0,4))

  p<- plot(g,
    vertex.size = vsize,
    vertex.label.family = "Helvetica",
    vertex.label.font = 2,
    vertex.label.cex=vcex,
    vertex.label.dist=.5,
    edge.color = ecol,
    edge.width = ew,
    mark.groups = by(seq_along(cg$membership), cg$membership, invisible))

  return(p)
}
