#' Creates and plot an igraph network from dataframe
#'
#' @importFrom igraph graph.edgelist E E<- set_edge_attr V V<- degree ecount
#' @importFrom  ggplot2 alpha
#'
#' @param df A dataframe as outputted by \code{\link{im_syng_tcga}} or \code{\link[imogimap]{im_syng}}
#' @param Immune_phenotype a character string indicating name of an immune feature.
#' @param icp_gene An optional character vector of immune checkpoints for color coding vertex. default is \code{\link[imogimap]{icp_gene_list}}
#' @param cohort An optional character string indicating a single TCGA disease.
#' @param cutoff A numeric indicating quantile cut-off for absolute values of synergy scores.
#' @param seed A single value, interpreted as an integer; Default: 1 (see set.seed for details).
#'
#' @return A plotted igraph object showing network of synergistic (red) and antagonistic (blue) associations.
#'
#' @examples
#' df <- im_syng_tcga(onco_gene = icp_gene_list[1:3], icp_gene = icp_gene_list[1:3], cohort = "acc")
#' im_netplot(df, Immune_phenotype = "EMTscore")
#'
#' @details
#'
#' im_netplot uses \pkg{igraph} to construct and plot network from synergy score data frame. Dark red (Dark blue) vertices identify over-expression (low-expression) of Immune-checkpoint genes and orange (blue) vertices identify over-expression (low-expression) of oncogenes. up-regulations (down-regulations) of the immune phenotype are depicted as red (blue) edges. Thickness of an edge is determined by the absolute value of the score, and the size of each vertex is determined by its degree. Edges with absolute values lower that the cut-off will be removed from the graph.
#'
#' cutoff=0 does not remove any edges. cutoff=1 removes all edges. Default is 0.95, which keeps scores with their absolute value  higher than 0.95 quantile.
#' seed is an integer that is the starting point from which random numbers are generated in igraph. Seed value is used for the reproducibility of plot layout. different seed numbers will generate different
#' @seealso [igraph]
#'
#' @export
#'
im_netplot<- function(df, icp_gene, cohort, Immune_phenotype, cutoff, seed=1) {

  if(missing(icp_gene)){
    icp_gene <- icp_gene_list
  }

  if(!missing(cohort)){
    df <- df[df$Disease==cohort,]
  }

  df <- df[complete.cases(df$Synergy_score),]
  df <- df[df$IAP== Immune_phenotype,]


  if(missing(cutoff)){
    cutoff <- as.numeric(quantile(abs(df$Synergy_score),probs=c(0.95)))
  }else{
    cutoff <- as.numeric(quantile(abs(df$Synergy_score),probs=c(cutoff)))
  }

  df <- df[abs(df$Synergy_score) >= cutoff,]
  if(nrow(df)==0){
    stop("ERROR: There is no data to show. Try decreasing cutoff.")
  }

  df$Synergy_sign <- ifelse(sign(df$Synergy_score)==1, "+", "-")
  df$Synergy_score <- abs(df$Synergy_score)
  df$Gene1_expression <- ifelse(df$Gene1_expression=="High","+","-")
  df$Gene2_expression <- ifelse(df$Gene2_expression=="High","+","-")
  df$gene1 <- paste0(df$Gene1,df$Gene1_expression)
  df$gene2 <- paste0(df$Gene2,df$Gene2_expression)


  df <- df[,c("gene1","Synergy_sign","Synergy_score","ligand_receptor_interaction","gene2")]
  colnames(df)<- c("gene_A", "Synergy_sign","Synergy_score","interaction" ,"gene_B")
  df <- df[!duplicated(t(apply(df, 1, sort))),]

  df_receptor_ligand <-df[df$interaction==TRUE,]
  if(nrow(df_receptor_ligand)>0){
    df_receptor_ligand$interaction<-FALSE
    df <- rbind(df,df_receptor_ligand)
  }
  sif <- df[, c("gene_A", "gene_B")]
  g <- igraph::graph.edgelist(as.matrix(sif), directed=FALSE)
  g <- igraph::set_edge_attr(g, "Synergy_sign", index=E(g), df[, "Synergy_sign"])
  g <- igraph::set_edge_attr(g, "Synergy_score", index=E(g), df[, "Synergy_score"])
  g <- igraph::set_edge_attr(g, "weight", index=E(g), abs(df[, "Synergy_score"]))
  g <- igraph::set_edge_attr(g, "group", index=E(g), df[, "interaction"])

  E(g)$color <- ifelse(E(g)$Synergy_sign=="+","red","blue")
  E(g)$color <- alpha(E(g)$color,0.5)
  E(g)$color <- ifelse(E(g)$group==TRUE,"black",E(g)$color)

  V(g)$color <- ifelse(substr(names(V(g)),1,nchar(names(V(g)))-1) %in% icp_gene,"grey",
                       ifelse(substr(names(V(g)),nchar(names(V(g))),nchar(names(V(g))))=="+","darkorange1", "skyblue"))
  V(g)$color <- ifelse(substr(names(V(g)),1,nchar(names(V(g)))-1) %in% icp_gene,ifelse(substr(names(V(g)),nchar(names(V(g))),nchar(names(V(g))))=="+","red4","skyblue4"), V(g)$color)


  ming <- min(E(g)$Synergy_score)
  E(g)$width <- 5*(E(g)$Synergy_score-ming)/(max(E(g)$Synergy_score)-ming)+2
  E(g)$width <- ifelse(E(g)$group==TRUE,5,E(g)$width)
  gcurve <- rep(-0.5, length = igraph::ecount(g))
  gcurve[E(g)$group==TRUE]<-0.5
  #par(omi=c(0,0,0,1))
  set.seed(seed)
  p <- plot(g,
            vertex.size = 8,
            vertex.label = substr(names(V(g)),1,nchar(names(V(g)))-1),
            vertex.label.family = "Helvetica",
            vertex.frame.color = "black",
            rescale=TRUE,add=FALSE,
            vertex.label.font = 1,
            vertex.label.color = "black",
            vertex.label.cex=0.80,
            vertex.label.dist=1,
            #edge.color = edge_color,
            #edge.width = 5*(E(g)$Synergy_score-ming)/(max(E(g)$Synergy_score)-ming)+2,
            edge.curved=gcurve)
  return(p)
}
