#' imogimap: Immuno-oncology gene expression assessment.
#'
#' The imogimap package provides three categories of functions:
#'
#' @section Functions to use public TCGA data:
#'
#' The following functions use \code{\link[curatedTCGAData]{curatedTCGAData}} to obtain publicly available RNASeq2GeneNorm data from The Cancer Genome Atlas (TCGA).
#'
#' \itemize{
##'  \item{\code{\link[imogimap]{im_syng_tcga}}:}{ Finds combinatorial association of immunotherapy co-targets with tumor intrinsic features as listed in TCGA_immune_features_list.}
##'  \item{\code{\link[imogimap]{im_boxplot_tcga}}:}{ Generates a stratified boxplot of immune feature values based on two genes}
#'   \item{\code{\link[imogimap]{im_cor_tcga}}:}{ Finds Spearman correlation of co-target genes with immuno-oncology feature}
#'}
#' @section  Functions to use with user-curated data:
#'
#'The following functions accept user-curated expression data as input.
#'
#' \itemize{
#'  \item{\code{\link[imogimap]{im_syng}}:}{ Finds combinatorial association of immunotherapy co-targets with tumor intrinsic features as listed in TCGA_immune_features_list.}
#'  \item{\code{\link[imogimap]{im_boxplot}}:}{ Generates a stratified boxplot of immune feature values based on two genes}
#'   \item{\code{\link[imogimap]{im_cor}}:}{ Finds Spearman correlation of co-target genes with immuno-oncology feature}
#'}
#'
#' @section Functions for output visualization and helper functions:
#' \itemize{
#' \item{\code{\link[imogimap]{im_netplot}}:}{ Creates and plot an igraph network from dataframe as outputed }
#' \item{\code{\link[imogimap]{find_a_synergy}}:}{ Evalues synergies for a data frame containing values of an immune feature and stratified expression levels of two genes.}
#' \item{\code{\link[imogimap]{calculate_syng_score}}:}{ Calculates a synergy score given median and standard error of means of four stratified groups}
#' \item{\code{\link[imogimap]{get_angio_score}}:}{ Calculates angiogenesis score}
#' \item{\code{\link[imogimap]{get_emt_score}}:}{ Calculates EMT score}
#' \item{\code{\link[imogimap]{get_ifng_score}}:}{ Calculates IFNG score}
#' \item{\code{\link[imogimap]{get_quantile_rank}}:}{ Stratify each column using quantile values.}
#'}
#'
#' @section Sample data sets and gene signature lists:
#'\itemize{
#'\item{\code{\link[imogimap]{AG_gene_list}}:}{ Angiogenesis signature genes}
#'\item{\code{\link[imogimap]{EMT_gene_list}}:}{ EMT signature genes}
#'\item{\code{\link[imogimap]{icp_gene_list}}:}{ Default list of immune checkpoints}
#'\item{\code{\link[imogimap]{sample_immune_cell_fraction_data}}:}{ Sample immune cell fraction data}
#'\item{\code{\link[imogimap]{sample_Leukocyte_fraction_data}}:}{ Sample Leukocyte fraction data}
#'\item{\code{\link[imogimap]{sample_mRNA_data}}:}{ Sample expression data}
#'\item{\code{\link[imogimap]{TCGA_Leukocyte_fraction}}:}{ Leukocyte fraction for TCGA data}
#'\item{\code{\link[imogimap]{TCGA_IMCell_fraction}}:}{ Immune cell type fractions for TCGA data}
#'\item{\code{\link[imogimap]{TCGA_TMB}}:}{ Tumor mutation burden for TCGA data}
#'\item{\code{\link[imogimap]{TCGA_immune_features_list}}:}{ Default immune features for TCGA data}
#'\item{\code{\link[imogimap]{TCGA_disease_list}}:}{ TCGA abbreviated disease names}
#
#'
#'}
#' @docType package
#' @name imogimap
NULL
#> NULL
