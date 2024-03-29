% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/im_syng.R
\name{im_syng}
\alias{im_syng}
\title{Find combinatorial association of immunotherapy co-targets with tumor intrinsic features.}
\usage{
im_syng(
  onco_gene,
  icp_gene,
  data_expression,
  data_feature,
  ndatamin = 8,
  add_features,
  method,
  specificity,
  N_iteration_specificity,
  sensitivity,
  N_iteration_sensitivity,
  add_receptor_ligand = FALSE
)
}
\arguments{
\item{onco_gene}{A character vector of gene IDs.}

\item{icp_gene}{An optional character vector of immune checkpoint gene/protein IDs.}

\item{data_expression}{A non-negative numeric matrix or data frame containing gene/protein expressions in linear scale.}

\item{data_feature}{An optional numeric matrix or data frame containing immune features.}

\item{ndatamin}{minimum number of samples. Synergy score calculation will be skipped for matrices with number of rows less than ndatamin}

\item{add_features}{An optional logical indicating if EMT score, angiogenesis score and IFNG expression should be added to immune features. Default is TRUE.}

\item{method}{A character string indicating which synergy score to be used. one of "max" or "independence". Default is "max".}

\item{specificity}{An optional logical indicating if specificity analysis should be done. Default is FALSE.}

\item{N_iteration_specificity}{Number of iterations for random sampling for specificity p.value calculation.
Default is 1000.}

\item{sensitivity}{An optional logical indicating if a sensitivity analysis should be done. Default is FALSE.}

\item{N_iteration_sensitivity}{Number of iterations for random sampling for sensitivity analysis.
Default is 1000.}

\item{add_receptor_ligand}{An optional logical indicating whether receptor_ligands pair should be added. Default is FALSE.}
}
\value{
a data.frames of synergy scores and bootstrapping p.values
}
\description{
Find combinatorial association of immunotherapy co-targets with tumor intrinsic features.
}
\details{
For details of synergy score calculations see get_syng_score function.

By default (if no icp_gene is specified), icp_gene_list will be used.

data_expression is formatted with genes/proteins as rows and samples/patients as columns.
For data_expression sample formats see sample_mRNA_data.

data_feature is formatted with samples/patients as rows and immune features as columns.
For data_feature sample format see sample_immune_cell_fraction_data.

For synergy score calculations all features are normalized to be on \code{[0,1]} range. For details of synergy score and significance pvalue calculations see \code{find_a_synergy} function.

A specificity p.value is computed using random sampling with replacement from two null models, generated from one of the two genes against a set of genes randomly selected from the genome. Two P-values are calculated for the synergistic interaction of the pair against the two null models. The highest of the two P-values is used to assess the specificity of the interaction against the whole genome. The number of randomly selected genes in each null model is determined by N_iteration_specificity.

Sensitivity (Robustness) score defined as normalized root mean square deviation of scores calculated over 70\% of samples, selected via random sampling. The number of sub-sample iterations is determined by N_iteration_sensitivity.
}
\examples{
im_syng(onco_gene =  c("BRAF"),
                  icp_gene= c("CD274","CTLA4"),
                  data_expression =  sample_mRNA_data,
                  data_feature = sample_Leukocyte_fraction_data)
}
\keyword{Synergy,}
\keyword{bootstrapping}
\keyword{checkpoint,}
\keyword{feature,}
\keyword{immune}
