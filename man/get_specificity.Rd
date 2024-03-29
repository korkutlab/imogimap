% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_specificity.R
\name{get_specificity}
\alias{get_specificity}
\title{Calculation and statistical assessment of synergistic associations using TCGA data.}
\usage{
get_specificity(
  df_syng,
  method = "max",
  ndatamin = 8,
  N_iteration_specificity = 10000
)
}
\arguments{
\item{df_syng}{a dataframe with 10 columns as outputed by im_syng_tcga.}

\item{method}{A character string indicating which synergy score to be used. one of "max" or "independence".
Default is "max".}

\item{ndatamin}{minimum number of samples. Synergy score calculation will be skipped for matrices with number of rows less than ndatamin}

\item{N_iteration_specificity}{Number of iterations for random sampling for specificity p.value calculation.
Default is 1000.}
}
\value{
Specificty pvalues for each row of dataframe
}
\description{
A helper function to allow user calculate specifity pvalue for a subset of its dataset. To calculate specificity pvalue for all data points set specificity = TRUE in im_syng_TCGA.
}
\details{
A specificity p.value is computed using random sampling with replacement from two null models, generated from one of the two genes against a set of genes randomly selected from the genome. Two P-values are calculated for the synergistic interaction of the pair against the two null models. The highest of the two P-values is used to assess the specificity of the interaction against the whole genome. The number of randomly selected genes in each null model is determined by N_iteration_specificity.
}
\examples{
df <- im_syng_tcga(onco_gene=c("TGFB1","SERPINB9"), cohort=c("ucec"),specificity = F)
df <- get_specificity(df)
}
\keyword{,bootstrapping}
\keyword{pvalue}
\keyword{specificity,}
