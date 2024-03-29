% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_TCI_score.R
\name{get_TCI_score}
\alias{get_TCI_score}
\title{Calculates a T cell inflammation score}
\usage{
get_TCI_score(pdf)
}
\arguments{
\item{pdf}{A numeric data frame or matrix with samples as columns and gene/proteins as rows.}
}
\value{
a dataframe of T cell inflammation scores for each sample
}
\description{
Calculates a T cell inflammation score
}
\details{
The T cell inflammation score is defined as mean value of gene expressions in tcell_inflamed_gene_list. Scores are mean-centered and normalized by standard deviation based on distribution of scores in all samples in pdf.
}
\examples{
get_TCI_score(pdf=sample_mRNA_data)
}
\keyword{T}
\keyword{cell}
\keyword{inflammation,}
\keyword{pathway}
\keyword{score}
