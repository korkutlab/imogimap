% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{TCGA_TMB}
\alias{TCGA_TMB}
\title{Tumor mutation burden}
\format{
A data frame with 10123 rows and 3 columns:
\describe{
\item{\code{Tumor_Sample_ID}}{character TCGA 15-charachter Tumor sample ID}
\item{\code{TMB_Non.silent_per_Mb}}{numeric non silent tumor mutation burden}
\item{\code{TMB_Silent_per_Mb}}{numeric silent tumor mutation burden}
}
}
\source{
\url{https://doi.org/10.1016/j.cell.2018.03.022}
}
\usage{
TCGA_TMB
}
\description{
Silent and non-silent tumor mutation burden per MB for TCGA samples
}
\keyword{datasets}
