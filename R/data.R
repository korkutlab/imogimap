#' @title TCGA Leukocyte fractions
#' @description Leukocyte fraction of samples in TCGA data
#' @format A data frame with 8775 rows and 3 variables:
#' \describe{
#'   \item{\code{Tumor_Sample_ID}}{character TCGA 15-charachter Tumor sample ID}
#'   \item{\code{Leukocyte_fraction}}{double Leukocyte fraction of cells}
#'}
#' @source \url{https://doi.org/10.1016/j.cell.2018.03.022}
"TCGA_Leukocyte_fraction"

#' @title EMT signature genes
#' @description  A set of 77 unique genes as the pan-cancer EMT signature
#' @format A data frame with 77 rows and 3 columns:
#' \describe{
#'   \item{\code{genes}}{character Gene hugo symbol}
#'   \item{\code{group}}{character Epithelial vs Mesenchymal}
#'   \item{\code{sign}}{numeric Negative for Epithelial and Positive for Mesenchymal}
#'}
#' @source \url{https://doi.org/10.1158/1078-0432.CCR-15-0876}
"EMT_gene_list"

#' @title Angiogenesis signature genes
#' @description  A set of 41 unique gene Hugo IDs as the core human Angiogenesis signature
#' @format A vector of characters with length 41:
#' \describe{
#'   Each character is a gene hugo symbol
#'}
#' @source \url{https://doi.org/10.1016/j.ccr.2013.06.004}
"AG_gene_list"

#' @title Tumor mutation burden
#' @description  Silent and non-silent tumor mutation burden per MB for TCGA samples
#' @format A data frame with 10123 rows and 3 columns:
#' \describe{
#'   \item{\code{Tumor_Sample_ID}}{character TCGA 15-charachter Tumor sample ID}
#'   \item{\code{TMB_Non.silent_per_Mb}}{numeric non silent tumor mutation burden}
#'   \item{\code{TMB_Silent_per_Mb}}{numeric silent tumor mutation burden}
#'}
#' @source \url{https://doi.org/10.1016/j.cell.2018.03.022}
"TCGA_TMB"

#' @title Immune cell type fractions for TCGA data
#' @description CIBERSORT fraction of each 21 immune cell type within immune population. Note that The values are not absolute values and are normalized by the sum of all 21 types as provided by CIBERSORT.
#' @format A data frame with 11080 rows and 23 variables:
#' \describe{
#'   \item{\code{PATIENT_BARCODE}}{character TCGA patient barcode}
#'   \item{\code{B.Cells.Memory}}{double fraction of Memory B Cells}
#'   \item{\code{B.Cells.Naive}}{double Naive B cells}
#'   \item{\code{Plasma.Cells}}{double Plasma cells}
#'   \item{\code{T.Cells.CD8}}{double CD8 T Cells}
#'   \item{\code{T.Cells.CD4.Naive}}{double Naive CD4 T Cells}
#'   \item{\code{T.Cells.CD4.Memory.Activated}}{double Activated Memory CD4 T cells}
#'   \item{\code{T.Cells.CD4.Memory.Resting}}{double Resting Memory CD4 T cells}
#'   \item{\code{T.Cells.Follicular.Helper}}{double Helper Follicular T Cells}
#'   \item{\code{T.Cells.Regulatory.Tregs}}{double Regulatory T cells}
#'   \item{\code{T.Cells.gamma.delta}}{double Gamma Delta T Cells}
#'   \item{\code{NK.Cells.Activated}}{double Activated NK Cells}
#'   \item{\code{NK.Cells.Resting}}{double Resting NK Cells}
#'   \item{\code{Monocytes}}{double Monocytes}
#'   \item{\code{Macrophages.M0}}{double M0 Macrophages}
#'   \item{\code{Macrophages.M1}}{double M1 Macrophages}
#'   \item{\code{Macrophages.M2}}{double M2 Macrophages}
#'   \item{\code{Dendritic.Cells.Resting}}{double Resting Dendritic Cells}
#'   \item{\code{Dendritic.Cells.Activated}}{double Activated Dendritic Cells}
#'   \item{\code{Mast.Cells.Activated}}{double Activated Mast Cells}
#'   \item{\code{Mast.Cells.Resting}}{double Resting Mast Cells}
#'   \item{\code{Eosinophils}}{double Eosinophils}
#'   \item{\code{Neutrophils}}{double Neutrophils}
#'}
#' @source \url{https://cibersort.stanford.edu}
"TCGA_IMCell_fraction"

#' @title Default list of immune checkpoints
#' @description List of 50 genes including 29 immune checkpoints and their 21 interacting ligands/receptors
#' @format A vector of characters
#' \describe{
#'   Each character is a gene hugo symbol
#'}
#' @source \url{https://doi.org/10.1101/2021.10.06.462889}
"icp_gene_list"

#' @title Sample data set for immune cell fraction data
#' @description Data frame containing immune cell fractions for a subset of TCGA cohorts to serve as formatting guide for user-curated data or in test runs.
#' @format A data frame with of with 78 rows and 22 columns. Row names are a patient IDs and column names are cell types.
#' @source \url{https://cibersort.stanford.edu}
"sample_immune_cell_fraction_data"


#' @title Sample Leukocyte fraction data
#' @description Data frame containing Leukocyte fractions for a subset of TCGA cohorts to serve as formatting guide for user-curated data or in test runs.
#' @format A data frame with of with 76 rows and 1 column. Row names are a patient IDs.
#' @source \url{https://doi.org/10.1016/j.cell.2018.03.022}
"sample_Leukocyte_fraction_data"


#' @title Sample mRNA expression data
#' @description Data frame containing RNASeq2GeneNorm expressions for a subset of TCGA cohorts to serve as formatting guide for user-curated data or in test runs.
#' @format A data frame with of with 19958 rows and 78 columns. Row names are gene Hugo IDs and column names are patient IDs.
#' @source \code{\link[curatedTCGAData]{curatedTCGAData}}
"sample_mRNA_data"



#' @title TCGA disease names
#' @description TCGA abbreviated disease names to be used as input in TCGA functions.
#' @format A character vector of length 32
#' @source \code{\link[curatedTCGAData]{curatedTCGAData}}
"TCGA_disease_list"


#' @title TCGA immune features list
#' @description Lists immune feature names that are default to TCGA functions.
#' @format A character vector of length 28
"TCGA_immune_features_list"

globalVariables(c("AG_gene_list", "EMT_gene_list", "TCGA_TMB","TCGA_immune_features_list",
                  "TCGA_disease_list","sample_mRNA_data","sample_Leukocyte_fraction_data",
                  "sample_immune_cell_fraction_data","sample_immune_cell_fraction_data",
                  "icp_gene_list","TCGA_IMCell_fraction","TCGA_Leukocyte_fraction"))
globalVariables(c("."))
