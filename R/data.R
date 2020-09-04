#' @title TCGA EMT scores
#' @description Epithelial–mesenchymal transition score for TCGA data. Higher EMT scores correspond to higher expression of mesenchymal genes.
#' @format A data frame with 9120 rows and 3 variables:
#' \describe{
#'   \item{\code{Tumor_Sample_ID}}{character TCGA 15-charachter Tumor sample ID}
#'   \item{\code{EMTscore}}{double  EMT score}
#'}
#' @source from Tong Pan
"TCGA_EMT"

#' @title TCGA Leukocyte fractions
#' @description Leukocyte fraction of samples in TCGA data
#' @format A data frame with 8775 rows and 3 variables:
#' \describe{
#'   \item{\code{Tumor_Sample_ID}}{character TCGA 15-charachter Tumor sample ID}
#'   \item{\code{Leukocyte_fraction}}{double Leukocyte fraction of cells}
#'}
#' @source \url{http://somewhere.important.com/}
"TCGA_Leukocyte_fraction"

#' @title EMT signature genes
#' @description  A set of 77 unique genes as the pan-cancer EMT signature
#' @format A data frame with 77 rows and 2 variables:
#' \describe{
#'   \item{\code{Gene}}{character Gene hugo symbol}
#'   \item{\code{Group}}{character Epithelial vs Mesenchymal}
#'   \item{\code{signe}}{numeric Negative for Epithelial and Positive for Mesenchymal}
#'}
#' @source \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4737991/#SD1}
"EMT_gene_list"

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

