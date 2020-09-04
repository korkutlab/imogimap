library(data.table)

setwd("imogene/")
icp_gene_list <-  scan("data-raw/genes_Immune_checkpoints.txt",what = "charachter")



dfCT  <- read.table("data-raw/pancan_samples.txt", sep = "\t", header=TRUE)
colnames(dfCT)[2] <- "Tumor_Sample_ID"
colnames(dfCT)[3] <- "Cohort"
sel <- as.character(dfCT$Tumor_Sample_ID)
dfCT$SUBTYPE <- ifelse(dfCT$SUBTYPE=="Not_Applicable",as.character("all"),as.character(dfCT$SUBTYPE))

dfmeth <- read.table("data-raw/DNA-methylation-based-immune-infiltration-scores.tsv",sep = '\t',header=F)
dfmeth <- dfmeth[,-1]
colnames(dfmeth) <- c("Tumor_Sample_ID_full","Leukocyte_fraction")
dfmeth$Tumor_Sample_ID <- substr(dfmeth$Tumor_Sample_ID_full, 1, 15)
dfmeth <- setDT(dfmeth)[, lapply(.SD, mean), by=c(names(dfmeth)[3]), .SDcols=2]
dfmeth <- dfmeth[which(dfmeth$Tumor_Sample_ID %in% sel),]
rownames(dfmeth)<-dfmeth$Tumor_Sample_ID
TCGA_Leukocyte_fraction <- dfmeth


EMT_rnaseq <- read.csv("data-raw/EMTscore_pancan_RNAseq_from_TongPan_original.csv")
EMT_rnaseq$Tumor_Sample_ID <- substr(EMT_rnaseq$PatientID, 1, 15)
EMT_rnaseq <- EMT_rnaseq[EMT_rnaseq$Tumor_Sample_ID %in% sel,]
EMT_rnaseq <- setDT(EMT_rnaseq)[, lapply(.SD, mean), by=c(names(EMT_rnaseq)[3]), .SDcols=2]
row.names(EMT_rnaseq) <- EMT_rnaseq$Tumor_Sample_ID
TCGA_EMT <- EMT_rnaseq

TCGA_IMCell_fraction <- read.csv("data-raw/ICT_fractions.csv")

lgn_ang <- scan("data-raw/genes_angiogenesis.txt",what = "charachter")

EMT_gene_list <- read.csv("data-raw/Pan-Cancer-EMT-Signature-Genes.csv")

EMT_gene_list$sign <- ifelse(EMT_gene_list$Group=="M",1,-1)
colnames(EMT_gene_list)[1]<-"genes"

dft<- read.csv("data-raw/sample_immune_cell_fraction.csv")
rownames(dft) <- substr(dft$X,1,12)
dft$X <- NULL
sample_immune_cell_fraction_data <- dft

dft<- read.csv("data-raw/sample_Leukocyte fraction.csv")
dft$ID <- substr(dft$X,1,12)
rownames(dft) <- substr(dft$X,1,12)
dft$X <- NULL
sample_Leukocyte_fraction_data <-dft

dft<- read.csv("data-raw/sample_mRNA.csv",check.names = F)
colnames(dft)<- substr(colnames(dft),1,12)
i<-which(dft[,1] %like% "SLC35E2")
dft[i,1] <- c("SLC35E2B","SLC35E2A")
dft <-dft[complete.cases(dft[,1]),]
rownames(dft) <- dft[,1]
dft <-dft[,-1]
sample_mRNA_data <- dft
#sample_mRNA_data <- log2(dft+1)
usethis::use_data(TCGA_Leukocyte_fraction,TCGA_EMT,EMT_gene_list,icp_gene_list, TCGA_IMCell_fraction, overwrite = T)
usethis::use_data(sample_mRNA_data,sample_Leukocyte_fraction_data,sample_immune_cell_fraction_data, overwrite = T)

dft <- read.csv("data-raw/TCGA_disease.csv")
dft <- dft[-1,]
TCGA_disease_list  <- dft$x

usethis::use_data(TCGA_disease_list, overwrite = T)

document()
setwd("..")
