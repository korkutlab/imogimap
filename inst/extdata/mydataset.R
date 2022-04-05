# Creates essential data sets from raw-data
library(data.table)

current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)

#List of TCGA diseases
dft <- read.csv("TCGA_disease.csv")
TCGA_disease_list  <- dft$x

#List of Immune Checkpoints
icp_gene_list <-  scan("genes_Immune_checkpoints.txt",what = "charachter")
icp_gene_list[1] <- "C10orf54"
icp_gene_list[12] <- "DKFZp686O24166"

#TCGA Immune Infiltration (Leukocyte fraction) scores
df_inf <- fread("TCGA_all_leuk_estimate.masked.20170107.tsv")
colnames(df_inf) <- c("disease","Tumor_Sample_ID_full","Leukocyte_fraction")
df_inf$Tumor_Sample_ID <- substr(df_inf$Tumor_Sample_ID_full, 1, 15)
df_inf <- df_inf[which(df_inf$Tumor_Sample_ID %in% sel),]
df_inf[df_inf$Leukocyte_fraction < 0,]$Leukocyte_fraction <- 0
df_inf$Tumor_Sample_ID_full <- NULL
df_inf$disease <- NULL
df_inf<- df_inf[,c("Tumor_Sample_ID","Leukocyte_fraction")]
df_inf <- setDT(df_inf)[, lapply(.SD, median), by=c(names(df_inf)[1]), .SDcols=2]
TCGA_Leukocyte_fraction <- df_inf

#EMT gene list
EMT_gene_list <- read.csv("Pan-Cancer-EMT-Signature-Genes.csv")
EMT_gene_list$sign <- ifelse(EMT_gene_list$Group=="M",1,-1)
colnames(EMT_gene_list)[1]<-"genes"
colnames(EMT_gene_list)[2]<-"group"
usethis::use_data( EMT_gene_list,overwrite = T)

#Angiogenesis gene list
AG_gene_list <- scan("genes_angiogenesis.txt",what = "charachter")

#Tumor mutation burden
TCGA_TMB <- read.table("mutation-load_updated.txt",header=T)
TCGA_TMB <- TCGA_TMB[,c("Tumor_Sample_ID","Non.silent_per_Mb","Silent_per_Mb")]
colnames(TCGA_TMB) <- c("Tumor_Sample_ID","TMB_Non.silent_per_Mb","TMB_Silent_per_Mb")


#TCGA Immune cell fractions from CIBERSORT
TCGA_IMCell_fraction <- read.csv("ICT_fractions.csv")


#Sample data sets
dft<- read.csv("sample_immune_cell_fraction.csv")
rownames(dft) <- substr(dft$X,1,12)
dft$X <- NULL
sample_immune_cell_fraction_data <- dft

dft<- read.csv("sample_Leukocyte_fraction.csv")
rownames(dft) <- substr(dft$X,1,12)
dft$X <- NULL
colnames(dft)<- "Leukocyte_fraction"
sample_Leukocyte_fraction_data <-dft

#sample mRNA data
sample_mRNA_data <- read.csv("sample_mRNA_data.csv",header = T,row.names = 1)
colnames(sample_mRNA_data) <- gsub(pattern = "\\.",replacement = "-",x=colnames(sample_mRNA_data))


dft <- read.csv("TCGA_immune_features.csv")
dft <- dft$x
TCGA_immune_features_list <- dft


# T cell dysfunction signature genes

lgn1 <-
readxl::read_xlsx("significant_dysfunction_scores.xlsx",sheet = 1)
lgn2 <-
readxl::read_xlsx("significant_dysfunction_scores.xlsx",sheet = 2)
lgn3 <-
readxl::read_xlsx("significant_dysfunction_scores.xlsx",sheet = 3)
lgn4 <-
readxl::read_xlsx("significant_dysfunction_scores.xlsx",sheet = 4)
lgn5 <-
readxl::read_xlsx("significant_dysfunction_scores.xlsx",sheet = 5)

lgn <-  rbind(lgn1,lgn2,lgn3,lgn4,lgn5) 
lgn <- lgn[lgn$FDR<0.1,]
library(plyr) 
lgn_count<- ddply(lgn,.(Symbol),nrow) 
lgn_count <-lgn_count[lgn_count$V1>1,] 
lgn <- lgn[lgn$Symbol %in% lgn_count$Symbol,] 
lgn <- unique(lgn$Symbol)

lgn5_sub <- lgn5[lgn5$Symbol %in% lgn,]
lgn5_sub <-lgn5_sub[lgn5_sub$FDR<0.05,]

lgn2_sub <- lgn2[lgn2$Symbol %in% lgn,]
lgn2_sub <-lgn2_sub[lgn2_sub$FDR<0.05,]
lgn2_sub <-lgn2[lgn2$FDR<0.05,]
lgn_tcell_tnbc <- lgn5_sub$Symbol

#Creates data

setwd("../")
usethis::use_data( sample_mRNA_data,overwrite = T)

usethis::use_data(
  TCGA_Leukocyte_fraction,
  EMT_gene_list,
  icp_gene_list,
  TCGA_IMCell_fraction,
  sample_mRNA_data,
  sample_Leukocyte_fraction_data,
  sample_immune_cell_fraction_data,
  TCGA_disease_list,
  TCGA_TMB,
  AG_gene_list,
  overwrite = T)

document()


