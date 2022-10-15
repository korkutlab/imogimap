# Creates essential data sets from raw-data
library(data.table)
library(imogimap)
current_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_dir)

#List of TCGA diseases
dft <- read.csv("TCGA_disease.csv")
TCGA_disease_list  <- dft$x
usethis::use_data( TCGA_disease_list,overwrite = T)

#List of Immune Checkpoints
icp_gene_list <-  scan("genes_Immune_checkpoints.txt",what = "charachter")
#icp_gene_list[1] <- "C10orf54"
#icp_gene_list[12] <- "DKFZp686O24166"
icp_gene_list <- c(icp_gene_list,"CSF1","IL34","NCR3",
              "FLT3LG","CSF2","CSF3","SLC7A1",
              "SIRPA","CELSR3","BMP10","CCL4L2",
              "TNF","TNFRSF13B","TNFRSF17","GPRC5B",
              "FAM3C","KLRG2","HLA-DPA1","PVR",
              "IL13RA2","ADGRG5")
usethis::use_data(icp_gene_list,overwrite = T)


#Ligand-receptor pairs
tmp1 <- readxl::read_xlsx("cellphoneDB_interaction_input.xlsx")
tmp1<-tmp1[complete.cases(tmp1$protein_name_a),]
tmp1<-tmp1[complete.cases(tmp1$protein_name_b),]
tmp2 <- readxl::read_xlsx("cellphoneDB_gene_name.xlsx")
tmp1$Gene_a <- plyr::mapvalues(tmp1$partner_a,tmp2$uniprot,tmp2$hgnc_symbol)
tmp1$Gene_b <- plyr::mapvalues(tmp1$partner_b,tmp2$uniprot,tmp2$hgnc_symbol)
lgn_receptor_ligand <- tmp1[,c("Gene_a","Gene_b")]
colnames(lgn_receptor_ligand)<-c("Gene1","Gene2")
lgn_receptor_ligand$ligand_receptor_interaction <- TRUE

# lgn_receptor_ligand <- read.csv("ligand_receptor_curated.csv")
# lgn_receptor_ligand$X<-NULL
# lgn_receptor_ligand <- apply(lgn_receptor_ligand, 1, function(x) data.frame(t(x[order(x)])))
#
# lgn_receptor_ligand <- data.table::rbindlist(lgn_receptor_ligand,use.names = FALSE)
# colnames(lgn_receptor_ligand)<-c("Gene1","Gene2")
# lgn_receptor_ligand$ligand_receptor_interaction <- TRUE

usethis::use_data(lgn_receptor_ligand,overwrite = T)


#TCGA Immune Infiltration (Leukocyte fraction) scores
df_inf <- fread("TCGA_all_leuk_estimate.masked.20170107.tsv")
colnames(df_inf) <- c("disease","Tumor_Sample_ID_full","Leukocyte_fraction")
df_inf$Tumor_Sample_ID <- substr(df_inf$Tumor_Sample_ID_full, 1, 15)
#df_inf <- df_inf[which(df_inf$Tumor_Sample_ID %in% sel),]
df_inf[df_inf$Leukocyte_fraction < 0,]$Leukocyte_fraction <- 0
df_inf$Tumor_Sample_ID_full <- NULL
df_inf$disease <- NULL
df_inf<- df_inf[,c("Tumor_Sample_ID","Leukocyte_fraction")]
df_inf <- setDT(df_inf)[, lapply(.SD, median), by=c(names(df_inf)[1]), .SDcols=2]
TCGA_Leukocyte_fraction <- df_inf
usethis::use_data( TCGA_Leukocyte_fraction,overwrite = T)

#EMT gene list
EMT_gene_list <- read.csv("Pan-Cancer-EMT-Signature-Genes.csv")
EMT_gene_list$sign <- ifelse(EMT_gene_list$Group=="M",1,-1)
colnames(EMT_gene_list)[1]<-"genes"
colnames(EMT_gene_list)[2]<-"group"
usethis::use_data( EMT_gene_list,overwrite = T)

#Angiogenesis gene list
AG_gene_list <- scan("genes_angiogenesis.txt",what = "charachter")
usethis::use_data( AG_gene_list,overwrite = T)

#Tumor mutation burden
TCGA_TMB <- read.table("mutation-load_updated.txt",header=T)
TCGA_TMB <- TCGA_TMB[,c("Tumor_Sample_ID","Non.silent_per_Mb","Silent_per_Mb")]
colnames(TCGA_TMB) <- c("Tumor_Sample_ID","TMB_Non.silent_per_Mb","TMB_Silent_per_Mb")
usethis::use_data( TCGA_TMB,overwrite = T)


#TCGA Immune cell fractions from CIBERSORT
TCGA_IMCell_fraction <- read.csv("ICT_fractions.csv")
usethis::use_data( TCGA_IMCell_fraction,overwrite = T)


#Sample data sets
dft<- read.csv("sample_immune_cell_fraction.csv")
rownames(dft) <- substr(dft$X,1,12)
dft$X <- NULL
sample_immune_cell_fraction_data <- dft
usethis::use_data( sample_immune_cell_fraction_data,overwrite = T)

dft<- read.csv("sample_Leukocyte_fraction.csv")
rownames(dft) <- substr(dft$X,1,12)
dft$X <- NULL
colnames(dft)<- "Leukocyte_fraction"
sample_Leukocyte_fraction_data <-dft
usethis::use_data( sample_Leukocyte_fraction_data,overwrite = T)


#sample mRNA data
sample_mRNA_data <- read.csv("sample_mRNA_data.csv",header = T,row.names = 1)
colnames(sample_mRNA_data) <- gsub(pattern = "\\.",replacement = "-",x=colnames(sample_mRNA_data))
usethis::use_data( sample_mRNA_data,overwrite = T)


dft <- read.csv("TCGA_immune_features.csv")
dft <- dft$x
TCGA_immune_features_list <- dft
usethis::use_data( TCGA_immune_features_list,overwrite = T)

# T cell inflamed
dft <- fread("t_inflamed_gene_expr.tsv")
tcell_inflamed_gene_list <- colnames(dft)[-c(1,20)]
usethis::use_data( tcell_inflamed_gene_list,overwrite = T)

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
usethis::use_data( lgn_tcell_tnbc,overwrite = T)

#Creates data

usethis::use_data( tcell_inflamed_gene_list,overwrite = T)
usethis::use_data( TCGA_immune_features_list,overwrite = T)

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


