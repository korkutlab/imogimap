# imogene
# imogene is an R package that assesses user-provided genomic profiles against a curated list of immune checkpoints and identifies combinations that have noteworthy synergistic associations with tumor intrinsic immune processes. In the absence of user-provided data, Imogene uses number of API functions to assess TCGA mRNA data via publicly available CBioPortal repositories. 

# Some imogene functions import cBioPortalData, dplyr, tidyr libraries.

#To install from github
install_github('bozorgui/imogene')



#----------Use with sample data --------------

df <- im_cor(cotarget = "BRAF",
  checkpoint = "CD276",
  data_expression = sample_mRNA_data,
  data_feature = sample_immune_cell_fraction_data, 
  add_features = T)

df <- im_syng( cotarget  =c("BRAF","TGFB1"),
  checkpoint = c("CD274","CTLA4"),
  data_expression =sample_mRNA_data,
  data_feature = sample_Leukocyte_fraction_data,
  add_features = T,
  add_pvalue = T,
  N_iteration = 1000)


im_plot(cotarget = "BRAF",
  checkpoint = "CD276",
  data_expression = sample_mRNA_data,
  data_feature = sample_Leukocyte_fraction_data)

 
im_boot("ACVR2B","CD274",Immune_Feature = "Leukocyte_fraction",
  df_mrna = sample_mRNA_data,N_iteration=10000,
  df_lf = sample_Leukocyte_fraction_data,
  df_ict = sample_immune_cell_fraction_data)
  
  
  
#-------------Use with cbioportal data --------------

#Get synergy scores
df <- im_syng_tcga(cotarget = "BRAF",checkpoint = "CD276",cohort = "acc",add_pvalue = T,N_iteration = 100)

#Get correlations between gene sets
df <-im_cor_tcga(cotarget = "BRAF",cohort = "acc")

#Immune feature boxplots
im_plot_tcga(cotarget = "BRAF",
  checkpoint =  "CD276",
  cohort = "acc",
  Immune_Feature = "TMB_Silent_per_Mb")

#Get bootstrapping probability
im_boot_tcga(gene1 = "ACVRL1",gene2="CD274",cohort="lihc", 
  Immune_Feature="TMB_Silent_per_Mb", N_iteration=1000)

