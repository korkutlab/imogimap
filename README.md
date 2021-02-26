
# imogene

<!-- badges: start -->
<!-- badges: end -->


The goal of imogene is to calculate statistical synergy scores based on mRNA expression profiles of multi-sample data to quantify combinatorial effects of single onco-genes or tumor intrinsic onco-genic pathways and immune checkpoints on immune related phenotype in tumor microenvironment.


## Installation

You can install the released version of imogene from  GitHub with:

``` r
install_github('bozorgui/imogene')
```

## Usage

This is a basic example of how to use imogene with TCGA data

``` r
library(imogene)
# List Hugo ID's for a single onco-gene or list or genes defining an onco-genic pathway  signature
my_onco <- c("TGFB1","TGFB2","TGFB3")

# List Hugo ID's for immune checkpoints that you are interested in. icp_gene_list can be used as default
my_icp <- icp_gene_list

#list TCGA abbreviated names and specify TCGA disease cohorts 
TCGA_disease_list
my_cohort <- c("luad","lusc")

#Calulate synergy scores, and optional pvalues, and variances, for combinatorial effects of all gene pairs on all immune phenotypes as listed in TCGA_immune_features_list
my_syng_df <-  im_syng_tcga(onco_gene  = my_onco,
         icp_gene = my_icp, 
         cohort = my_cohort, 
         method = "independence",
         add_pvalue = F, 
         N_iteration = 1000, 
         sensitivity = F)

#Generate stratified boxplot that represents data used to get a single synergy score.
im_boxplot_tcga(onco_gene = "TGFB1", icp_gene = "CD270",cohort = "luad",
Immune_Feature = "IFNGscore")
```
This is a basic example of how to use imogene with user's data
``` r
library(imogene)
# Use sample_mRNA_data directly or as formatting guide for expression data
my_expressions <- sample_mRNA_data

# Use sample_Leukocyte_fraction_data directly or as formatting guide for immune feature/phenotype data
my_features <- sample_Leukocyte_fraction_data

# List gene ID's for onco-genes that you are interested in, as listed in your data
my_onco <- c("TGFB1","TGFB2","TGFB3")

# List gene ID's for immune checkpoints that you are interested in. icp_gene_list can be used as default.
my_icp <- icp_gene_list

# Calulate synergy scores, and optional pvalues, and variances, for combinatorial effects of all gene pairs on immune features. 
df <- im_syng( onco_gene = my_onco,
             icp_gene = my_icp,
             data_expression = my_expressions ,
             data_feature = my_features,
             add_features = F,
             add_pvalue = T,
             N_iteration = 1000,
            sensitivity = T)

 
#Generate stratified boxplot that represents data used to get a single synergy score.
im_boxplot(cotarget = "TGFB1", 
  icp_gene = "CD276",
  data_expression = sample_mRNA_data,
  data_feature = sample_Leukocyte_fraction_data)
```

This is a basic example of generating graphical network for synergy scores.

```r       
#Generate graphical network based on synergy scores for a single immune feature
jpeg("syng_network.jpeg",width=1000,height=1000)
im_netplot(df = my_syng_df ,
           Immune_Feature = "IFNGscore",
           cutoff = 0.35,
           cohort = "luad" ,
           icp_gene = my_icp,
           seed = 123)
dev.off()
```
