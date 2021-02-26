
# imogene

<!-- badges: start -->
<!-- badges: end -->

The goal of imogene is to calculates statistical synergy scores based on gene expression profiles in multi-sample data to quantify combinatorial effects of oncogenes and immune checkpoints on immune related phenotype.


## Installation

You can install the released version of imogene from  GitHub with:

``` r
install_github('bozorgui/imogene')
```

## Usage

This is a basic example of how to use imogene with TCGA data

``` r
library(imogene)
# List Hugo ID's for onco-genes that you are interested in
my_onco <- c("TGFB1","TGFB2","TGFB3")

# List Hugo ID's for immune checkpoints that you are interested in. icp_gene_list can be used as default
my_icp <- icp_gene_list

#list TCGA abbreviated names and specify TCGA disease cohorts 
TCGA_disease_list
my_cohort <- c("luad","lusc")

#Calulate synergy scores, and optional pvalues, and variances
my_syng_df <-  im_syng_tcga(onco_gene  = my_onco,
         icp_gene = my_icp, 
         cohort = my_cohort, 
         method = "independence",
         add_pvalue = F, 
         N_iteration = 1000, 
         sensitivity = F)
```
