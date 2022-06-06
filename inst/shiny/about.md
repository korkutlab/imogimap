# About

Immuno-oncology gene interaction Maps, ImogiMap, is a method and
informatics tool to automate the combinatorial searches for interactions between tumor-associated
and immune checkpoint processes. ImogiMap enables network analysis
of combinatorial interactions between Immune checkpoint (ICP) genes, genes that are related to
tumor-associated processes (TAP), and immune associated processes (IAPs). 

TAP genes may be any set of genes that are expressed in tumor ecosystems (likely but not necessarily in the tumor cells) and with functions involving hallmarks of cancer (e.g., proliferation, apoptosis, DNA repair, tumor metabolism, immune evasion). IAPs may involve processes that are related to immunotherapy responses and immune states in tumor ecosystems such as Immune cell infiltration (leukocyte fraction), Immune cell type fractions, tumor mutation burden, mRNA-based pathway scores for epithelial mesenchymal transition (EMT) status, vascularization, and IFNG expression.By mapping the ICP-TAP interactions, ImogiMap is able to inform on underlying mechanisms that jointly relate TAPs to therapeutically actionable ICPs and IAPs that are potential therapy response predictors.

## How to use imogimap
 
Use the "Input Gene" field to insert Hugo IDs of your curated list of TAP genes. Use the default ICP genes that are already included or insert your own set of ICP genes. Choose an IAP from the immune phenotype drop-down menu, and specify the TCGA disease cohort that you are interested in. Hit the  "Submit" key and and go to the "Results" tab. If you wish to include specificity and sensitivity statistical assessments, check the corresponding boxes and specify the number of iterations.


## How to interpret the results

Under the "Results" tab you will see a synergy score, and three statistical assessment measures for each gene pair. A positive synergy score with both genes marked as "Expressed" means synergistic increase of IAP level in response to expression of both genes. Similarly a negative synergy score with both genes marked as "Expressed" means synergistic decrease of IAP level in response to expression of both genes. 

A color-coded graphical network of combinatorial associations is shown in the "Network" section. Red (Blue) edges represent synergistic up-regulation (down-regulation) of IAP. Dark red (Dark blue) vertices identify over-expression (inhibition) of ICP genes, and orange (blue) vertices identify over-expression (low expression) of input genes. A cutoff parameter is included to subset the list to those pairs with absolute value of synergy scores higher than the cutoff. 

Choose a gene pair from the drop-down menu in the boxplot section to see boxplots of IAP values, stratified by gene expression levels. IAP and expression data are taken from the specified TCGA  disease cohort. 

For more information see  https://doi.org/10.1101/2021.10.06.462889

