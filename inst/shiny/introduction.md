## Getting started with Imogimap:

1.  Use the "Input Gene" field to insert case sensitive Hugo ID(s) of a single gene or a curated gene list representing a tumor-associated process (TAP).

2.  Use the default Immune checkpoint genes that are already included or insert your own set of ICP Hugo IDs.

3.  Check "add receptor/ligand" box to automatically add genes corresponding to the known interacting molecules using cellphoneDB receptor-ligand database.\

4.  Choose an immune associate phenotype (IAP) from the "immune phenotype"" drop-down menu, and specify the TCGA disease cohort that you are interested in.

5.  Choose a combined action score calculation method.

6.  Check sensitivity box and specify number of iterations, to statistically measure the sensitivity of scores to the exact data configuration. Note that sensitivity calculations are time-consuming and may exceed Imogimap time limit. See tips below for recommendations.

7.  Check specificity box and specify number of iterations, to measure the probability of finding a score higher than the observed combined action score. Note that specificity calculations are time-consuming and may exceed Imogimap time limit. See tips below for recommendations.

8.  Hit the "Submit" key and and go to the "Results" tab.

#### Tips:

The standard workflow for sensitivity and specificity assessments can sometimes be prohibitively time-consuming for gene lists longer than 10. To improve efficiency and run times we recommend to follow one of these alternative workflows:

    -   Evaluate in steps! Conduct an initial analysis with specificity and sensitivity boxes un-checked; Use the results to select gene pairs of interest based on the reported combined action scores and Wilcoxon p.values; Finally, for each selected gene pair, re-conduct separate analyses with specificity and sensitivity boxes checked.

    -   Decrease number of iterations! The higher number of iterations we choose, the more accurate specificity pvalue and sensitivityR we will have. If your calculations are time consuming you can choose to decrease number of iterations for each assessment. We recommend number of iterations to be between 10-1000.

    -   Run locally! Use the Imogimap github link to download the R software package and run it locally on your machine.
