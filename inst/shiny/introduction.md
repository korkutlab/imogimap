### How to use imogimap:

1. Use the "Input Gene" field to insert Hugo ID(s) of a single gene or a curated gene list representing a tumor-associated process (TAP). 
2. Use the default Immune checkpoint (ICP) genes that are already included or insert your own set of ICP genes.
3. Check "add receptor/ligand" box to automatically add genes corresponding to the known interacting molecules using cellphoneDB receptor-ligand database.  
4. Choose an immune associate phenotype (IAP) from the "immune phenotype"" drop-down menu, and specify the TCGA disease cohort that you are interested in. 
5. If you wish to include specificity and sensitivity statistical assessments, check the corresponding boxes and specify number of iterations. 
4. Hit the  "Submit" key and and go to the "Results" tab. 


Synergy score calculations and statistical analysis may be time-consuming for long gene lists. Imogimap has currently an hour time limit. For gene sets larger than 5 we highly recommend to  follow one of these workflows:

<ul><li> Conduct an initial analysis without specificity and sensitivity assessments, select TAP-ICP gene pairs of interest based on synergy scores and Wilcoxon p.values, and re-conduct separate analyses for each selected pair with specificity and sensitivity boxes checked. 

</li><li>Use the github link on top of this page to download our R software package and run the R software locally on your machine.

</li><li>Decrease number of iterations. The higher number of iterations we choose, the more accurate specificity pvalue and sensitivityR we will have. If your calculations are time consuming you can choose to decrease number of iterations for each assessment. We recommend number of iterations to be between 10-1000. 


See the "About" section for details on how to interpret the results 
