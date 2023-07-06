# pca_network
Analysis of public &amp; local  PCa (circRNA, miRNA, mRNA) datasets 

### circRNA analysis

Overlapping circRNAs detected by GSE113153, GSE118959 are detected using `circrna_intersection.R`. This script requires circrnas to be present in at least two experimental groups - clone1 vs control, clone9 vs control, gleason high vs gleason low. circrnas must also have the same fold change direction. 

circrnas detected via this method were overlapped with their mirna targets using `circrna_mirna_database.R`. 

### miRNA analysis

Overlapping mirnas detected by GSE21036, GSE23022 and GSE36803 were detected using `mirna_intersection.R`. The same filtering rules above applied here. 

### circrna-mirna merging

The output of `circrna_mirna_databse.R` containing the DE circrna seeds and their predicted mirna targets was read in R along with the results of `mirna_intersection.R`. In this way, we overlapped the predicted targets with the DE mirnas found in the experiments using `circrna_mirna_network.R`. This script outputs circrna - mirna pairs detected in our analysis, along with fold change values. 


