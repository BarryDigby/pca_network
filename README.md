# pca_network
Analysis of public &amp; local  PCa (circRNA, miRNA, mRNA) datasets 

### circRNA analysis

Overlapping circRNAs detected by GSE113153, GSE118959 are detected using `circrna_intersection.R`:

- Clone 9 vs. Control: up-regulated 20, down-regulated 79.
- Clone 1 vs. Control: up-regulated 649, down-regulated 771
- Gleason >8 vs. Gleason <8: up-regulated 2761, down-regulated 2245

Using Clone 1 results, we extracted circRNAs that are dysregulated in high enzalutamide resistance that are also present in moderate enzalutamide resistance, and Gleason scores of above 8 using the following R logic: `intersection = master %>% group_by(circbaseID) %>% filter(length(unique(experiment))>=2 & "CLONE1" %in% experiment) %>% ungroup()`. Subsequent filtering to enforce fold change directionality produced 289 unique circRNAs. 


[circ_intersection_pdf](results/circrna_intersection.pdf)

### miRNA targets of DE circRNAs

Circbank provided 16844373 predicted circRNA-miRNA interactions, whilst CDSC (after mapping the genomic coordicnates to circBase IDs) provided 4773949 predicted circRNA-miRNA interactions. After intersection of differentially expressed circRNAs with predicted circRNA-miRNA interactions, a total of 36513 circRNA-miRNA interactions were returned (278 unique circRNAs, 2754 unique miRNAs).


### miRNA analysis

Significant reformatting was required for the miRNA microarray probes which utilised outdated miRNA IDs. The miRBase aliases file was used to perform this map ([https://www.mirbase.org/ftp.shtml](https://www.mirbase.org/ftp.shtml)) and correctly annotate the differentially expressed miRNAs. Pre-processing and alignment steps performed on raw small RNA-Seq data can be found in `scripts/alignment.steps.sh`.

- Clone 9 vs control: up-regulated 4, down-regulated 1
- Clone 1 vs control: up-regulated 146, down-regulated 110
- Tumor vs. normal (GSE21036): up-regulated 155, down-regulated 149
- Matched Tumor vs. normal (GSE23022): up-regulated 17, down-regulated 8
- Matched Tumor vs. normal (GSE36803): up-regulated 36, down-regulated 92
- Tumor vs. Normal (GSE45604): up-regulated 22, down-regulated 67
- Tumor vs. Normal (GSE46738): up-regulated 68, down-regulated 68
- Tumor vs. Normal (TCGA): up-regulated 130, down-regulated 128

Once again the logic `intersection = intersection %>% group_by(miRNA) %>% filter(length(unique(experiment))>=2 & "CLONE1" %in% experiment) %>% ungroup()` was used to find miRNAs in high enzalutamide resistance that are present in other datasets. After filtering for common fold-change direction, the intersection returned 42 unique miRNAs. 


### circRNA - miRNA network

A circRNA - miRNA network was constructed by filtering the predicted circRNA - miRNA interactions given by CircBank and CDSC using the differentially expressed miRNAs (42) resulting in a circRNA - miRNA network coupled with fold change values that has been derived from our differential expression analysis.

The resulting network hosts 195 circRNAs and 41 miRNAs which represents a loss of circRNA information from the initial DE results (278 circRNAs).

### mRNA analysis

- Clone1 vs control: up-regulated 4686, down-regulated 4379
- Clone 9 vs control: up-regulated 1420, down-regulated 1303
- Primary (castration) treatment vs control (GSE88752): up-regulated 2745, down-regulated 3013
- Secondary (enzalutamide) treatment vs control (GSE88752): up-regulated 3564, down-regulated 2942
- Day 7 enzalutamide vs. control (GSE143408): up-regulated 3829, down-regulated 3458
- Day 14 enzalutamide vs. control (GSE143408): up-regulated 5693, down-regulated 4792
- Day 21 enzalutamide vs. control (GSE143408): up-regulated 6128, down-regulated 4977
- 6 month enzalutamide vs. control (GSE78201): up-regulated 518, down-regulated 571

All of the above datasets involve enzalutamide treatment - really robust signature derived that is common to enzalutamide resistance in prostate cancer. Overlap these genes with TCGA PRAD data to identify if these signatures are present in canonical prostate cancer (like circrna, mirna overlaps).

- Tumor vs normal (TCGA): up-regulated 3860, down-regulated 4506 

We now required genes to be present in more than three experiments, however they must be present in the TCGA and clone1 results (i.e validated as enzalutamide resistance signature by one of other experiments). In this way we can identify genes involved in enzalutamide resistance that are present in canonical tumor samples (TCGA). This results in 359 unique genes. 

### mRNA targets of DE miRNAs

Mirbase (3340308), miRNet (43847 - user submitted list of DE miRNAs) and mirtarbase (937083) all provided predicted miRNA - mRNA interactions. mirBase genes were provided as refseq identifiers and were converted accordingly using biomart. The resulting intersection with differentially expressed miRNAs resulted in 57179 predicted miRNA - mRNA interactions (41 unique miRNAs, 15249 unique genes) which will be subset downstream to contain differentially expressed mRNAs returned by the mRNA analysis section above. 

### circRNA miRNA mRNA network

The previously generated circRNA - miRNA network was loaded (195 circRNAs and 41 miRNA) along with the predicted miRNA - mRNA interactions. Differentially expressed mRNAs were used to subset the predicted miRNA - mRNA interactions, and along with average mrna logFC, were added back to the circRNA - miRNA network. 

The final network has 195 circRNAs, 41 miRNAs and 307 genes. 

Further filtering is required before finalising the ceRNA network:

1. the higher the expression of circRNA, the lower the expression of targeted miRNAs and the higher the expression of downstream mRNAs (2481 results - 98 circRNAs, 24 miRNAs, 129 genes)

2. the lower the expression of circRNA, the higher the expression of targeted miRNAs and the lower the expression of downstream mRNAs. (861 results - 43 circRNAs, 16 miRNAs, 112 genes)

 

