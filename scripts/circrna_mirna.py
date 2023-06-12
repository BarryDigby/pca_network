#!/usr/bin/env python3

import pandas as pd

""" 
placeholder code, add to this once clone1_vs_all or clone1_vs_control is discussed.
"""

# Open GSE113153 limma results and read its contents
GSE113153 = pd.read_csv(
    "/data/github/pca_network/circrna/GSE113153/high_vs_low.txt", sep="\t", header=0
)

# Open GSE118959 limma results and read its contents
GSE118959_clone1 = pd.read_csv(
    "/data/github/pca_network/circrna/GSE118959/clone1_vs_control.txt",
    sep="\t",
    header=0,
)

# Open GSE118959 limma results (Clone 9) and read contents
GSE118959_clone9 = pd.read_csv(
    "/data/github/pca_network/circrna/GSE118959/clone9_vs_control.txt",
    sep="\t",
    header=0
)

# Open circbank file and read its contents
circbank = pd.read_csv(
    "/data/github/pca_network/data/circbank_miRNA_all_v1.txt", sep="\t", header=0
)

# Identify common circRNAs between studies. Use hsa_circ_0001 nomenclature.
keep = ["logFC", "adj.P.Val", "B", "Alias", "GeneSymbol"]

GSE113153 = GSE113153[keep].copy()
GSE118959_C1 = GSE118959_clone1[keep].copy()
GSE118959_C9 = GSE118959_clone9[keep].copy()

# Print % rows with missing Alias values (needed for circbank).

print(
    f"Dataset GSE113153 (high vs. low) contains {GSE113153['Alias'].isna().sum()} NaN values. The dataset contains a total of {GSE113153.shape[0]} values.\n\nMissing percentage:{(GSE113153['Alias'].isna().sum()/GSE113153.shape[0])*100}%"
)

print(
    f"Dataset GSE118959 (clone1 vs. control) contains {GSE118959_C1['Alias'].isna().sum()} NaN values. The dataset contains a total of {GSE118959_C1.shape[0]} values.\n\nMissing percentage:{(GSE118959_C1['Alias'].isna().sum()/GSE118959_C1.shape[0])*100}%"
)

print(
    f"Dataset GSE118959 (clone9 vs. control) contains {GSE118959_C9['Alias'].isna().sum()} NaN values. The dataset contains a total of {GSE118959_C9.shape[0]} values.\n\nMissing percentage:{(GSE118959_C9['Alias'].isna().sum()/GSE118959_C9.shape[0])*100}%"
)

# Drop these NaN columns
GSE113153.dropna(subset=['Alias'], inplace=True)
GSE118959_C1.dropna(subset=['Alias'], inplace=True)
GSE118959_C9.dropna(subset=['Alias'], inplace=True)

GSE113153.set_index('Alias', inplace=True)
GSE118959_C1.set_index('Alias', inplace=True)
GSE118959_C9.set_index('Alias', inplace=True)

GSE113153['experiment'] = "GSE113153"
GSE118959_C1['experiment'] = "GSE118959_Clone1"
GSE118959_C9['experiment'] = "GSE118959_Clone9"

dataframes = [GSE113153, GSE118959_C1, GSE118959_C9]

column_join = pd.concat(dataframes, join="outer").sort_index().reset_index()

column_join = column_join.groupby('Alias').filter(lambda x : x['Alias'].shape[0]>2)

column_join.to_csv("/data/github/pca_network/scripts/test.csv", index=False)

print(column_join.groupby('Alias').apply(lambda x: (x['logFC'] > 0).all() or (x['logFC'] < 0).all()))


# decide how to handle differing LFC vals? 66 circs overlapping all three exps... 

## 

""" merged_df = (
    pd.concat([GSE113153, GSE118959_C1, GSE118959_C9])
)

## prints duplicates rows
print(merged_df[merged_df.index.duplicated(keep=False)].sort_index())

# keep circs that have same fold change direction in both experiments... 
merged_sort = merged_df[merged_df.index.duplicated(keep=False)]
x = merged_sort.groupby(merged_sort.index)
same_sign = x.filter(lambda x: (x['logFC'] > 0).all() or (x['logFC'] < 0).all())
same_sign = same_sign.sort_index()

print(type(same_sign))

same_sign.to_csv("/data/github/pca_network/scripts/test.csv") """