# Zoospore-RNA-Analysis

Outline:
Part 1: "Zoospore_RNA_data_analysis_pt1_DESeq2.py"
Part 2: "Zoospore_RNA_data_analysis_pt2_Format.py"
Part 3: "Zoospore_RNA_data_analysis_pt3_CAZymes.py"

Part 1
Script: "Zoospore_RNA_data_analysis_pt1_DESeq2.py"

Required Inputs:
1) DESeq2 Output: "deseq2_output_20230315.csv" from DESeq2 analysis in R script.
Columns: proteinID, baseMean, mat_vs_zoosp_log2FC, lfcSE, stat, pvalue, mat_vs_zoosp_Padj
R v4.2.2
Bioconductor v3.16
DESeq2 v1.36.0

2) DESeq2 Normalized counts: "deseq2_normalized_counts_labeled_20230315.csv"

3) RNAseq counts: "updated_counts.csv"

4) 
Functions:

Outputs:
