# Zoospore-RNA-Analysis

Outline:
Part -1: Clean up raw RNAseq data
Part 0: Perform differential gene expression (DGE) analysis of raw RNAseq data using DESeq2 in R
Part 1: Process DESeq2 data and define KOG, GO, IPR, KEGG gene annotations
Part 2: Format DESeq2 data, define additional gene annotations, and generate Fisher's Exact Statistics for DE of gene annotations
Part 3: Define additional CAZyme and cellulosome gene annotations




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
