# Zoospore-RNA-Analysis

Outline:
Part -1: Clean up raw RNAseq data
Part 0: Perform differential gene expression (DGE) analysis of raw RNAseq data using DESeq2 in R
Part 1: Process DESeq2 data and define KOG, GO, IPR, KEGG gene annotations
Part 2: Format DESeq2 data, define additional gene annotations, and generate Fisher's Exact Statistics for DE of gene annotations
Part 3: Define additional CAZyme and cellulosome gene annotations



Part -1: Clean up raw RNAseq data

Required Inputs:
1) FASTA file for genes in Neocallimastix californiae G1 genome,'Neosp1_FilteredModels5_deflines_pre.fasta'
Downloaded from JGI Mycocosm, under Annotation --> Filtered Models ("best") --> Genes --> "Neosp1_FilteredModels5_deflines.gff.gz" for Project 1029446, 20171122.
The FASTA lists gene keys, identified by Mycocosm proteinID, followed by amino acid sequences

2) counts_in_filename = 'counts_RNAseq_original.txt'
Note, raw RNA-Seq reads are available at NCBI BioProject: PRJNA982907 to PRJNA982924 

3) tpm_counts_in_filename = 'tpm_counts_RNAseq_original.txt'
Transcripts-per-million normalization is considered in this study.

Outputs to temp_outputs:
counts_RNAseq_updated.xlsx
Neosp1_FilteredModels5_deflines_duped_proteinIDs_sorted.fasta
Neosp1_FilteredModels5_deflines_post.fasta
tpm_counts_RNAseq_updated


Part 0: DESeq2 in R




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
