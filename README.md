#  Transcriptomic Analysis of Anaerobic Gut Fungal Life Stages: Zoospores vs Mats

## Background
These scripts generate **global gene expression profiles** and perform **differential gene expression analysis** for zoospore vs mat sample groups.
In this project, **RNA-Seq data** was acquired for two sample groups from the **anaerobic gut fungi** (phylum Neocallimastigomycota) *Neocalliamstix californiae* G1 [[link]](https://mycocosm.jgi.doe.gov/Neosp1/Neosp1.home.html): 
1. Culture cell pellets enriched in **zoospores,** the young life stage of anaerobic gut fungi
2. **Fungal mats** with mixed life stages, including sporangia, the mature life stage of anaerobic gut fungi

The main steps are as follows:
1. Prepare RNA-Seq data for DESeq2 differential gene expression analysis.
2. Perform DESeq2 for log2 fold-change and p-values for differential expression of each gene between the two sample groups.
3. Align functional gene annotations, downloaded from the Joint Genome Institute (JGI) Mycocosm portal [[link]](https://mycocosm.jgi.doe.gov/Neosp1/Neosp1.home.html).
4. Define cutoffs for significant differential regulation of individual genes.
5. For each annotation, generate gene counts for upregulation/downregulation of the annotated genes in zoospores vs. mats. Additionally, generate Fisher's Exact Statistics to describe differential regulation of annotations (distinct from differential regulation of individual genes).
6. Align additional gene annotations, specific to system of study (anaerobic gut fungi).
7. Generate figures to visualize statistics, ie: volcano plots.

## Associated Publication
Please see the associated publication for more project-specific details:

Lazarina V. Butkovich, Patrick A. Leggieri, Stephen P. Lillington, Tejas A. Navaratna, Candice L. Swift, Nikola G. Malinov, Thea R. Zalunardo, Oliver B. Vining, Anna Lipzen, Mei Wang, Juying Yan, Vivian Ng, Igor Grigoriev, Michelle A. O’Malley, **"Separation of life stages within anaerobic fungi highlights differences in global transcription and metabolism."** (in preparation).

## Installation
### Dependencies:
R: rtools, jsonlite, rlang, BiocManager (DESeq2)
Python: pandas, scipy, bioinfokit, openpyxl, xlsxwriter
Download Roboto font from Google fonts, OR change fonts in Python plots

## Description of Scripts
This project includes multiple scripts:
| Script Name                                       | Description                               |
| ------------------------------------------------- | ----------------------------------------- |
| Zoospore_RNA_data_analysis_pt1_RNAseq_Cleanup.py | RNA-Seq raw data contains some transcripts with identical amino acid sequences. For the purposes of describing putative functions of the corresponding genes (ID'd by proteinID) to these transcripts, the RNA-Seq raw counts and TPM counts data are consolidated for groups of transcripts with identical amino acid sequence in this script. |
| Zoospore_RNA_data_analysis_pt2_DESeq2_in_R.R | In order to perform differential gene expression (DGE) analysis, this script runs DESeq2 in R to generate log2-fold change data and p-values for each proteinID.  |
| Zoospore_RNA_data_analysis_pt3_DGE_Main_Annotations.py | This script aligns putative functional annotations from JGI Mycocosm to the proteinIDs considered in this dataset. The main annotations considered are KOG, GO, IPR, and KEGG. Significant upregulation of a gene in the zoospores is defined as follows: DESeq2 p-value < 0.05, log2 fold-change < 1, and average TPM in zoospore samples > 1. Significant upregulation of a gene in the mats is defined as follows: DESeq2 p-value < 0.05, log2 fold-change > 1, and average TPM in mat samples > 1. For each gene annotation group, differential expression data is organized, revealing the number of genes in the annotation group that are significantly upregulated in zoospores and the number that are significantly upregulated in mats. These data are necessary for the Fisher's Exact Statistical Tests (see below). This script generates Fisher's Exact Statistical Test data for gene annotation groups of sufficient size (greater than 10 genes). The Fisher's Exact Test generates 2 adjusted p-values for each annotation group to answer two questions: (1) "Is this gene annotation group upregulated in zoospores?" and (2) "Is this gene annotation group significantly upregulated in mats?" Additional annotations are included: secondary metabolites, Orthologous groups with other anaerobic gut fungal strains using results from the OrthoFinder tool, and CAZymes. |
| Zoospore_RNA_data_analysis_pt4_DGE_Additional_Annotations.py | This script takes the results from part 3 and adds Excel sheet outputs to focus on additional specific annotations or format specific annotations introduced in part 3. |
| Zoospore_RNA_data_analysis_pt5_dbCAN2_and_Cellulosomes.py | In order to roughly describe how cellulosome components are differentially regulated in this transcriptomic dataset, proteomics data with previously detected, likely cellulosome components (dockerin or scaffoldin). The tool dbCAN2 predicts additional CAZyme annotations to supplement CAZyme annotations from JGI MycoCosm. |
| Zoospore_RNA_data_analysis_pt6_Volcano_Plots.py | Additional script for generating volcano plots in the publication. |

### Script Details
*Part 1: RNA Data Analysis RNAseq Cleanup*
Zoospore_RNA_data_analysis_pt1_RNAseq_Cleanup.py

Required Inputs by Script:
1) 'Neosp1_FilteredModels5_deflines_pre.fasta'

FASTA file for genes in Neocallimastix californiae G1 genome,'Neosp1_FilteredModels5_deflines_pre.fasta' downloaded from JGI MycoCosm, under Annotation --> Filtered Models ("best") --> Genes --> "Neosp1_FilteredModels5_deflines.gff.gz" for Project 1029446, 20171122.
The FASTA lists gene keys, identified by MycoCosm proteinID, followed by amino acid sequences

2) 'counts_RNAseq_original.txt'

Note, raw RNA-Seq reads are available at NCBI BioProject: PRJNA982907 to PRJNA982924 

3) 'tpm_counts_RNAseq_original.txt'

Transcripts-per-million normalization is considered in this study. 

Outputs to temp_folder:

1) 'Neosp1_FilteredModels5_deflines_post.fasta'

Modified input FASTA files

2) 'counts_RNAseq_updated.csv'

Counts combined for proteinIDs with identical amino acid sequences.

3) 'tpm_counts_RNAseq_updated.csv'

Transcripts-per-million normalized counts with counts combined for proteinIDs with identical amino acid sequences.

4) 'Neosp1_FilteredModels5_deflines_duped_proteinIDs_sorted.fasta'

This output is a dictionary describing the proteinIDs with identical amino acid sequences.

""""""""""""""""""""""""""""""""""""""""""""""""""""""

*Part 2: DESeq2 Analysis in R*
R v4.2.2
Bioconductor v3.16
DESeq2 v1.36.0

Manually Added Inputs:
1) 'coldata_DESeq2.txt'

DESeq2 Metadata input, describing samples and sample groups

Inputs from Temp Folder (previous script outputs):
1) 'counts_RNAseq_updated.csv'

Outputs to temp_folder:
1) 'deseq2_output.csv'

DESeq2 output, including log2fold-change and p-value data for each proteinID.

2) 'deseq2_normalized_counts.csv'

3) 'deseq2_normalized_counts_labeled.csv'

DESeq2-normazlied counts data with proteinIDs

""""""""""""""""""""""""""""""""""""""""""""""""""""""

*Part 3: Differential Gene Expression Analysis*

Changeable values:
1) mat_samples = list of mat sample names

2) zoosp_samples = list of zoospore sample names

3) log2FC_cutoff = log2 fold-change cutoff for describing significant differential regulation (set value to 1)

4) pval_cutoff = DESeq2 p-value cutoff describing significant differential regulation (set value to 0.05)

5) tpm_cutoff = Transcripts-per-million normalized counts cutoff for describing if a gene is "expressed" or "not expressed." ProteinIDs with TPM values below this cutoff may still be described as expressed, but I implement this cutoff to roughly estimate how many genes for any given gene annotation are likely expressed versus lowly/not expressed. (set value to 1)

Manually Added Inputs:
1) Downloaded JGI MycoCosm Annotations for KOG, GO, IPR, KEGG (file versions generated in 9/18/2017)
KOG_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_KOG.tsv"
GO_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_GO.tsv"
IPR_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_IPR.tsv"
KEGG_annot_filename = "Neosp1_GeneCatalog_proteins_20170918_KEGG.tsv"

2) Secondary Metabolite Annotations from JGI MycoCosm (https://github.com/lbutkovich/Mycocosm_SM_Annotations)
SM_annot_filename = "Neosp1_SMs_orthologs.xlsx"

3) OrthoFinder data from OrthoFinder tool (https://github.com/davidemms/OrthoFinder):
Orthologs_filename = "Orthogroups_GF.tsv"

4) CAZyme annotations from MycoCosm with dbCAN2 predictions considered (https://doi.org/10.1093/nar/gky418):
CAZyme_annot_filename = "G1_cazymes_with_dbCAN2.csv"
Cellulosome-associated proteinIDs:
cellulosome_annot_filename = "G1_cellulosomes_proteomics.csv"

Inputs from Previous Scripts (deposited in temp folder):
1) DESeq2 Output: "deseq2_output.csv" 

2) DESeq2 Normalized counts: "deseq2_normalized_counts_labeled.csv"

3) RNA-Seq counts (cleaned up): "counts_RNAseq_updated.csv"

4) RNA-Seq TPM counts (cleaned up): "tpm_counts_RNAseq_updated.csv" 

Outputs to Temp Folder:
1) 'DGE_Summary_output_main.xlsx'

(see description under version in Output Folder)

2) 'All_annotation_dataframes.xlsx'

In this script, pandas dataframes are generated for the various main annotations. For any given annotation in a main annotation class, the excel data sheets list the proteinIDs in the annotation and the number of those proteinIDs that are classified as "upregulated in zoospores" and "upregulated in mats." 

Outputs to Output Folder:
1) 'volcano.png'

Volcano plot using (log2FC vs. -log10(padj))

2) 'DGE_Statistics_and_Data_Tables.xlsx'

Excel output for general statistics for the transcriptomic dataset

3) 'DGE_Summary_output_main.xlsx'

DGE for main annotations (KOG, GO, IPR, KEGG), including Fisher's Exact Test statistics for gene annotations with 10< genes. This file is stored in both temp and output folder, and downstream scripts take the version from the temp folder and add onto it.

4) 'DGE_Summary_CAZymes_output.xlsx'

This output is similar to "DGE_Summary_output_main.xlsx" except that it covers CAZyme annotations from JGI MycoCosm.

""""""""""""""""""""""""""""""""""""""""""""""""""""""

*Part 4: RNA Data Analysis for Specific Annotations*

Manually Added Inputs:
1) 'Neosp1_SMs_orthologs.xlsx'

also input for part 3

2) 'G1_hydrogenosomes.xlsx'

This excel contains a list of proteinIDs with putative hydrogenosome-related function based on JGI MycoCosm annotations. The excel also contains DeepLoc (https://doi.org/10.1093/bioinformatics/btx431) cellular localization notes, including an amino acid sequence-based prediction score for how a protein localizes to the hydrogenosome ("Mitochondria" and "Plastid" proxy for hydrogenosomes).

3) 'G1_SWEETs.xlsx'

This excel describes Sugars Will Eventually be Exported Transporters annotations for N. californiae (https://doi.org/10.1016/j.ymben.2021.04.009).

4) 'G1_transcription_factors.xlsx'

This excel contains putative annotations for N. californiae transcription factors, based on JGI MycoCosm annotations with either GO ID 3700 ("DNA-binding transcription factor activity") or PF00096 ("Zinc finger, C2H2 type") annotations, and this selection is further narrowed down by considering only genes that contain the strings “transcription” or “Zn-finger” in the KOG defline or KOG class descriptors. 

5) 'G1_UPR_HSR.csv'

This excel contains a list of proteinIDs with unfolded protein response genes and heat shock response genes from S. Seppala et al. (https://doi.org/10.1186/s12934-016-0611-7). These annotations are considered but not focused on in the manuscript.

6) TCDB (Transporter Classification Database, https://doi.org/10.1093/nar/gkaa1004) related files
- 'output_TCDB_BLASTp_filtered.csv' This .csv contains information from local BLASTp (loose cutoff for matches, e-value < 0.01)of N. californiae genes to the genes in the TCDB. The purpose of the local BLASTp is to assign putative transporter functions to N. californiae proteinIDs.
- 'TC_specific_family_defs.csv'
- 'TC_ChEBI_IDs.csv'
- 'TC_superfamily_defs.csv'

Inputs from Previous Scripts (deposited in temp folder):
1) 'DGE_Summary_output_main.xlsx'

from script Part 3

Outputs to Output Folder:
1) 'DGE_summary_output_formatted.xlsx'

Excel DGE analysis for genes with specific annotations per excel sheet (see input annotations).

""""""""""""""""""""""""""""""""""""""""""""""""""""""

Part 5: RNA Data Analysis for CAZyme and Cellulosome Annotations

Manually Added Inputs:
1) dbCAN2 CAZyme predictions (multiple parts, with filenames as 'G1_partX_dbCAN_output.xlsx', where X=1 through 6)

These excel sheets were generated by running the CAZyme prediction tool dbCAN2 (https://doi.org/10.1093/nar/gky418).

Inputs from Previous Scripts (deposited in temp folder):
1) 'DGE_Summary_output_main.xlsx'

2) 'cellulosomes_BLASTp_results.csv'

(see description for 'G1_cellulosomes_proteomics.csv')

3) 'G1_cellulosomes_proteomics.csv'

This .csv describes proteomics data results for N. californiae grown in Medium C with reed canary grass. The proteomics data was originally aligned to one set of N. californiae proteinIDs and filtered for proteinIDs with putative cellulosome component function. In order to match with the rest of this data analysis, I performed local BLASTp of the proteinIDs I am considering (listed in 'Neosp1_FilteredModels5_deflines_post.fasta') against the proteinIDs with proteomics data. The resulting analysis aligns the transcriptomic data (zoospore vs. mats) with proteomics data with putative cellulosome component annotations.

Outputs to Output Folder:
1) 'DGE_summary_dbCAN2_and_cellulosomes.xlsx'

""""""""""""""""""""""""""""""""""""""""""""""""""""""

Part 6: Volcano Plot Generation for Publication

(to be filled in)

## Support
For support with using these scripts, please contact lbutkovich@ucsb.edu.

## Authors and Acknowledgements
Primary author: Lazarina Butkovich (University of California, Santa Barbara)

Thank you to Fred Krauss for feedback and assistance in writing and formatting these scripts. 

Additionally, thank you to Dr. Patrick Leggieri for assistance in determining appropriate statistical analysis. 