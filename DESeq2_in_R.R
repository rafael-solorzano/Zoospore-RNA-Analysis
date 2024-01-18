################################################################################                      
### Zoospore RNA Manuscript: DESEQ2 Analysis in R by Lazarina Butkovich
################################################################################
# https://lashlock.github.io/compbio/R_presentation.html
# Citation for DESeq2: Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

# Helpful website: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
############################################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
############################################################

############################################################
#Necessary libraries
library("DESeq2")
# library('GenomicRanges')
# ##library('rtracklayer')
# require(gplots)
# require(RColorBrewer)
# library(data.table)
############################################################
############################################################
# install.packages("htmltools")
# library(htmltools)
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# 
# library( "DESeq2" )
# library(ggplot2)
############################################################
############################################################

#DESeq2 Script:
# This script is designed for 3 zoospore replicates and 15 mat replicates.
# Double check counts data and coldata have the appropriate columns/rows

#Read counts data from file
#cts = read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\counts.txt",header=TRUE)
cts = read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\counts_2023_2_zoosp_reps.txt", header=TRUE, )
# cts <- as.matrix(read.csv("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\counts.txt"),sep=",",row.names="GeneID")


# Make the first column values be the string after the last '|'
cts$GeneID = gsub(".*\\|", "", cts$GeneID, perl=TRUE)

# Change GeneID column name to proteinID
colnames(cts)[1] <- "proteinID"

# Make the values in the proteinID column be integers
#cts$proteinID <- as.integer(cts$proteinID)

# Make the index values of cts be the values in proteinID column
rownames(cts) <- cts$proteinID

# Remove the proteinID column because # cts columns needs to match # rows coldata
cts <- cts[,-1]

# Import the metadata, or sample information table
coldata = read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\coldata_final_20230126_2_zoosp_reps.txt",header=TRUE)
# coldata = read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\coldata_treated_untreated_20230130.txt",header=TRUE)
rownames(coldata) <- coldata$sample
# condition (and type if applicable) are factor variables, so use factor function on those columns
coldata$condition <- factor(coldata$condition, levels = c('zoosp','mat'))
#coldata$type <- factor(coldata$type)

#Rows of coldata sample names need to match the order of cts column sample names, ignoring proteinID column, so check:
if (all(rownames(coldata) != colnames(cts))) {
  print("The rownames of coldata do not match the columns of cts")
}

dds = DESeqDataSetFromMatrix(countData =cts, colData = coldata, design =~condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, alpha=0.05)
head(results(dds, tidy=TRUE))
summary(res)

# Make the index values into column values for a new column proteinID
res$proteinID <- rownames(res)
res <- res[, c(ncol(res), 1:(ncol(res)-1))]
# Rename 'log2FoldChange' column to 'zoosp_vs_mat_log2FC'
colnames(res)[colnames(res) == 'log2FoldChange'] <- 'mat_vs_zoosp_log2FC'
# Rename 'padj' column to 'zoosp_vs_mat_Padj'
colnames(res)[colnames(res) == 'padj'] <- 'mat_vs_zoosp_Padj'

# Export results file to .csv for data analysis in Python script
write.csv(res,"C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\deseq2_output_20230315_2_zoosp_reps.csv", row.names = FALSE)

# Get DESeq2-normalized counts data used for log2FC calculations. https://support.bioconductor.org/p/66067/ 
DESeq2_cts <- counts(dds,normalized=TRUE)

# Export DESEq2-normalized counts
write.csv(DESeq2_cts,"C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\deseq2_normalized_counts_20230315_2_zoosp_reps.csv", row.names = FALSE)
# Re-import DESeq2-normalized counts and add proteinIDs column info
DESeq2_cts_unlabeled = read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\deseq2_normalized_counts_20230315_2_zoosp_reps.csv", header=TRUE, sep=",")
DESeq2_cts_unlabeled$proteinID <- rownames(res)
DESeq2_cts_labeled <- DESeq2_cts_unlabeled[, c(ncol(DESeq2_cts_unlabeled), 1:(ncol(DESeq2_cts_unlabeled)-1))]
# Make columns for DESeq2-normalized counts averages in zoosp and mats
cols_zoosp = coldata[coldata$condition == 'zoosp', 'sample']
cols_mat = coldata[coldata$condition == 'mat', 'sample']
DESeq2_cts_labeled$zoosp_avg_DESeq2_norm_cts <- apply(DESeq2_cts_labeled[, cols_zoosp], 1, mean)
DESeq2_cts_labeled$mat_avg_DESeq2_norm_cts <- apply(DESeq2_cts_labeled[, cols_mat], 1, mean)
# Calculate variances
DESeq2_cts_labeled$zoosp_var_DESeq2_norm_cts <- apply(DESeq2_cts_labeled[, cols_zoosp], 1, var)
DESeq2_cts_labeled$mat_var_DESeq2_norm_cts <- apply(DESeq2_cts_labeled[, cols_mat], 1, var)
# Calculate log2FC from avg DESeq2-normalized counts (not taking into account shrinkage)
DESeq2_cts_labeled$log2FC_check <- log2(DESeq2_cts_labeled$mat_avg_DESeq2_norm_cts / DESeq2_cts_labeled$zoosp_avg_DESeq2_norm_cts)
# 1/30/23: The log2FC_check values are all mostly slightly off from the DESeq2 reported log2FC

# Re-export DESEq2-normalized counts
write.csv(DESeq2_cts_labeled,"C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\deseq2_normalized_counts_labeled_20230315_2_zoosp_reps.csv", row.names = FALSE)

# Volcano Plot #1: without tpm counts cutoff
#reset par, which describes how many plots will appear at once
par(mfrow=c(1,1))
# Make a basic volcano plot, just labeling |log2FC|>1 and p_adj<0.05
with(res, plot(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, main="Volcano plot"))
with(subset(res, mat_vs_zoosp_Padj<.05 & mat_vs_zoosp_log2FC>1), points(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col="blue"))
with(subset(res, mat_vs_zoosp_Padj<.05 & mat_vs_zoosp_log2FC<(-1)), points(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col="red"))
legend("topright", legend=c("Upregulated in Mats, *p_adj<0.05", "Upregulated in Zoospores, *p_adj<0.05"), pch=20, col=c("blue", "red"), bty="n")

# Manually export if you want to keep the plot


# Volcano Plot #2: with tpm counts cutoff
# Import and generate zoosp and mat tpm average values
tpm_cts = read.delim("C:\\Users\\lazab\\Desktop\\python_scripts\\workspace\\Zoospore_RNA\\tpm_counts_updated_sorted_2_zoosp_reps.csv", header=TRUE, sep=",")
# Make the first column values be the string after the last '|'
tpm_cts$proteinID_str = gsub(".*\\|", "", tpm_cts$proteinID_str)
# Change GeneID column name to proteinID
colnames(tpm_cts)[1] <- "proteinID"
# Make proteinID values numbers
tpm_cts$proteinID = as.numeric(tpm_cts$proteinID)
# Order by proteinID
tpm_cts <- tpm_cts[order(tpm_cts$proteinID),]
# reset index for tpm_cts
rownames(tpm_cts) <- NULL
# create average mat and zoosp tpm_cts columns
tpm_cts$zoosp_avg_tpm <- apply(tpm_cts[, cols_zoosp], 1, mean)
tpm_cts$mat_avg_tpm <- apply(tpm_cts[, cols_mat], 1, mean)

res_easy = read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\deseq2_normalized_counts_labeled_20230315_2_zoosp_reps.csv",header=TRUE,sep=",")
res_easy <- res_easy[order(res_easy$proteinID),]
# reset index for res_easy
rownames(res_easy) <- NULL

# Add tpm_cts data to res_easy dataframe
res_easy$zoosp_avg_tpm <- tpm_cts$zoosp_avg_tpm
res_easy$mat_avg_tpm <- tpm_cts$mat_avg_tpm

# Re-import DESeq2 results (log2FC and padj)
deseq2_output_easy <- read.delim("C:\\Users\\lazab\\Desktop\\R_scripts\\DESeq2 zoospores\\deseq2_output_20230315_2_zoosp_reps.csv",sep=",")
deseq2_output_easy$proteinID = as.numeric(deseq2_output_easy$proteinID)
deseq2_output_easy <- deseq2_output_easy[order(deseq2_output_easy$proteinID),]
rownames(deseq2_output_easy) <- NULL

# Add mat_vs_zoosp_log2FC from deseq2_output_easy to res_easy
res_easy$mat_vs_zoosp_log2FC <- deseq2_output_easy$mat_vs_zoosp_log2FC
# Add mat_vs_zoosp_Padj from deseq2_output_easy to res_easy
res_easy$mat_vs_zoosp_Padj <- deseq2_output_easy$mat_vs_zoosp_Padj

# Make another volcano plot that includes tpm count upreg'd >1 cutoff, using res_easy dataframe
with(res_easy, plot(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, main="Volcano plot"))
with(subset(res_easy, mat_vs_zoosp_Padj<.05 & mat_vs_zoosp_log2FC>1 & mat_avg_tpm>1), points(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col="lightblue3"))
with(subset(res_easy, mat_vs_zoosp_Padj<.05 & mat_vs_zoosp_log2FC<(-1) & zoosp_avg_tpm>1), points(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col="indianred3"))
legend("topright", legend=c("Upregulated in Mats", "Upregulated in Zoospores"), pch=20, col=c("lightblue3", "indianred3"), bty="n")

# For transparent points, use the below code with library scales:
library(scales)
with(res_easy, plot(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col=alpha("gray40", 0.5),xlab="log2-fold change", ylab="-log10(q)",cex.lab=1.25))
with(subset(res_easy, mat_vs_zoosp_Padj<.05 & mat_vs_zoosp_log2FC>1 & mat_avg_tpm>1), points(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col=alpha("lightblue3", 0.5)))
with(subset(res_easy, mat_vs_zoosp_Padj<.05 & mat_vs_zoosp_log2FC<(-1) & zoosp_avg_tpm>1), points(mat_vs_zoosp_log2FC, -log10(mat_vs_zoosp_Padj), pch=20, col=alpha("indianred3", 0.5)))
legend("topright", legend=c("Upregulated in Mats", "Upregulated in Zoospores"), cex=1.25, pch=20, col=c("lightblue3", "indianred3"), bty="n")

# Manually export if you want to keep the plot