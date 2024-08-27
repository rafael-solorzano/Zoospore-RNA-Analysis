############################################################

# Zoospore RNA Manuscript: DESeq2 Analysis in R by Lazarina Butkovich

# https://lashlock.github.io/compbio/R_presentation.html
# Citation for DESeq2: Love, M.I., Huber, W., Anders, S. (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15:550. 10.1186/s13059-014-0550-8

# Helpful website: http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start


############################################################
# Install DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library("DESeq2")

############################################################


###
#Import files
###
# Indicate input folder with manually prepared inputs
setwd('input')
# Import the metadata, or sample information table
coldata = read.delim("coldata_DESeq2.txt",header=TRUE)

# Indicate temp folder with inputs from previous script
setwd('../')
setwd('temp')
# Read counts data from file
cts = read.delim("counts_RNAseq_updated.csv", header=TRUE, sep=",")


###
#Format counts and metadata
###
# Format cts data
# Make the first column values be the string after the last '|'
cts$GeneID = gsub(".*\\|", "", cts$GeneID, perl=TRUE)
# Change GeneID column name to proteinID
colnames(cts)[1] <- "proteinID"
# Make the index values of cts be the values in proteinID column
rownames(cts) <- cts$proteinID
# Remove the proteinID column because # cts columns needs to match # rows coldata
cts <- cts[,-1]

# Format coldata metadata
rownames(coldata) <- coldata$sample
# condition (and type if applicable) are factor variables, so use factor function on those columns
coldata$condition <- factor(coldata$condition, levels = c('zoosp','mat'))
#Rows of coldata sample names need to match the order of cts column sample names, ignoring proteinID column, so check:
if (all(rownames(coldata) != colnames(cts))) {
  print("The rownames of coldata do not match the columns of cts")
}


###
#Run DESeq2
###
dds = DESeqDataSetFromMatrix(countData =cts, colData = coldata, design =~condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, alpha=0.05)
head(results(dds, tidy=TRUE))
summary(res)


###
#Format DESeq2 results
###
# Make the index values into column values for a new column proteinID
res$proteinID <- rownames(res)
res <- res[, c(ncol(res), 1:(ncol(res)-1))]
# Rename 'log2FoldChange' column to 'mat_vs_zoosp_log2FC'
colnames(res)[colnames(res) == 'log2FoldChange'] <- 'mat_vs_zoosp_log2FC'
# Rename 'padj' column to 'zoosp_vs_mat_Padj'
colnames(res)[colnames(res) == 'padj'] <- 'mat_vs_zoosp_Padj'


###
#Export results
###
# Deposit results in temp folder
# Export results file to .csv for data analysis in Python script
write.csv(res,"deseq2_output.csv", row.names = FALSE)

# Get DESeq2-normalized counts data used for log2FC calculations. https://support.bioconductor.org/p/66067/ 
DESeq2_cts <- counts(dds,normalized=TRUE)

# Export DESEq2-normalized counts
write.csv(DESeq2_cts,"deseq2_normalized_counts.csv", row.names = FALSE)

# Re-import DESeq2-normalized counts and add proteinIDs column info
DESeq2_cts_unlabeled = read.delim("deseq2_normalized_counts.csv", header=TRUE, sep=",")
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

# Re-export formatted DESEq2-normalized counts
write.csv(DESeq2_cts_labeled,"deseq2_normalized_counts_labeled.csv", row.names = FALSE)