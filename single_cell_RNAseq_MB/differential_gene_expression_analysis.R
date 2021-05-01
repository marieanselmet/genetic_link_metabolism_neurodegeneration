# To clear the RStudio console: cat("\f")

#########################################################################################################

################################# LOAD DATA #################################
library(DESeq2)
# Set path
setwd("~/Desktop/Lab/git/genetic_link_metabolism_neurodegeneration/single_cell_RNAseq_MB/data")
# Load preprocessed counts data (all genes and all samples)
readcounts_all_genes <- read.delim("processed_read_counts_cognition.txt", row.names="genes")
# Load all CPM corrected and normalized transcription data
normalized_data_all_genes <- read.delim("GSE74989_NormalizedCPMcorrectedData.txt", row.names="genes")
# Load preprocessed cognition data
normalized_data_cognition_genes <- read.delim("cognition_expression_data.txt", row.names="genes")


# Make a data.frame with meta-data where row.names should match the individual sample names
sample_info <- data.frame(condition = gsub("_[0-9]+", "", names(readcounts_all_genes)), 
                          row.names = names(readcounts_all_genes) )
names(sample_info) <- "condition" # rename column as "condition"


################################# DATA NORMALIZATION #################################
# Generate the DESeq DataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts_all_genes, colData = sample_info,
                                   design = ~ condition)
# Check the DESeq object:
colData(DESeq.ds) %>% head
assay(DESeq.ds, "counts") %>% head
rowData(DESeq.ds) %>% head

# Test what counts () returns
counts(DESeq.ds) %>% str

# Remove genes without any counts
DESeq.ds <- DESeq.ds[ rowSums(counts(DESeq.ds)) > 0, ]

# Investigate different library sizes
colSums(counts(DESeq.ds)) # should be the same as colSums(readcounts)

# Calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# If you check colData() again, you see that this now contains the sizeFactors
colData(DESeq.ds)

# counts() allows you to immediately retrieve the _normalized_ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE)

# Transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)


par(mfrow=c(2,1)) # to plot the following two images underneath each other
# First, boxplots of non-transformed read counts (one per sample)
boxplot(counts.sf_normalized, notch = TRUE,
        main = "untransformed read counts", ylab = "read counts")
# Box plots of log2-transformed read counts
boxplot(log.norm.counts, notch = TRUE,
        main = "log2-transformed read counts",
        ylab = "log2(read counts)")



str(colData(DESeq.ds)$condition)
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition , "DAL")
DESeq.ds <- DESeq(DESeq.ds)
DGE.results <- results(DESeq.ds, independentFiltering = TRUE, alpha = 0.05)
summary(DGE.results)
head(DGE.results)
