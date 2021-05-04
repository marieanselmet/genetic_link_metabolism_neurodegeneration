# To clear the RStudio console: cat("\f")

#########################################################################################################

################################# LOAD DATA #################################
# Load libraries
library(DESeq2)
library(vsn)
library(ggplot2)
library(NMF)
library(edgeR)


# Set path
setwd("~/Desktop/Lab/git/genetic_link_metabolism_neurodegeneration/single_cell_RNAseq_MB/data")
# Load preprocessed counts data (cognition genes and neuron samples for the moment)
readcounts_cognition_genes <- read.delim("processed_read_counts_cognition.txt", row.names="genes")
# Load all CPM corrected and normalized transcription data
#normalized_data_all_genes <- read.delim("GSE74989_NormalizedCPMcorrectedData.txt", row.names="genes")
# Load preprocessed cognition data
#normalized_data_cognition_genes <- read.delim("cognition_expression_data.txt", row.names="genes")


# Make a data.frame with meta-data where row.names should match the individual sample names
sample_info <- data.frame(condition = gsub("_[0-9]+", "", names(readcounts_cognition_genes)), 
                          row.names = names(readcounts_cognition_genes) )
names(sample_info) <- "condition" # rename column as "condition"


################################# DATA NORMALIZATION #################################
# Generate the DESeq DataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts_cognition_genes, colData = sample_info,
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
#### -> We can see here that the log2 distribution representation (transformed counts) is 
####   far more simpler to visualize

# Many statistical tests and analyses assume that data is homoskedastic -> check data homoskedasticity:
msd_plot <- meanSdPlot(log.norm.counts,
                       ranks=FALSE, # show the data on the original scale plot = FALSE)
                       plot = FALSE)
msd_plot$gg + ggtitle("sequencing depth normalized log2(read counts)") + ylab("standard deviation")
#### -> Here there is a dependence of the variance on the mean, which violates the assumption of homoskedasticity.

# Variance shrinkage to reduce the amount of heteroskedasticity
# Obtain regularized log-transformed values: rlog() function returns values that are both normalized for sequencing depth
# and transformed to the log2 scale where the values are adjusted to fit the experiment-wide trend of the variance-mean
# relation-ship
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)

# Mean-sd plot for rlog-transformed data
msd_plot <- meanSdPlot(rlog.norm.counts, ranks=FALSE, # show the data on the original scale plot = FALSE)
                       plot = FALSE)
msd_plot$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation")
#### even worse results?????


##### Determine whether samples display greater variability between experimental conditions than between replicates
# Unsupervised hierarchical clustering to check the similarity between replicates
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson" ))
plot( hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
      main = "rlog transformed read counts\ndistance: Pearson correlation")

# PCA
P <- plotPCA(DESeq.rlog)
P <- P + theme_bw() + ggtitle("Rlog transformed counts")
print(P)


################################# DGE Analysis #################################
####### DESeq #######
str(colData(DESeq.ds)$condition) # factor with 7 levels -> 7 cell types

# Set DAL as the first -level -factor
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "DAL")
model_matrix <- model.matrix(~colData(DESeq.ds)$condition)
print(model_matrix)

# DAL vs V2
# Run the DGE analysis
DESeq.ds <- DESeq(DESeq.ds)
resultsNames(DESeq.ds)
# Walt test
DGE.results <- results(DESeq.ds, contrast=c("condition","DAL","V2"), alpha = 0.05)
summary(DGE.results)
print(DGE.results)

# LRT test
#dds <- DESeq(DESeq.ds, test="LRT", reduced=~1)
#res <- results(dds)
#summary(res)
#print(res)

# Put the DESeq results in a data.frame object
table(DGE.results$padj < 0.05)
rownames(subset(DGE.results, padj < 0.05)) # list of significative (after correction) diffentially expressed genes

# Histogram distribution of the obtained p-values
hist(DGE.results$pvalue ,
     col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies of p-values")
# MA plot: genes that pass the significance threshold (adjusted p-value <0.05) are colored in blue
plotMA(DGE.results, alpha = 0.05, main = "DAL vs. V2")

### Heatmaps ###
# Sort the results according to the adjusted p-value
DGE.results.sorted <- DGE.results[order(DGE.results$padj), ]
# Genes under the desired adjusted p-value significance threshold
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
# Extract the normalized read counts for DE genes into a matrix (aheatmap needs a matrix of values)
hm.mat_DGEgenes <- log.norm.counts[DGEgenes , ]
aheatmap(hm.mat_DGEgenes, Rowv = NA, Colv = NA) # heatmap with sorted p-values
# Combine the heatmap with hierarchical clustering
png("rplot.jpg", width = 350, height = "350", units = "px")
aheatmap(hm.mat_DGEgenes,
         Rowv = TRUE, Colv = TRUE, # add dendrograms to rows and columns
         distfun = "euclidean", hclustfun = "average")
dev.off()
# Scale the read counts per gene to emphasize the sample-type-specific differences
aheatmap(hm.mat_DGEgenes ,
         Rowv = TRUE , Colv = TRUE ,
         distfun = "euclidean", hclustfun = "average",
         scale = "row") # values are transformed into distances from the center of the row-specific average: (actual value - mean of the group) / standard deviation



################# 
####### edgeR #######
# edgeR requires a matrix of read counts where the row names = gene IDs and the column names = sample IDs
sample_info.edger <- factor (c(rep("DAL", 10), rep("V2", 10)))
sample_info.edger <- relevel(sample_info.edger, ref = "DAL")
# Convert the matrix to an edgeR object
edgeR.DGElist <- DGEList(counts = readcounts_cognition_genes, group = sample_info) # why an error? same number of columns and rows yet??
