# To clear the RStudio console: cat("\f")

#########################################################################################################

################################# LOAD DATA #################################
# Load libraries
library(DESeq2)
library(vsn)
library(ggplot2)
library(NMF)
library(writexl)
library(data.table)

# Set path
setwd("~/Desktop/Lab/git/genetic_link_metabolism_neurodegeneration/single_cell_RNAseq_MB/data")
# Load pre-processed counts data (cognition genes and neuron samples for the moment)
readcounts_cognition_genes <- read.delim("processed_read_counts_cognition.txt", row.names="genes")

# Make a data.frame with meta-data where row.names should match the individual sample names
sample_info <- data.frame(neuron_type = gsub("_[0-9]+", "", names(readcounts_cognition_genes)), 
                          row.names = names(readcounts_cognition_genes))
sample_info$neuron_type <- as.factor(sample_info$neuron_type) # make the condition a factor
treatment_ <- c(rep("paired", times=5, each=1), rep("unpaired", times=5, each=1), rep("paired", times=5, each=1), rep("unpaired", times=5, each=1), 
                rep("paired", times=5, each=1), rep("unpaired", times=5, each=1), rep("paired", times=6, each=1), rep("unpaired", times=6, each=1),
                rep("paired", times=5, each=1), rep("unpaired", times=5, each=1), rep("paired", times=4, each=1), rep("unpaired", times=4, each=1),
                rep("paired", times=5, each=1), rep("unpaired", times=5, each=1))
sample_info$treatment <- as.factor(treatment_)


################################# DATA NORMALIZATION #################################
# Generate the DESeq DataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts_cognition_genes, colData = sample_info,
                                   design = ~ treatment + neuron_type )

# Remove genes without any counts
DESeq.ds <- DESeq.ds[ rowSums(counts(DESeq.ds)) > 0, ]

# Investigate different library sizes
colSums(counts(DESeq.ds)) # should be the same as colSums(readcounts_cognition_genes)
colSums(readcounts_cognition_genes)

# Calculate the size factors and add them to the data set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# colData() now contains the sizeFactors
colData(DESeq.ds)

# counts() allows to immediately retrieve the normalized read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE)

# Transform size-factor normalized read counts to log2 scale using a pseudocount of 1
log.norm.counts <- log2(counts.sf_normalized + 1)

# First, boxplots of non-transformed read counts (one per sample)
boxplot(counts.sf_normalized, notch = TRUE,
        main = "Untransformed read counts", ylab = "read counts", xlab = "Cells", axes=FALSE)

# Then, boxplots of the log2-transformed read counts
boxplot(log.norm.counts, notch = TRUE,
        main = "log2-transformed read counts",
        ylab = "log2(read counts)", xlab = "Cells", axes=FALSE)
#### -> We can see here that the log2 distribution representation (transformed counts) is 
####   far more simple to visualize

# Many statistical tests and analyses assume that data is homoskedastic -> check data homoskedasticity:
msd_plot <- meanSdPlot(log.norm.counts,
                       ranks=FALSE, # show the data on the original scale plot = FALSE)
                       plot = FALSE)
msd_plot$gg + ggtitle("Sequencing depth normalized log2(read counts)") + ylab("standard deviation")
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
#### -> Little reduction of the scale of the std variation vs mean


##### Determine whether samples display greater variability between different cell types than between replicates
# Unsupervised hierarchical clustering to check the similarity between replicates
distance.m_rlog <- as.dist(1 - cor(rlog.norm.counts, method = "pearson" ))
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts),
      main = "rlog transformed read counts\ndistance: Pearson correlation", xlab = "Cells");

# PCA
P <- plotPCA(DESeq.rlog,intgroup=c("neuron_type","sizeFactor"))
P <- P + theme_bw() + ggtitle("Rlog transformed counts") +  theme(legend.position = "none") 
print(P) ## Améliorer ce plot


################################# DGE Analysis #################################
####### DESeq #######
str(colData(DESeq.ds)$neuron_type) # factor with 7 levels -> 7 cell types

# Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.58 # Corresponds to a fold change of ~1.5 


# 1) First perform a LRT test to check if there is any difference when the neuron type is included or not in the model while
# controling for the treatment (paired and unpaired)
# -> identify any genes that show changes in expression across one of the 7 neuron types based on the neuron type only

LRT_test <- function(data, padj.cutoff, lfc.cutoff, histogram_p_values=FALSE, MA_plot=FALSE, histo_log2FoldChange=FALSE)
{
  ### Perform LRT test 
  dds_LRT <- data
  design(dds_LRT) <- ~ treatment + neuron_type # you want to control for the treatment to only see the effect of neuron_type
  dds_LRT <- DESeq(dds_LRT, test="LRT", full = design(dds_LRT), reduced = ~ treatment) # LRT test (controlling for the treatment)
  res_LRT_viz <- results(dds_LRT, test="LRT", tidy = FALSE) # Result table for vizualisation
  res_LRT <- results(dds_LRT, test="LRT", tidy = TRUE) # Create an other result table with gene names as a column for saving results in an xlsx file (tidy=TRUE)
  colnames(res_LRT)[1] <- "gene" # Gene names in the 1st column
  
  res_LRT_viz$threshold <- as.logical(res_LRT_viz$padj < padj.cutoff) # Column to denote which genes are significant
  res_LRT$threshold <- as.logical(res_LRT$padj < padj.cutoff) 
  res_LRT_df <- data.frame(res_LRT)
  res_LRT_df <-  res_LRT_df[complete.cases(res_LRT_df), ] # Remove NA values
  res_LRT_thresh <- subset(res_LRT_df, threshold == TRUE) # Remove non significant genes
  
  ### Visualize results
  # Histogram distribution of the obtained p-values
  if (histogram_p_values) {
    pdf("DESeq_results/LRT/histogram_p_values.pdf")
    par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
    hist(res_LRT_viz$pvalue ,
         col = "blue", border = "white", xlab = "", ylab = "", main = "Frequencies of p-values")
    dev.off()
  }
  
  # MA plot: genes that pass the significance threshold (adjusted p-value <0.05) are colored in blue
  if (MA_plot) {
    pdf("DESeq_results/LRT/MA_plot.pdf")
    par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
    plotMA(res_LRT_viz, alpha = 0.05, main = "Full model vs. reduced model")  # shrink the log2 fold changes? -> Seems not necessary for the moment
    dev.off()
  }
  
  
  # Boxplot of expression
  if (histo_log2FoldChange) {
    pdf("DESeq_results/LRT/histo_log2FoldChange.pdf")
    par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
    sorted <- res_LRT_df[order(res_LRT_df$log2FoldChange),]
    print(sorted)
    #hist(log2FoldChange, data = sorted,  notch = TRUE,
            #main = "log2FoldChange across genes",
            #ylab = "log2FoldChange", xlab = "Genes", axes=FALSE, names=res_LRT$gene)
    #plot(x=sorted$log2FoldChange, y=sorted$gene, xlim=c(0,150) , ylim=c(-20,20))
    barplot(sorted$log2FoldChange, main="log2FoldChange", horiz=FALSE, names.arg=sorted$gene, col="#69b3a2", las=2, cex.axis=1.5)
    dev.off()
    
  }
  #print(res_LRT)
  ### Save and return results
  write_xlsx(res_LRT_thresh[,1:6],"DESeq_results/LRT/res_LRT.xlsx") # Save results of significant genes in an xlsx file
  return(res_LRT_viz) # Return all genes (significant and non significant) results 
  
}

# Call LRT test function and visualization of DESeq analysis results
res_LRT <- LRT_test(DESeq.ds, padj.cutoff, lfc.cutoff, TRUE, TRUE, TRUE)
print(res_LRT)





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


# 2) Then perform LRT tests to check if there is any difference between one cell type and the mean of all neurons
# -> identify any genes that show change in expression between one specific cell type and the mean of neurons
dds_mean <- DESeq.ds
design(dds_mean) <- ~ 1 + rep + neuron_type
dds_mean <- DESeq(dds_mean, betaPrior = TRUE)
print(resultsNames(dds_mean))

LRT_test_vs_mean <- function(data, contrast_list, padj.cutoff, lfc.cutoff, file_name)
{
  # Compares mean vs one specified neuron type and writes the significant results in an xlsx file
  res_mean <- results(data, contrast=c(0,1,0,0,0,0,0,0))
  print(res_mean)
  res_mean$threshold <- as.logical(res_mean$padj < padj.cutoff) 
  res_mean_df <- data.frame(res_mean)
  res_mean_thresh <- subset(res_mean_df, threshold == TRUE)
  write_xlsx(res_mean_thresh[,1:6], file_name)
  return(res_mean_df)
}


# DAL vs mean:
res_mean_DAL <- LRT_test_vs_mean(dds_mean, c(0,1,0,0,0,0,0,0), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_DAL.xlsx")
#Histogram distribution of the obtained p-values
hist(res_mean_DAL$pvalue ,
     col = "grey", border = "white", xlab = "", ylab = "", main = "frequencies of p-values")
# MA plot: genes that pass the significance threshold (adjusted p-value <0.05) are colored in blue
plotMA(res_mean_DAL, alpha = 0.05, main = "DAL vs. V2")  # shrink the log2 fold changes?

# Boxplot of expression

# V2 vs mean:
res_mean_V2 <- LRT_test_vs_mean(dds_mean, c(0,0,1,0,0,0,0,0), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_V2.xlsx")

# AB_KCs vs mean:
res_mean_AB_KCs <- LRT_test_vs_mean(dds_mean, c(0,0,0,1,0,0,0,0), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_AB_KCs.xlsx")

# G_KCs vs mean:
res_mean_G_KCs <- LRT_test_vs_mean(dds_mean, c(0,0,0,0,1,0,0,0), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_G_KCs.xlsx")

# V3 vs mean:
res_mean_V3 <- LRT_test_vs_mean(dds_mean, c(0,0,0,0,0,1,0,0), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_V3.xlsx")

# R27 vs mean:
res_mean_R27 <- LRT_test_vs_mean(dds_mean, c(0,0,0,0,0,0,1,0), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_R27.xlsx")

# G386 vs mean:
res_mean_G386 <- LRT_test_vs_mean(dds_mean, c(0,0,0,0,0,0,0,1), padj.cutoff, lfc.cutoff, "DESeq_results/res_mean_G386.xlsx")



# 3) Wald test comparisons 
Wald_test <- function(data, compared_type, reference_type, padj.cutoff, lfc.cutoff, file_name)
{
  # Compares mean vs one specified neuron type and writes the significant results in an xlsx file
  dds_Wald <- data
  colData(dds_Wald)$neuron_type <- relevel(colData(dds_Wald)$neuron_type, reference_type)
  design(dds_Wald) <- ~ 1 + rep + neuron_type   ### include rep? -> 3 beta values do not converge in this case
  dds_Wald <- DESeq(dds_Wald) # Wald test
  res_Wald <- results(dds_Wald, contrast=c("neuron_type", compared_type, reference_type), alpha = padj.cutoff, test="Wald")
  print(res_Wald)
  
  res_Wald$threshold <- as.logical(res_Wald$padj < padj.cutoff) 
  res_Wald_df <- data.frame(res_Wald)
  res_Wald_thresh <- subset(res_Wald_df, threshold == TRUE)
  write_xlsx(res_Wald_thresh[,1:6], file_name)
  return(res_Wald_df)
}


# DAL:
res_DAL_V2 <- Wald_test(DESeq.ds, "V2", "DAL", padj.cutoff, lfc.cutoff, "DESeq_results/res_DAL_V2.xlsx")
res_DAL_V3 <- Wald_test(DESeq.ds, "V3", "DAL", padj.cutoff, lfc.cutoff, "DESeq_results/res_DAL_V3.xlsx")
res_DAL_AB_KCs <- Wald_test(DESeq.ds, "AB_KCs", "DAL", padj.cutoff, lfc.cutoff, "DESeq_results/res_DAL_AB_KCs.xlsx")
res_DAL_G_KCs <- Wald_test(DESeq.ds, "G_KCs", "DAL", padj.cutoff, lfc.cutoff, "DESeq_results/res_DAL_G_KCs.xlsx")
res_DAL_R27 <- Wald_test(DESeq.ds, "R27", "DAL", padj.cutoff, lfc.cutoff, "DESeq_results/res_DAL_R27.xlsx")
res_DAL_G386 <- Wald_test(DESeq.ds, "G386", "DAL", padj.cutoff, lfc.cutoff, "DESeq_results/res_DAL_G386.xlsx")

# V2:
res_V2_DAL <- Wald_test(DESeq.ds, "DAL", "V2", padj.cutoff, lfc.cutoff, "DESeq_results/res_V2_DAL.xlsx")
res_V2_V3 <- Wald_test(DESeq.ds, "V3", "V2", padj.cutoff, lfc.cutoff, "DESeq_results/res_V2_V3.xlsx")
res_V2_AB_KCs <- Wald_test(DESeq.ds, "AB_KCs", "V2", padj.cutoff, lfc.cutoff, "DESeq_results/res_V2_AB_KCs.xlsx")
res_V2_G_KCs <- Wald_test(DESeq.ds, "G_KCs", "V2", padj.cutoff, lfc.cutoff, "DESeq_results/res_V2_G_KCs.xlsx")
res_V2_R27 <- Wald_test(DESeq.ds, "R27", "V2", padj.cutoff, lfc.cutoff, "DESeq_results/res_V2_R27.xlsx")
res_V2_G386 <- Wald_test(DESeq.ds, "G386", "V2", padj.cutoff, lfc.cutoff, "DESeq_results/res_V2_G386.xlsx")

# V3:
res_V3_DAL <- Wald_test(DESeq.ds, "DAL", "V3", padj.cutoff, lfc.cutoff, "DESeq_results/res_V3_DAL.xlsx")
res_V3_V2 <- Wald_test(DESeq.ds, "V2", "V3", padj.cutoff, lfc.cutoff, "DESeq_results/res_V3_V2.xlsx")
res_V3_AB_KCs <- Wald_test(DESeq.ds, "AB_KCs", "V3", padj.cutoff, lfc.cutoff, "DESeq_results/res_V3_AB_KCs.xlsx")
res_V3_G_KCs <- Wald_test(DESeq.ds, "G_KCs", "V3", padj.cutoff, lfc.cutoff, "DESeq_results/res_V3_G_KCs.xlsx")
res_V3_R27 <- Wald_test(DESeq.ds, "R27", "V3", padj.cutoff, lfc.cutoff, "DESeq_results/res_V3_R27.xlsx")
res_V3_G386 <- Wald_test(DESeq.ds, "G386", "V3", padj.cutoff, lfc.cutoff, "DESeq_results/res_V3_G386.xlsx")

# AB_KCs:
res_AB_KCs_DAL <- Wald_test(DESeq.ds, "DAL", "AB_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_AB_KCs_DAL.xlsx")
res_AB_KCs_V2 <- Wald_test(DESeq.ds, "V2", "AB_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_AB_KCs_V2.xlsx")
res_AB_KCs_V3 <- Wald_test(DESeq.ds, "V3", "AB_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_AB_KCs_V3.xlsx")
res_AB_KCs_G_KCs <- Wald_test(DESeq.ds, "G_KCs", "AB_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_AB_KCs_G_KCs.xlsx")
res_AB_KCs_R27 <- Wald_test(DESeq.ds, "R27", "AB_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_AB_KCs_R27.xlsx")
res_AB_KCs_G386 <- Wald_test(DESeq.ds, "G386", "AB_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_AB_KCs_G386.xlsx")

# G_KCs:
res_G_KCs_DAL <- Wald_test(DESeq.ds, "DAL", "G_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_G_KCs_DAL.xlsx")
res_G_KCs_V2 <- Wald_test(DESeq.ds, "V2", "G_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_G_KCs_V2.xlsx")
res_G_KCs_V3 <- Wald_test(DESeq.ds, "V3", "G_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_G_KCs_V3.xlsx")
res_G_KCs_AB_KCs <- Wald_test(DESeq.ds, "AB_KCs", "G_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_G_KCs_AB_KCs.xlsx")
res_G_KCs_R27 <- Wald_test(DESeq.ds, "R27", "G_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_G_KCs_R27.xlsx")
res_G_KCs_G386 <- Wald_test(DESeq.ds, "G386", "G_KCs", padj.cutoff, lfc.cutoff, "DESeq_results/res_G_KCs_G386.xlsx")

# R27:
res_R27_DAL <- Wald_test(DESeq.ds, "DAL", "R27", padj.cutoff, lfc.cutoff, "DESeq_results/res_R27_DAL.xlsx")
res_R27_V2 <- Wald_test(DESeq.ds, "V2", "R27", padj.cutoff, lfc.cutoff, "DESeq_results/res_R27_V2.xlsx")
res_R27_V3 <- Wald_test(DESeq.ds, "V3", "R27", padj.cutoff, lfc.cutoff, "DESeq_results/res_R27_V3.xlsx")
res_R27_AB_KCs <- Wald_test(DESeq.ds, "AB_KCs", "R27", padj.cutoff, lfc.cutoff, "DESeq_results/res_R27_AB_KCs.xlsx")
res_R27_G_KCs <- Wald_test(DESeq.ds, "G_KCs", "R27", padj.cutoff, lfc.cutoff, "DESeq_results/res_R27_G_KCs.xlsx")
res_R27_G386 <- Wald_test(DESeq.ds, "G386", "R27", padj.cutoff, lfc.cutoff, "DESeq_results/res_R27_G386.xlsx")

# G386:
res_G386_DAL <- Wald_test(DESeq.ds, "DAL", "G386", padj.cutoff, lfc.cutoff, "DESeq_results/res_G386_DAL.xlsx")
res_G386_V2 <- Wald_test(DESeq.ds, "V2", "G386", padj.cutoff, lfc.cutoff, "DESeq_results/res_G386_V2.xlsx")
res_G386_V3 <- Wald_test(DESeq.ds, "V3", "G386", padj.cutoff, lfc.cutoff, "DESeq_results/res_G386_V3.xlsx")
res_G386_AB_KCs <- Wald_test(DESeq.ds, "AB_KCs", "G386", padj.cutoff, lfc.cutoff, "DESeq_results/res_G386_AB_KCs.xlsx")
res_G386_G_KCs <- Wald_test(DESeq.ds, "G_KCs", "G386", padj.cutoff, lfc.cutoff, "DESeq_results/res_G386_G_KCs.xlsx")
res_G386_R27 <- Wald_test(DESeq.ds, "R27", "G386", padj.cutoff, lfc.cutoff, "DESeq_results/res_G386_R27.xlsx")


############################ Diagnostic plots ############################


# Put the DESeq results in a data.frame object
#table(DGE.results$padj < 0.05)
#rownames(subset(DGE.results, padj < 0.05)) # list of significative (after correction) diffentially expressed genes

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





