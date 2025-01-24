# Load required library
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Step 1: Read the count matrix
# Assuming "reformatted_counts.txt" is the reformatted FeatureCounts output
counts <- read.table("C:\\Sem 3\\RNA _seq\\count\\reformatted_counts.txt", header = TRUE, row.names = 1)

# Step 2: Define metadata (colData)
# Update the 'condition' column as per your experimental design
colData <- data.frame(
  row.names = colnames(counts),
  condition = factor(c(
    "HER2", "HER2", "HER2", 
    "NonTNBC", "NonTNBC", "NonTNBC", 
    "Normal", "Normal", "Normal", 
    "TNBC", "TNBC", "TNBC"
  ))
)

# Step 3: Create the DESeqDataSet object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ condition
)

# Step 4: Run DESeq2 analysis
dds <- DESeq(dds)

# Step 5: Variance Stabilizing Transformation (VST)
# We set blind=TRUE as we are focusing on visualization
vsd <- vst(dds, blind = TRUE)

# Step 6: Assess clustering using PCA
# PCA plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Custom PCA Plot using ggplot2
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA of Samples Based on Gene Expression")

# Step 7: Heatmap of sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colData$condition
colnames(sampleDistMatrix) <- colData$condition

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         main = "Sample Distance Heatmap")

# Save outputs (optional)
# Save PCA data
write.csv(pcaData, "PCA_data.csv")
# Save vst-transformed data
write.csv(as.data.frame(assay(vsd)), "vst_transformed_data.csv")

