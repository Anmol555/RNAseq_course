# Load required libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # For human annotations (use org.Mm.eg.db for mouse, etc.)
library(DESeq2)
library(ggplot2)

setwd("C:/Sem 3/RNA _seq/count")

# Define pairwise comparisons
comparisons <- list(
  "HER2_vs_Normal" = c("HER2", "Normal"),
  "TNBC_vs_Normal" = c("TNBC", "Normal"),
  "NonTNBC_vs_Normal" = c("NonTNBC", "Normal")
)

# Loop through each comparison
for (comparison_name in names(comparisons)) {
  # Extract comparison details
  contrast <- comparisons[[comparison_name]]
  numerator <- contrast[1]
  denominator <- contrast[2]
  
  # Step 1: Define Differentially Expressed Genes (DEGs)
  # Extract results for the current pairwise comparison
  res <- results(dds, contrast = c("condition", numerator, denominator))
  
  # Filter DE genes (padj < 0.05)
  DE_genes <- res[which(res$padj < 0.05), ]
  
  # Extract Ensembl IDs of DE genes
  DE_genes <- rownames(DE_genes)
  
  # Step 2: Define the "universe" (all genes measured in your analysis)
  # Extract Ensembl IDs of all genes from your dataset
  universe_genes <- rownames(res)
  
  # Step 3: Perform GO Enrichment Analysis
  go_enrichment <- enrichGO(
    gene = DE_genes,                # Ensembl IDs of DE genes
    universe = universe_genes,      # Ensembl IDs of all genes measured
    OrgDb = org.Hs.eg.db,           # Annotation package for the organism (human)
    keyType = "ENSEMBL",            # Input gene IDs are in Ensembl format
    ont = "BP",                     # Ontology: Biological Process (can also use "MF", "CC", or "ALL")
    pAdjustMethod = "BH",           # Adjust p-values using Benjamini-Hochberg
    pvalueCutoff = 0.05,            # Adjusted p-value cutoff for significance
    qvalueCutoff = 0.2,             # Q-value cutoff for significance
    readable = TRUE                 # Convert Ensembl IDs to gene symbols in results
  )
  
  # Save enrichment results to CSV
  csv_path <- paste0("GO_enrichment_results_", comparison_name, ".csv")
  write.csv(as.data.frame(go_enrichment), csv_path)
  print(paste("GO enrichment results saved to:", csv_path))
  
  # Step 4: Visualize the Results
  # Dot plot
  dotplot_path <- paste0("GO_enrichment_dotplot_", comparison_name, ".png")
  dot_plot <- dotplot(go_enrichment, showCategory = 20) +
    ggtitle(paste("GO Enrichment: Biological Process\n", comparison_name))
  ggsave(dotplot_path, plot = dot_plot, width = 10, height = 6)
  print(paste("Dot plot saved to:", dotplot_path))
  
  # Bar plot
  barplot_path <- paste0("GO_enrichment_barplot_", comparison_name, ".png")
  bar_plot <- barplot(go_enrichment, showCategory = 20) +
    ggtitle(paste("GO Enrichment: Biological Process\n", comparison_name))
  ggsave(barplot_path, plot = bar_plot, width = 10, height = 6)
  print(paste("Bar plot saved to:", barplot_path))
}

