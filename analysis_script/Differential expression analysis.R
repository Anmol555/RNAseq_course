# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)

setwd("C:/Sem 3/RNA _seq/count")

# Step 1: Check DESeq2 setup
if (!exists("dds")) stop("The DESeqDataSet object 'dds' does not exist. Ensure 'dds' is created and DESeq() has been run.")

# Define a function to extract results with error handling
extract_results <- function(dds, group1, group2) {
  contrast_name <- paste(group1, "vs", group2)
  res <- tryCatch(
    {
      results(dds, contrast = c("condition", group1, group2))
    },
    error = function(e) {
      cat(paste("Error extracting results for", contrast_name, ":", e$message, "\n"))
      return(NULL)
    }
  )
  return(res)
}

# Step 2: Extract pairwise contrasts
contrasts <- list(
  c("HER2", "Normal"),
  c("HER2", "TNBC"),
  c("HER2", "NonTNBC"),
  c("Normal", "TNBC"),
  c("Normal", "NonTNBC"),
  c("TNBC", "NonTNBC")
)

results_list <- lapply(contrasts, function(contrast) {
  extract_results(dds, contrast[1], contrast[2])
})
names(results_list) <- sapply(contrasts, function(contrast) paste(contrast[1], "vs", contrast[2]))

# Step 3: Filter significant DE genes
filter_results <- function(res) {
  if (is.null(res)) return(NULL)
  res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
}

filtered_results <- lapply(results_list, filter_results)

# Step 4: Summarize results
summarize_results <- function(res, contrast_name) {
  if (is.null(res)) {
    cat(paste("No results for", contrast_name, "\n"))
    return()
  }
  upregulated <- sum(res$log2FoldChange > 0, na.rm = TRUE)
  downregulated <- sum(res$log2FoldChange < 0, na.rm = TRUE)
  cat(paste0(contrast_name, ": Upregulated = ", upregulated, ", Downregulated = ", downregulated, "\n"))
}

cat("Summary of Differentially Expressed Genes:\n")
invisible(lapply(names(filtered_results), function(name) {
  summarize_results(filtered_results[[name]], name)
}))

# Step 5: Save results to CSV
save_results <- function(res, contrast_name) {
  if (is.null(res)) return()
  filename <- paste0("sig_", gsub(" ", "_", contrast_name), ".csv")
  write.csv(as.data.frame(res), filename)
  cat(paste("Results saved for:", contrast_name, "as", filename, "\n"))
}

invisible(lapply(names(filtered_results), function(name) {
  save_results(filtered_results[[name]], name)
}))

# Step 6: Generate MA Plots
generate_MA_plot <- function(res, contrast_name) {
  if (is.null(res)) return()
  png_filename <- paste0("MA_plot_", gsub(" ", "_", contrast_name), ".png")
  png(png_filename)
  plotMA(res, main = paste("MA Plot:", contrast_name), ylim = c(-5, 5))
  dev.off()
  cat(paste("MA Plot saved for:", contrast_name, "as", png_filename, "\n"))
}

invisible(lapply(names(results_list), function(name) {
  generate_MA_plot(results_list[[name]], name)
}))

# Step 7: Generate Improved Volcano Plots
generate_improved_volcano_plot <- function(res, contrast_name, fold_change_threshold = 1, pvalue_threshold = 0.05) {
  if (is.null(res)) {
    cat(paste("No results available for", contrast_name, "\n"))
    return()
  }
  
  # Add columns for -log10(p-value) and significance classification
  res$logp <- -log10(res$padj)
  res$significance <- ifelse(res$padj < pvalue_threshold & abs(res$log2FoldChange) > fold_change_threshold,
                             ifelse(res$log2FoldChange > 0, "Upregulated", "Downregulated"),
                             "Not Significant")
  
  # Add gene labels for only significant genes
  res$gene_label <- ifelse(res$significance != "Not Significant", rownames(res), NA)
  
  # Create volcano plot
  volcano_plot <- ggplot(res, aes(x = log2FoldChange, y = logp, color = significance)) +
    geom_point(alpha = 0.7, size = 2.5) +  # Enhanced visibility for points
    geom_text_repel(
      aes(label = gene_label),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.35,
      point.padding = 0.3,
      segment.color = "black"  # Lines connecting labels to points
    ) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = "white"),
      plot.background = element_rect(fill = "white", color = "white"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 12),  # Improved axis text visibility
      axis.title = element_text(color = "black", size = 14, face = "bold"),
      legend.title = element_text(color = "black", size = 12),
      legend.text = element_text(color = "black", size = 10),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "black")
    ) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted p-value") +
    ggtitle(paste("Enhanced Volcano Plot:", contrast_name)) +
    geom_hline(yintercept = -log10(pvalue_threshold), linetype = "dashed", color = "black", size = 0.8) +
    geom_vline(xintercept = c(-fold_change_threshold, fold_change_threshold), linetype = "dashed", color = "black", size = 0.8) +
    coord_cartesian(xlim = c(-max(abs(res$log2FoldChange), na.rm = TRUE), max(abs(res$log2FoldChange), na.rm = TRUE)),
                    ylim = c(0, max(res$logp, na.rm = TRUE) + 5))
  
  # Save plot as PNG
  filename <- paste0("Improved_Enhanced_Volcano_plot_", gsub(" ", "_", contrast_name), ".png")
  ggsave(filename = filename, plot = volcano_plot, width = 12, height = 8)
  cat(paste("Improved Volcano Plot saved for:", contrast_name, "as", filename, "\n"))
}

# Generate improved volcano plots for all contrasts
invisible(lapply(names(results_list), function(name) {
  generate_improved_volcano_plot(results_list[[name]], name)
}))
