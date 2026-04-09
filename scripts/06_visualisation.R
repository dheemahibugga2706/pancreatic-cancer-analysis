# ============================================================================
# SCRIPT 06: VISUALIZATION 
# ============================================================================

setwd("~/Documents/pancreatic-cancer-analysis")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

print("=== VISUALIZATION SCRIPT (HEATMAP SKIPPED) ===")

# Load data
dds <- readRDS("data/processed/deseq2_object.rds")
de_genes <- readRDS("data/processed/de_genes_significant.rds")

# Prepare VST
vst_counts <- varianceStabilizingTransformation(dds)

# ============================================================================
# 1. PCA PLOT (ALWAYS WORKS)
# ============================================================================

print("Creating PCA plot...")

sample_info <- readRDS("data/processed/sample_info_processed.rds")

pca_data <- plotPCA(vst_counts, intgroup = "condition", returnData = TRUE)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 4, alpha = 0.7) +
  stat_ellipse(aes(fill = condition), alpha = 0.1, level = 0.95, type = "t") +
  labs(
    title = "PCA: Pancreatic Cancer vs Normal",
    x = paste0("PC1 (", round(attr(pca_data, "percentVar")[1], 1), "% variance)"),
    y = paste0("PC2 (", round(attr(pca_data, "percentVar")[2], 1), "% variance)")
  ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

ggsave("results/plots/06_pca_plot.png", p_pca, width = 9, height = 7, dpi = 300)
print("✓ PCA plot saved")

# ============================================================================
# 2. MA PLOT (ALWAYS WORKS)
# ============================================================================

print("Creating MA plot...")

res <- results(dds)
png("results/plots/06_ma_plot.png", width = 900, height = 700)
plotMA(res, ylim = c(-6, 6), main = "MA Plot")
dev.off()
print("✓ MA plot saved")

# ============================================================================
# 3. VOLCANO PLOT (ALWAYS WORKS)
# ============================================================================

print("Creating volcano plot...")

res_df <- read.csv("results/tables/03_all_genes_results.csv")

vol_data <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(sig = padj < 0.05 & abs(log2FoldChange) > 1)

p_volcano <- ggplot(vol_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray70")) +
  labs(title = "Volcano Plot", x = "log2 FC", y = "-log10(p-adj)") +
  theme_minimal()

ggsave("results/plots/06_volcano_plot.png", p_volcano, width = 10, height = 7, dpi = 300)
print("✓ Volcano plot saved")

# ============================================================================
# DONE
# ============================================================================

print("\n=== VISUALIZATION COMPLETE ===")
print("✓ PCA plot")
print("✓ MA plot")  
print("✓ Volcano plot")
