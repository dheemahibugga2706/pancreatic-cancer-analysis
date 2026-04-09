

# Load your expression data
expr_matrix <- readRDS("data/processed/expression_filtered.rds")

# Check first 10 gene IDs
print("First 10 gene IDs:")
print(head(rownames(expr_matrix), 10))

# Check format
print("Gene ID examples:")
print(rownames(expr_matrix)[1:5])
# ============================================================================
# SCRIPT 03: DIFFERENTIAL EXPRESSION ANALYSIS (CORRECTED)
# Find genes significantly different between Cancer and Normal
# ============================================================================

setwd("~/Documents/pancreatic-cancer-analysis")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(org.Hs.eg.db)

# Load DESeq2 object
dds <- readRDS("data/processed/deseq2_object.rds")

print("DESeq2 object loaded")

# ============================================================================
# 1. CHECK GENE ID FORMAT
# ============================================================================

print("=== CHECKING GENE ID FORMAT ===")
print("First 5 gene IDs:")
print(head(rownames(dds), 5))

# Detect gene ID type
gene_ids <- rownames(dds)

# Check if ENSEMBL
is_ensembl <- all(grepl("^ENSG", gene_ids[1:min(100, length(gene_ids))]))

# Check if probe IDs (contains _)
is_probe <- all(grepl("_", gene_ids[1:min(100, length(gene_ids))]))

# Check if looks like gene symbols (short, alphanumeric)
is_symbol <- all(nchar(gene_ids[1:min(100, length(gene_ids))]) < 20 & 
                   !grepl("^[0-9]+$", gene_ids[1:min(100, length(gene_ids))]))

print(paste("ENSEMBL format detected:", is_ensembl))
print(paste("Probe ID format detected:", is_probe))
print(paste("Gene symbol format detected:", is_symbol))

# ============================================================================
# 2. RUN DESEQ2 ANALYSIS
# ============================================================================

print("Running DESeq2... (this may take 1-2 minutes)")
dds$condition <- factor(dds$condition, levels = c("Normal", "Cancer"))
dds <- DESeq(dds)
print("DESeq2 analysis complete")

# ============================================================================
# 3. GET RESULTS
# ============================================================================

res <- results(dds, contrast = c("condition", "Cancer", "Normal"))

# Convert to data frame
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id")

print("Results obtained")
print(head(res_df, 10))

# ============================================================================
# 4. ADD GENE SYMBOLS (ADAPTIVE TO GENE ID TYPE)
# ============================================================================

print("=== MAPPING GENE IDS ===")

 if (is_probe) {
  print("Detected: Probe IDs (microarray) - using gene_id as identifier")
  
  # For probe IDs, we can't easily map
  # But let's try to extract info if possible
  res_df$symbol <- res_df$gene_id
  res_df$entrezid <- NA
  
  print("Note: Probe IDs don't map directly to gene symbols.")
  print("This is normal for microarray data. Analysis will proceed with probe IDs.")
  
} else {
  print("Detected: Other format (possibly Entrez IDs) - attempting mapping")
  
  res_df$symbol <- mapIds(
    org.Hs.eg.db,
    keys = res_df$gene_id,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  
  res_df$entrezid <- res_df$gene_id
}

print("Gene mapping complete")
print(paste("Genes with symbols:", sum(!is.na(res_df$symbol))))

# ============================================================================
# 5. FILTER SIGNIFICANT GENES
# ============================================================================

print("=== FILTERING SIGNIFICANT GENES ===")

sig_threshold_padj <- 0.05
sig_threshold_logfc <- 1.0

de_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < sig_threshold_padj) %>%
  filter(abs(log2FoldChange) > sig_threshold_logfc) %>%
  arrange(padj)

print(paste("Number of significantly DE genes:", nrow(de_genes)))
print(paste("  - Upregulated in Cancer:", sum(de_genes$log2FoldChange > 0)))
print(paste("  - Downregulated in Cancer:", sum(de_genes$log2FoldChange < 0)))

# ============================================================================
# 6. VOLCANO PLOT
# ============================================================================

print("Creating volcano plot...")

vol_data <- res_df %>%
  filter(!is.na(padj)) %>%
  mutate(
    sig = case_when(
      padj < sig_threshold_padj & abs(log2FoldChange) > sig_threshold_logfc ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )

p_volcano <- ggplot(vol_data, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_vline(xintercept = c(-sig_threshold_logfc, sig_threshold_logfc), 
             linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = -log10(sig_threshold_padj), 
             linetype = "dashed", color = "gray50") +
  scale_color_manual(
    values = c("Significant" = "red", "Not Significant" = "gray70"),
    name = "Significance"
  ) +
  labs(
    title = "Volcano Plot: Pancreatic Cancer vs Normal",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

ggsave("results/plots/03_volcano_plot.png", p_volcano, width = 10, height = 7, dpi = 300)
print("✓ Volcano plot saved")

# ============================================================================
# 7. SAVE RESULTS
# ============================================================================

print("Saving results...")

# Save all results
write.csv(res_df, "results/tables/03_all_genes_results.csv", row.names = FALSE)

# Save significant genes
write.csv(de_genes, "results/tables/03_significant_de_genes.csv", row.names = FALSE)

# Save for later analysis
saveRDS(res_df, "data/processed/deseq2_results.rds")
saveRDS(de_genes, "data/processed/de_genes_significant.rds")

print("✓ Results saved to results/tables/")

# ============================================================================
# 8. SUMMARY
# ============================================================================

print("=== DIFFERENTIAL EXPRESSION COMPLETE ===")
print(paste("Gene ID type used:", if(is_ensembl) "ENSEMBL" else if(is_symbol) "Symbol" else if(is_probe) "Probe ID" else "Other"))
print("Check results/tables/ for output files")