# ============================================================================
# SCRIPT 04: CHECK KEY PANCREATIC CANCER GENES
# ============================================================================

setwd("~/Documents/pancreatic-cancer-analysis")

library(tidyverse)

# Load results
de_genes <- readRDS("data/processed/de_genes_significant.rds")

# Key genes in pancreatic cancer
pancreatic_cancer_genes <- c(
  "KRAS",      # Most commonly mutated in pancreatic cancer
  "TP53",      # Tumor suppressor
  "CDKN2A",    # Tumor suppressor (p16)
  "PDAC",      # Pancreatic ductal adenocarcinoma marker
  "SMAD4",     # TGF-beta pathway
  "EGFR",      # Growth factor receptor
  "HER2",      # Growth factor receptor
  "CA19-9",    # Tumor marker
  "MUC1",      # Tumor antigen
  "PTEN",      # PI3K pathway inhibitor
  "AKT1"       # Survival pathway
)

print("=== KEY PANCREATIC CANCER GENES ===")

for (gene in pancreatic_cancer_genes) {
  if (gene %in% de_genes$symbol) {
    gene_data <- de_genes %>% filter(symbol == gene)
    print(paste(
      gene, ":",
      "logFC =", round(gene_data$log2FoldChange, 2),
      "| p-adj =", sprintf("%.2e", gene_data$padj)
    ))
  } else {
    print(paste(gene, ": Not significantly DE"))
  }
}

print("Analysis complete")