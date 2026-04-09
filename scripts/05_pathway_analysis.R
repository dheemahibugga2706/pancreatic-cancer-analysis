# ============================================================================
# SCRIPT 05: PATHWAY ENRICHMENT ANALYSIS
# Identify biological pathways dysregulated in pancreatic cancer
# ============================================================================

setwd("~/Documents/pancreatic-cancer-analysis")

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
library(ggplot2)

# Load DE genes
de_genes <- readRDS("data/processed/de_genes_significant.rds")

print(paste("Analyzing", nrow(de_genes), "significant DE genes"))

# ============================================================================
# 1. PREPARE GENE LIST FOR ENRICHMENT
# ============================================================================

# Get Entrez IDs for DE genes
gene_list <- de_genes %>%
  filter(!is.na(entrezid)) %>%
  pull(entrezid)

print(paste("Genes for enrichment analysis:", length(gene_list)))

# ============================================================================
# 2. KEGG PATHWAY ENRICHMENT
# ============================================================================

print("Running KEGG pathway enrichment...")
kegg_enrich <- enrichKEGG(
  gene = gene_list,
  organism = "hsa",  # Human
  pvalueCutoff = 0.05
)

kegg_results <- as.data.frame(kegg_enrich)

if (nrow(kegg_results) > 0) {
  print("Top 10 KEGG pathways:")
  print(kegg_results[1:min(10, nrow(kegg_results)), c("Description", "pvalue", "GeneRatio")])
  
  # Save results
  write.csv(kegg_results, "results/tables/05_kegg_pathways.csv", row.names = FALSE)
  
  # Plot top pathways
  p_kegg <- dotplot(kegg_enrich, showCategory = 10, title = "Top 10 KEGG Pathways")
  ggsave("results/plots/05_kegg_pathways_dotplot.png", p_kegg, width = 10, height = 7, dpi = 300)
  print("KEGG plot saved")
} else {
  print("No significant KEGG pathways found")
}

# ============================================================================
# 3. GO ENRICHMENT (BIOLOGICAL PROCESS)
# ============================================================================

print("Running GO enrichment (Biological Process)...")
go_enrich_bp <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",  # Biological Process
  pvalueCutoff = 0.05
)

go_results_bp <- as.data.frame(go_enrich_bp)

if (nrow(go_results_bp) > 0) {
  print("Top 10 GO Biological Process terms:")
  print(go_results_bp[1:min(10, nrow(go_results_bp)), c("Description", "pvalue", "GeneRatio")])
  
  write.csv(go_results_bp, "results/tables/05_go_bp_terms.csv", row.names = FALSE)
  
  p_go <- dotplot(go_enrich_bp, showCategory = 10, title = "Top 10 GO Biological Process Terms")
  ggsave("results/plots/05_go_bp_dotplot.png", p_go, width = 10, height = 7, dpi = 300)
  print("GO BP plot saved")
} else {
  print("No significant GO terms found")
}

# ============================================================================
# 4. CANCER-RELATED PATHWAY ANALYSIS
# ============================================================================

print("Checking cancer-related pathways...")

# Known pathways in pancreatic cancer
cancer_pathways <- c(
  "KRAS signaling",
  "TGF-beta signaling",
  "PI3K/AKT/mTOR signaling",
  "p53 pathway",
  "Wnt signaling"
)

# Check if these appear in results
print("Cancer pathways found in KEGG:")
for (pathway in cancer_pathways) {
  if (any(grepl(pathway, kegg_results$Description, ignore.case = TRUE))) {
    match_rows <- gregg_results %>% filter(grepl(pathway, Description, ignore.case = TRUE))
    for (i in 1:nrow(match_rows)) {
      print(paste(" ", match_rows$Description[i]))
    }
  }
}

print("=== PATHWAY ANALYSIS COMPLETE ===")