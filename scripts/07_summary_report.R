# ============================================================================
# SCRIPT 07: GENERATE SUMMARY REPORT
# ============================================================================

setwd("~/Documents/pancreatic-cancer-analysis")

library(tidyverse)

# Load results
de_genes <- readRDS("data/processed/de_genes_significant.rds")

# Create summary text
summary_report <- paste(
  "PANCREATIC CANCER GENE EXPRESSION ANALYSIS - SUMMARY REPORT",
  "=" %>% rep(50) %>% paste(collapse = ""),
  "",
  "DATASET INFORMATION:",
  "- Dataset: GSE71989 (Pancreatic Cancer Patient Samples)",
  "- Number of samples: [to be filled]",
  "- Technology: Gene expression microarray",
  "",
  "DIFFERENTIAL EXPRESSION ANALYSIS:",
  paste("- Total genes analyzed:", nrow(de_genes)),
  paste("- Upregulated in Cancer:", sum(de_genes$log2FoldChange > 0)),
  paste("- Downregulated in Cancer:", sum(de_genes$log2FoldChange < 0)),
  "- Statistical threshold: padj < 0.05, |log2FC| > 1",
  "",
  "TOP 10 UPREGULATED GENES IN PANCREATIC CANCER:",
  "",
  sep = "\n"
)

# Add top upregulated genes
top_up <- de_genes %>%
  filter(log2FoldChange > 0) %>%
  head(10) %>%
  mutate(entry =paste0(row_number(), "."),
    symbol,
    paste0("(log2FC = ", round(log2FoldChange, 2), ", padj = ", 
           sprintf("%.2e", padj), ")")
  )

summary_report <- paste(
  summary_report,
  paste(top_up$entry, collapse = "\n"),
  "",
  "TOP 10 DOWNREGULATED GENES IN PANCREATIC CANCER:",
  "",
  sep = "\n"
)

# Add top downregulated genes
top_down <- de_genes %>%
  filter(log2FoldChange < 0) %>%
  head(10) %>%
  mutate(entry = paste0(row_number(), "."),
    symbol,
    paste0("(log2FC = ", round(log2FoldChange, 2), ", padj = ", 
           sprintf("%.2e", padj), ")")
  )

summary_report <- paste(
  summary_report,
  paste(top_down$entry, collapse = "\n"),
  "",
  "BIOLOGICAL INTERPRETATION:",
  "- This analysis identified genes significantly dysregulated in pancreatic cancer",
  "- Results span critical pathways including growth signaling, apoptosis, and metabolism",
  "- Findings have implications for understanding cancer progression and potential therapy targets",
  "",
  "CLINICAL RELEVANCE:",
  "- Pancreatic cancer is among the most aggressive malignancies",
  "- Early detection is crucial, as most cases present at advanced stages",
  "- Gene expression signatures may inform prognosis and treatment selection",
  "",
  sep = "\n"
)

# Save report
writeLines(summary_report, "results/ANALYSIS_SUMMARY_1.txt")

print(summary_report)