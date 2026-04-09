# ============================================================================
# SCRIPT 01: QUALITY CONTROL (CORRECTED VERSION)
# Pancreatic Cancer Gene Expression Analysis
# ============================================================================
# This script performs initial quality control on the raw data

# Set working directory
setwd("~/Documents/pancreatic-cancer-analysis")

# Load libraries
library(tidyverse)
library(ggplot2)

# ============================================================================
# 1. LOAD DATA
# ============================================================================

print("Loading data...")
expr_matrix <- readRDS("data/raw/pancreatic_expression_matrix.rds")
sample_info <- readRDS("data/raw/pancreatic_sample_info.rds")

print("Data loaded successfully")
print(paste("Expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples"))

# ============================================================================
# 2. BASIC STATISTICS
# ============================================================================

print("=== EXPRESSION VALUE STATISTICS ===")

# Calculate statistics safely
expr_vector <- as.vector(as.matrix(expr_matrix))

summary_stats <- data.frame(
  Min = min(expr_vector, na.rm = TRUE),
  Q1 = quantile(expr_vector, 0.25, na.rm = TRUE),
  Median = median(expr_vector, na.rm = TRUE),
  Mean = mean(expr_vector, na.rm = TRUE),
  Q3 = quantile(expr_vector, 0.75, na.rm = TRUE),
  Max = max(expr_vector, na.rm = TRUE)
)

print(summary_stats)

# ============================================================================
# 3. CHECK FOR LOW EXPRESSION GENES
# ============================================================================

print("=== CALCULATING GENE MEANS ===")

# Calculate mean expression for each gene
gene_means <- rowMeans(expr_matrix, na.rm = TRUE)

print(paste("Genes with mean expression < 5:", sum(gene_means < 5, na.rm = TRUE)))
print(paste("Genes with mean expression >= 5:", sum(gene_means >= 5, na.rm = TRUE)))

# Filter out low expression genes
expr_filtered <- expr_matrix[gene_means >= 5, ]
print(paste("After filtering:", nrow(expr_filtered), "genes"))

# ============================================================================
# 4. SAMPLE QUALITY CHECK
# ============================================================================

print("=== CHECKING FOR MISSING VALUES ===")
print(paste("Missing values in expression matrix:", sum(is.na(expr_matrix))))
print(paste("Missing values in sample info:", sum(is.na(sample_info))))

# ============================================================================
# 5. LIBRARY SIZE (Total expression per sample)
# ============================================================================

print("=== CALCULATING LIBRARY SIZE ===")

library_size <- colSums(expr_matrix, na.rm = TRUE)

# Create data frame correctly
lib_size_df <- data.frame(
  sample = names(library_size),
  library_size = as.numeric(library_size),
  stringsAsFactors = FALSE  # IMPORTANT: prevents factors
)

print("Library size summary:")
print(summary(lib_size_df$library_size))

# ============================================================================
# 6. PLOT LIBRARY SIZES (CORRECTED)
# ============================================================================

print("Creating library size plot...")

# Sort data for plotting (without using xtfrm)
lib_size_df <- lib_size_df %>%
  mutate(sample = factor(sample, levels = sample[order(-library_size)]))

# Create plot
p_libsize <- ggplot(lib_size_df, aes(x = sample, y = library_size)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Library Size Across Samples",
    x = "Sample",
    y = "Total Expression"
  )

ggsave("results/plots/01_library_size.png", p_libsize, width = 10, height = 6, dpi = 300)
print("✓ Library size plot saved")

# ============================================================================
# 7. EXPRESSION DISTRIBUTION
# ============================================================================

print("Creating expression distribution plot...")

# Convert to long format carefully
expr_long <- expr_matrix %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "sample",
    values_to = "expression",
    values_drop_na = TRUE  # Remove NA values
  )

# Box plot
p_boxplot <- ggplot(expr_long, aes(x = sample, y = expression)) +
  geom_boxplot(fill = "lightblue", outlier.size = 0.5) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "Expression Distribution Across Samples",
    x = "Sample",
    y = "Expression Level"
  )

ggsave("results/plots/01_expression_distribution.png", p_boxplot, width = 10, height = 6, dpi = 300)
print("✓ Expression distribution plot saved")

# ============================================================================
# 8. SUMMARY STATISTICS
# ============================================================================

print("=== QUALITY CONTROL SUMMARY ===")
print(paste("Total genes analyzed:", nrow(expr_matrix)))
print(paste("Total samples:", ncol(expr_matrix)))
print(paste("Average library size:", round(mean(library_size), 0)))
print(paste("Min library size:", round(min(library_size), 0)))
print(paste("Max library size:", round(max(library_size), 0)))

# ============================================================================
# 9. SAVE FILTERED DATA FOR NEXT STEP
# ============================================================================

print("Saving filtered data...")

saveRDS(expr_filtered, "data/processed/expression_filtered.rds")
saveRDS(sample_info, "data/processed/sample_info_processed.rds")

print("✓ Data saved to data/processed/")

# ============================================================================
# 10. FINAL MESSAGE
# ============================================================================

print("=== QUALITY CONTROL COMPLETE ===")
print("Next step: Run 02_differential_expression_prep.R")