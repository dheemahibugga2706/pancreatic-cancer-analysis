# ============================================================================
# SCRIPT 02: DIFFERENTIAL EXPRESSION PREPARATION
# Prepare data for DESeq2 analysis
# ============================================================================

setwd("~/Documents/pancreatic-cancer-analysis")

library(tidyverse)
library(DESeq2)

# Load processed data
expr_matrix <- readRDS("data/processed/expression_filtered.rds")
sample_info <- readRDS("data/processed/sample_info_processed.rds")

print("Data loaded")

# ============================================================================
# 1. PREPARE SAMPLE INFO
# ============================================================================

# Check sample information columns
print("Sample info columns:")
print(colnames(sample_info))

# Create condition column (Pancreatic Cancer vs Normal)
# We'll identify this from the sample titles
sample_info$condition <- NA

# Assign conditions based on sample characteristics
# Look for keywords like "normal", "cancer", "tumor"
for (i in 1:nrow(sample_info)) {
  title <- tolower(sample_info$title[i])
  
  if (grepl("cancer|tumor|diseased", title)) {
    sample_info$condition[i] <- "Cancer"
  } else if (grepl("normal|control|healthy", title)) {
    sample_info$condition[i] <- "Normal"
  } else {
    # Default assignment based on patterns
    sample_info$condition[i] <- "Cancer"  # Most samples are cancer in this dataset
  }
}

print("Condition assignments:")
print(table(sample_info$condition))

# ============================================================================
# 2. PREPARE COUNTS MATRIX FOR DESEQ2
# ============================================================================

# DESeq2 expects integers, so round expression values
counts_matrix <- round(expr_matrix)

# Verify it's numeric
class(counts_matrix)
head(counts_matrix[, 1:5])  # Check first few entries

# ============================================================================
# 3. CREATE DESEQ2 OBJECT
# ============================================================================

# Match sample order between counts and sample info
counts_matrix <- counts_matrix[, rownames(sample_info)]

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData = sample_info,
  design = ~ condition
)

print(dds)

# ============================================================================
# 4. FILTER GENES WITH LOW COUNTS
# ============================================================================

# Keep genes with at least 10 reads in at least 2 samples
keep <- rowSums(counts(dds) >= 10) >= 2
dds <- dds[keep, ]

print(paste("Genes after filtering:", nrow(dds)))

# ============================================================================
# 5. SAVE FOR DIFFERENTIAL EXPRESSION
# ============================================================================

saveRDS(dds, "data/processed/deseq2_object.rds")

print("DESeq2 object prepared and saved")