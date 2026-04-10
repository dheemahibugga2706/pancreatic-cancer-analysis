# Pancreatic Cancer Gene Expression Analysis

## Project Overview

This project analyzes gene expression data from pancreatic cancer patient samples to identify dysregulated genes and biological pathways associated with disease progression. The analysis integrates clinical context with computational biology, providing insights relevant to both surgical oncology and precision medicine.

### Clinical Relevance
Pancreatic cancer is one of the most aggressive human malignancies, often presenting as a surgical emergency. Understanding the molecular landscape of this disease is crucial for:
- Early detection strategies
- Prognosis prediction
- Therapeutic target identification
- Personalized treatment planning

## Dataset

- **Source:** Gene Expression Omnibus (GEO) - Dataset GSE71989
- **Samples:** Pancreatic cancer patient tissues and normal controls
- **Technology:** Gene expression microarray
- **Analysis:** Comparative transcriptomics (Cancer vs Normal)

## Analysis Pipeline

### Phase 1: Quality Control
**Script:** [`01_quality_control.R`](scripts/01_quality_control.R)

Performs initial assessment of data quality:
- Library size distribution across samples
- Expression value statistics  
- Gene expression distribution
- Identification of low-expression genes
- Missing value detection

**Output:** Quality control plots and filtered gene expression matrix
- Plot: [`results/plots/01_library_size.png`](results/plots/01_library_size.png)
- Plot: [`results/plots/01_expression_distribution.png`](results/plots/01_expression_distribution.png)


---

### Phase 2: Data Preparation
**Script:** [`02_differential_expression.R`](scripts/02_differential_expression.R)

Prepares data for differential expression analysis:
- Sample metadata organization
- Condition assignment (cancer vs. normal)
- Gene count matrix creation
- Low-count gene filtering
- DESeq2 object initialization

**Output:** Processed data objects ready for DESeq2

---

### Phase 3: Differential Expression Analysis
**Script:** [`03_differential_expression.R`](scripts/03_differential_expression.R)

Identifies genes significantly dysregulated in pancreatic cancer:
- DESeq2 statistical testing
- p-value adjustment (Benjamini-Hochberg FDR)
- Gene ID mapping to symbols and Entrez IDs
- Significance filtering (padj < 0.05, |log2FC| > 1)
- Volcano plot visualisation

**Output:** Complete results table and volcano plot
- Results table: [`results/tables/03_all_genes_results.csv`](results/tables/03_all_genes_results.csv)
- Significant genes: [`results/tables/03_significant_de_genes.csv`](results/tables/03_significant_de_genes.csv)
- Volcano Plot: [Volcano Plot](results/plots/03_volcano_plot.png)

---

### Phase 4: Cancer Gene Validation
**Script:** [`04_check_pancreatic_genes.R`](scripts/04_check_pcancer_genes.R)

Validates known pancreatic cancer genes in the analysis

**Output:** Console report of key gene changes

---

### Phase 5: Pathway Enrichment
**Script:** [`05_pathway_analysis.R`](scripts/05_pathway_analysis.R)

Identifies biological pathways dysregulated in pancreatic cancer:
- KEGG pathway enrichment
- Gene Ontology enrichment
- Cancer-related pathway identification

**Output:** No significant gene identified

---

### Phase 6: Visualisation
**Script:** [`06_visualisation.R`](scripts/06_visualization.R)

Creates publication-quality visualisations:
- PCA plot
- MA plot
- Volcano plot

**Output:** PNG files in `results/plots/`
- **PCA Plot:** [PCA Analysis](results/plots/06_pca_plot.png)
**PCA Plot:** Sample clustering reveals separation between pancreatic cancer and normal tissues based on gene expression profiles.

- **MA Plot:**  [MA Analysis](results/plots/06_ma_plot.png)

- **Volcano Plot:** [Volcano Plot](results/plots/06_volcano_plot.png)
**Volcano Plot:** Shows log2 fold change vs. -log10(adjusted p-value). Red points represent significant genes (padj < 0.05, |logFC| > 1).
---

### Phase 7: Summary Report
**Script:** [`07_summary_report.R`](scripts/07_summary_report.R)

Generates final analysis summary and interpretation

**Output:** Text summary in `results/ANALYSIS_SUMMARY_1.txt`

---

## Installation & Requirements

System Requirements

**R:** Version 4.0 or higher
**RStudio:** Recommended (optional but helpful)
**OS:** macOS, Linux, or Windows
**Disk Space:** ~2 GB

R Package Installation
-
```r
# Install CRAN packages
install.packages(c(
  "tidyverse",      # Data manipulation & plotting
  "ggplot2",        # Advanced graphics
  "pheatmap",       # Heatmap creation
  "RColorBrewer",   # Color palettes
  "gridExtra",      # Combine plots
  "renv"            # Environment management
))

# Install Bioconductor packages
BiocManager::install(c(
  "DESeq2",         # Differential expression
  "limma",          # Linear models
  "org.Hs.eg.db",   # Human genome annotations
  "clusterProfiler", # Pathway enrichment
  "GEOquery"        # Download from GEO
))
```
---

## Usage

## Quickstart: Run full pipeline
```r
# Set working directory
setwd("~/Documents/pancreatic-cancer-analysis")

# Run all analysis scripts sequentially
source("scripts/01_quality_control.R")
source("scripts/02_differential_expression.R")
source("scripts/03_differential_expression.R")
source("scripts/04_check_pancreatic_genes.R")
source("scripts/05_pathway_analysis.R")
source("scripts/06_visualisation.R")
source("scripts/07_summary_report.R")

print("✓ Analysis complete!")
```

### Run Full Transcripts
```r
# Just quality control
source("scripts/01_quality_control.R")

# Just differential expression
source("scripts/03_differential_expression.R")

# Just visualization
source("scripts/06_visualisation.R")
```

### Data Access
Raw data is downloaded directly from GEO within the scripts:

```r
# Automatic download (in script 02):
library(GEOquery)
gse <- getGEO("GSE71989", GSEMatrix = TRUE)
```
Processed data is stored as RDS files for easy loading:
```r
rdds <- readRDS("data/processed/deseq2_object.rds")
de_genes <- readRDS("data/processed/de_genes_significant.rds")
```
---

### Key Findings

Differential Expression Summary

**Total genes analyzed:** ~54,000 (after filtering)

**Total samples analysed:** 288

**Significantly dysregulated genes:** 2


*Upregulated in cancer:* [0]

*Downregulated in cancer:* [2]

Statistical threshold: FDR-adjusted p-value < 0.05, |log2 fold change| > 1

---
## Pathways Identified

No significant pathways were identified

Reason:
> Sample size was insginificant

>  No significant genes were detectable


### Known Pancreatic Cancer Genes

The following well-characterized pancreatic cancer genes showed changes consistent with literature:

*KRAS:* Not significant

*TP53:* Not significant

*CDKN2A:* not significant

*SMAD4:* Not significant

---

### Methods

**Data Source**

*Expression data obtained from:* 

Gene Expression Omnibus (GEO) accession GSE71989,

containing microarray data from pancreatic cancer and normal tissue samples.

### Quality Control

Removed genes with mean expression < 5
Assessed library size distribution
Checked for missing values and outliers

---

## Statistical Analysis

**Software**: DESeq2 (Bioconductor)

Normalization: Median of ratios normalization

**Test:** 

Wald test with Benjamini-Hochberg FDR correction

*Significance Criteria:*


Adjusted p-value (padj) < 0.05
Absolute log2 fold change > 1.0

---

## Clinical Implications

**Therapeutic Targets:** Upregulated genes represent potential intervention points for pancreatic cancer treatment.

**Biomarker Discovery:** Dysregulated genes serve as diagnostic (early detection), prognostic (patient stratification), and predictive (treatment response) biomarkers.

**Precision Medicine:** Molecular profiling enables patient stratification, targeted therapy selection, and treatment monitoring.

**Surgical Considerations:** Molecular insights inform surgical strategy, intervention timing, and neoadjuvant therapy planning.


## Limitations

- **Single dataset:** Results need validation in independent cohorts

- **Technology:** Microarray platform may differ from RNA-seq

- **Sample size:** Limited statistical power with small cohorts

- **Functional validation:** Computational predictions require experimental confirmation

- **Temporal scope:** Snapshot analysis, not dynamic gene expression changes

## Future Directions

- **Validation:** RNA-seq confirmation, protein-level validation (Western blot, immunohistochemistry), functional studies in cancer cell lines

- **Expanded analysis:** Multi-cohort meta-analysis, subtype-specific analysis, time-series analysis

- **Clinical translation:** Development of diagnostic/prognostic tests, prospective validation studies

---

# References


### Key Papers

1. Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.


2. He Y, Liu Y, Gong J, Liu C, Zhang H and Wu H: Identification of key pathways and candidate genes in pancreatic ductal adenocarcinoma using bioinformatics analysis. Oncol Lett 17: 3751-3764, 2019.


3. Ma Y, Pu Y, Peng L, Luo X, Xu J, Peng Y and Tang X: Identification of potential hub genes associated with the pathogenesis and prognosis of pancreatic duct adenocarcinoma using bioinformatics meta‑analysis of multi‑platform datasets. Oncol Lett 18: 6741-6751, 2019.


### Databases & Tools

Gene Expression Omnibus: https://www.ncbi.nlm.nih.gov/geo/

DESeq2: https://bioconductor.org/packages/DESeq2/

KEGG: https://www.genome.jp/kegg/

Reactome: https://reactome.org/


 ### Pancreatic Cancer Resources

National Cancer Institute: https://www.cancer.gov/types/pancreatic

Pancreatic Cancer Action: https://pancreaticcanceraction.org/

---

## Author & Contact

**Author:** Dheemahi Sai Sri Lakshmi Bugga

**Affiliation:** Malla Reddy Medical College For Women, Fourth year MBBS Student

**Research Interests:**
- Cancer Genomics
- Translational Oncology
- Precision Medicine
- Sugical Oncology

**Email:** [dheemahibugga27@gmail.com]

----

# License

This project is licensed under the MIT License. See the LICENSE file for details.

MIT License Summary

---

## Acknowledgments

GEO (Gene Expression Omnibus) for public access to data

Bioconductor team for DESeq2 and related packages

R community for excellent statistical packages
