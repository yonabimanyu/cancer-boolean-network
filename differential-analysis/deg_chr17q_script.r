# ==============================================================================
# DIFFERENTIAL EXPRESSION GENE (DEG) ANALYSIS - Chr17q Copy Number Status
# ==============================================================================
#
# Title: Differential Expression Analysis Based on Chr17q Copy Number Status
#
# Description: 
#   This script performs differential expression analysis comparing tumor samples
#   with Chr17q copy number alterations (GAIN vs Disomic) against normal controls.
#   The analysis uses limma package for statistical testing and identifies genes
#   significantly affected by Chr17q copy number changes in breast cancer.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Input Files:
#   - data/exp_combined_2.rds          : Combined RNA-seq expression data (SummarizedExperiment)
#   - data/copy_number_prepared.rds    : Copy number status annotations
#   - data/gene_info_final_3.5.rds     : Gene annotation reference (Ensembl to Symbol)
#
# Output Files:
#   - results/deg_chr17q_combined.csv          : All DEG results (3 contrasts)
#   - results/deg_chr17q_significant.csv       : Significant DEGs (|logFC|>1, FDR<0.05)
#   - results/deg_gain_vs_dis.csv              : Primary contrast (GAIN vs DIS)
#   - results/volcano_plot_gain_vs_control.pdf : Volcano plot visualization
#
# Dependencies:
#   - data.table (>= 1.14.0)
#   - dplyr (>= 1.0.0)
#   - limma (>= 3.48.0)
#   - EnhancedVolcano (>= 1.10.0)
#
# ==============================================================================

# ==============================================================================
# LIBRARY IMPORTS
# ==============================================================================

library(data.table)
library(dplyr)
library(limma)
library(EnhancedVolcano)

# ==============================================================================
# SECTION 1: DATA LOADING AND METADATA INTEGRATION
# ==============================================================================

# Load pre-processed expression data
# Note: This assumes exp_combined_2 object has been created in previous steps
# and contains combined tumor/normal RNA-seq data from TCGA BRCA

# Extract metadata from RNA-seq SummarizedExperiment object
meta_data_exp <- colData(exp_combined_2)
meta_data_exp_df <- as.data.frame(meta_data_exp)

# Create truncated sample IDs (15 characters) for matching with copy number data
# Rationale: TCGA barcodes follow hierarchical structure, first 15 chars identify patient-sample
meta_data_exp_df$Sample_15 <- substr(meta_data_exp_df$sample, 1, 15)

# Join with copy number annotations to assign Chr17q status
# Methodology: Integrate genomic alteration data (copy number) with transcriptomic data
meta_data_exp_df <- meta_data_exp_df %>%
  left_join(copy_number_prepared[, c("SAMPLE_ID_15", "Chr17q_Status")], 
            by = c("Sample_15" = "SAMPLE_ID_15"))

# Assign "Control" label to normal tissue samples
# Keep GAIN/DIS labels for tumor samples based on copy number status
# Scientific rationale: Normal tissue serves as biological baseline for comparison
meta_data_exp_df$Chr17q_Status_Final <- ifelse(
  meta_data_exp_df$sample_type == "Solid Tissue Normal", 
  "Control",
  meta_data_exp_df$Chr17q_Status
)

# Convert to factor with explicit ordering: Control < DIS < GAIN
# Rationale: Ordered factor ensures consistent contrast interpretation
meta_data_exp_df$Chr17q_Status_Final <- factor(
  meta_data_exp_df$Chr17q_Status_Final, 
  levels = c("Control", "DIS", "GAIN")
)

# Convert back to DFrame for compatibility with SummarizedExperiment
meta_data_chr17q_exp <- as(meta_data_exp_df, "DFrame")

# Display sample distribution across Chr17q status groups
cat("Sample distribution by Chr17q Status:\n")
print(table(meta_data_chr17q_exp$Chr17q_Status_Final))

# ==============================================================================
# SECTION 2: GENE EXPRESSION FILTERING BY VARIANCE
# ==============================================================================

# Calculate variance for each gene across all samples
# Rationale: Genes with zero variance provide no discriminatory information
# and can cause numerical instability in linear modeling
variances_chr17q <- apply(log2_combined_matrix_2, 1, var)

# Filter out genes with zero variance
filtered_matrix_chr17q <- log2_combined_matrix_2[variances_chr17q > 0, ]

cat("Genes before variance filtering:", nrow(log2_combined_matrix_2), "\n")
cat("Genes after variance filtering:", nrow(filtered_matrix_chr17q), "\n")

# Filter metadata to include only samples with complete Chr17q status
# Scientific rationale: Exclude samples with ambiguous copy number status
# to maintain statistical rigor
filtered_meta_data_chr17q <- meta_data_chr17q_exp[
  meta_data_chr17q_exp$Chr17q_Status_Final %in% c("Control", "DIS", "GAIN"), 
]
rownames(filtered_meta_data_chr17q) <- filtered_meta_data_chr17q$barcode

# Align expression matrix with filtered metadata
# Ensures sample order consistency between phenotype and expression data
filtered_matrix_chr17q <- filtered_matrix_chr17q[, rownames(filtered_meta_data_chr17q)]

cat("Expression matrix dimensions:", paste(dim(filtered_matrix_chr17q), collapse = " x "), "\n")
cat("Metadata dimensions:", paste(dim(filtered_meta_data_chr17q), collapse = " x "), "\n")
cat("Sample distribution by Chr17q Status (after filtering):\n")
print(table(filtered_meta_data_chr17q$Chr17q_Status_Final))

# ==============================================================================
# SECTION 3: GENE ANNOTATION - ENSEMBL ID TO GENE SYMBOL MAPPING
# ==============================================================================

# Remove version numbers from Ensembl IDs (e.g., ENSG00000012048.23 -> ENSG00000012048)
# Rationale: Version numbers can cause mismatches with annotation databases
ensembl_ids_clean <- gsub("\\..*", "", rownames(filtered_matrix_chr17q))

# Create data frame for merging with gene annotation reference
matrix_df <- data.frame(
  ensembl_id = ensembl_ids_clean,
  row_index = 1:nrow(filtered_matrix_chr17q),
  stringsAsFactors = FALSE
)

# Merge with gene information to obtain gene symbols
# Scientific rationale: Gene symbols are more interpretable than Ensembl IDs
# for biological interpretation and literature comparison
matrix_gene_info <- merge(
  matrix_df,
  gene_info_final_3.5[, c("ensembl_id", "gene_symbol")],
  by = "ensembl_id",
  all.x = TRUE
)

# Quality control: Identify genes without valid annotations
cat("Genes without Ensembl ID:", 
    sum(is.na(matrix_gene_info$ensembl_id) | matrix_gene_info$ensembl_id == ""), "\n")
cat("Genes without gene symbol:", 
    sum(is.na(matrix_gene_info$gene_symbol) | matrix_gene_info$gene_symbol == ""), "\n")

# Identify Ensembl IDs that failed to map to gene symbols
missing_symbol <- matrix_gene_info[
  !is.na(matrix_gene_info$ensembl_id) & 
  matrix_gene_info$ensembl_id != "" & 
  (is.na(matrix_gene_info$gene_symbol) | matrix_gene_info$gene_symbol == ""), 
]

cat("Ensembl IDs without gene_symbol mapping:", nrow(missing_symbol), "\n")

# Remove genes without valid gene symbols
# Rationale: Unannotated genes cannot be interpreted biologically
matrix_gene_info <- matrix_gene_info[
  !(is.na(matrix_gene_info$gene_symbol) | matrix_gene_info$gene_symbol == ""), 
]

cat("Genes after annotation filtering:", nrow(matrix_gene_info), "\n")

# ==============================================================================
# SECTION 4: REPRESENTATIVE GENE SELECTION FOR DUPLICATE SYMBOLS
# ==============================================================================

# Convert to data.table for efficient grouped operations
filtered_matrix_dt <- data.table(
  ensembl_id = matrix_gene_info$ensembl_id,
  gene_symbol = matrix_gene_info$gene_symbol,
  row_index = matrix_gene_info$row_index,
  filtered_matrix_chr17q[matrix_gene_info$row_index, ]
)

# Verify specific gene presence (quality control check)
cat("BRCA1 entries found:", sum(filtered_matrix_dt$gene_symbol == "BRCA1", na.rm = TRUE), "\n")

# Identify sample columns (exclude metadata columns)
sample_cols_exp <- names(filtered_matrix_dt)[
  !names(filtered_matrix_dt) %in% c("ensembl_id", "gene_symbol", "row_index")
]

# Select one representative Ensembl ID per gene symbol based on highest MAD
# Methodology: Multiple Ensembl IDs can map to the same gene symbol due to
# alternative transcripts or annotation ambiguity. We select the most variable
# transcript using Median Absolute Deviation (MAD) as a robust variance measure.
representative_genes_2 <- filtered_matrix_dt[, {
  if(length(sample_cols_exp) > 0) {
    gene_vars <- apply(.SD[, sample_cols_exp, with = FALSE], 1, function(x) {
      mad(as.numeric(x), na.rm = TRUE, constant = 1.4826)
    })
    .SD[which.max(gene_vars)]
  }
}, by = gene_symbol]

cat("Genes before duplicate removal:", nrow(filtered_matrix_dt), "\n")
cat("Genes after selecting representatives:", nrow(representative_genes_2), "\n")

# Verify uniqueness: Each gene symbol should have exactly one Ensembl ID
duplicate_check <- representative_genes_2 %>%
  distinct(gene_symbol, ensembl_id) %>%
  dplyr::count(gene_symbol) %>%
  filter(n > 1)

if(nrow(duplicate_check) > 0) {
  cat("WARNING: Duplicate gene symbols still present!\n")
} else {
  cat("SUCCESS: All gene symbols are unique\n")
}

# ==============================================================================
# SECTION 5: FINAL EXPRESSION MATRIX CONSTRUCTION
# ==============================================================================

# Construct final expression matrix with representative genes
filtered_matrix_chr17q_final <- as.matrix(
  representative_genes_2[, sample_cols_exp, with = FALSE]
)
rownames(filtered_matrix_chr17q_final) <- representative_genes_2$ensembl_id

cat("Final expression matrix dimensions:", 
    paste(dim(filtered_matrix_chr17q_final), collapse = " x "), "\n")

# ==============================================================================
# SECTION 6: STATISTICAL DESIGN - LINEAR MODEL SETUP
# ==============================================================================

# Create design matrix without intercept (~ 0 + Chr17q_Status)
# Methodology: No-intercept model allows direct comparison between all groups
# Each coefficient represents the mean expression for that group
design_chr17q <- model.matrix(~ 0 + Chr17q_Status_Final, data = filtered_meta_data_chr17q)
colnames(design_chr17q) <- levels(filtered_meta_data_chr17q$Chr17q_Status_Final)

cat("Design matrix columns:\n")
print(colnames(design_chr17q))

# ==============================================================================
# SECTION 7: LINEAR MODEL FITTING WITH LIMMA
# ==============================================================================

# Fit linear model for each gene
# Methodology: Limma uses empirical Bayes shrinkage to improve variance estimates
# particularly beneficial for RNA-seq data with limited replicates
fit_chr17q <- lmFit(filtered_matrix_chr17q_final, design_chr17q)

cat("Genes in fitted model:", nrow(fit_chr17q$coefficients), "\n")
cat("Model coefficients:", colnames(fit_chr17q$coefficients), "\n")

# ==============================================================================
# SECTION 8: CONTRAST DEFINITIONS FOR PAIRWISE COMPARISONS
# ==============================================================================

# Define contrast matrix for three key comparisons:
# 1. GAIN vs DIS: Primary contrast - effect of Chr17q gain on disomic background
# 2. GAIN vs Control: Tumor with gain vs normal tissue
# 3. DIS vs Control: Disomic tumor vs normal tissue
# Scientific rationale: These contrasts isolate the specific effect of Chr17q
# gain while accounting for general tumor vs normal differences
contrast_matrix_chr17q <- makeContrasts(
  "GAIN_vs_DIS" = GAIN - DIS,    
  "GAIN_vs_Control" = GAIN - Control,
  "DIS_vs_Control" = DIS - Control,
  levels = design_chr17q
)

cat("Contrast matrix:\n")
print(contrast_matrix_chr17q)

# ==============================================================================
# SECTION 9: CONTRAST APPLICATION AND EMPIRICAL BAYES MODERATION
# ==============================================================================

# Apply contrast matrix to linear model fits
fit_chr17q <- contrasts.fit(fit_chr17q, contrast_matrix_chr17q)

# Apply empirical Bayes shrinkage to variance estimates
# Methodology: eBayes borrows information across genes to stabilize variance
# estimates, improving statistical power for differential expression detection
fit_chr17q <- eBayes(fit_chr17q)

# ==============================================================================
# SECTION 10: DIFFERENTIAL EXPRESSION RESULTS EXTRACTION
# ==============================================================================

# Initialize list to store results from each contrast
deg_results_chr17q_list <- list()

# Loop through each contrast and extract DEG results
for (contrast_name in colnames(contrast_matrix_chr17q)) {
  # Extract results with FDR correction for multiple testing
  # Methodology: Benjamini-Hochberg FDR controls expected proportion of false positives
  deg_result <- topTable(fit_chr17q, coef = contrast_name, adjust = "fdr", number = Inf)
  
  # Add contrast identifier and Ensembl ID
  deg_result$contrast <- contrast_name
  deg_result$ensembl_gene_id <- rownames(deg_result)
  
  # Store in list
  deg_results_chr17q_list[[contrast_name]] <- deg_result
}

# Combine all contrast results into single data frame
deg_results_chr17q_combined <- do.call(rbind, deg_results_chr17q_list)

cat("Total DEG results across all contrasts:", nrow(deg_results_chr17q_combined), "\n")

# Verify presence of key gene (quality control)
cat("BRCA1 present in results:", "BRCA1" %in% deg_results_chr17q_combined$gene_symbol, "\n")

# ==============================================================================
# SECTION 11: GENE ANNOTATION ENRICHMENT
# ==============================================================================

# Merge DEG results with complete gene annotation
# Adds gene symbol, Entrez ID, and chromosome information
deg_results_chr17q_combined <- merge(
  deg_results_chr17q_combined,
  gene_info_final_3.5[, c("ensembl_id", "gene_symbol", "entrez_id", "chromosome")],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_id",
  all.x = TRUE
)

# Quality control: Check annotation completeness
cat("Genes without Ensembl ID:", 
    sum(is.na(deg_results_chr17q_combined$ensembl_gene_id) | 
        deg_results_chr17q_combined$ensembl_gene_id == ""), "\n")
cat("Genes without gene symbol:", 
    sum(is.na(deg_results_chr17q_combined$gene_symbol) | 
        deg_results_chr17q_combined$gene_symbol == ""), "\n")

# ==============================================================================
# SECTION 12: UNIQUE IDENTIFIER CREATION
# ==============================================================================

# Create unique row identifier by combining gene symbol and contrast
# Rationale: Allows unambiguous identification of each gene-contrast pair
deg_results_chr17q_combined$unique_id <- paste(
  deg_results_chr17q_combined$gene_symbol,
  deg_results_chr17q_combined$contrast,
  sep = "_"
)

# Verify uniqueness
unique_check <- length(unique(deg_results_chr17q_combined$unique_id)) == 
                nrow(deg_results_chr17q_combined)
cat("Unique ID verification:", ifelse(unique_check, "PASSED", "FAILED"), "\n")

# Set as row names for easy access
rownames(deg_results_chr17q_combined) <- deg_results_chr17q_combined$unique_id

# ==============================================================================
# SECTION 13: SIGNIFICANT GENE FILTERING
# ==============================================================================

# Define significance criteria: |log2FC| > 1 and FDR < 0.05
# Methodology: Log2FC threshold of 1 corresponds to 2-fold change, ensuring
# biological relevance. FDR < 0.05 controls false discovery rate at 5%.
significant_chr17q_combined <- deg_results_chr17q_combined[
  abs(deg_results_chr17q_combined$logFC) > 1 & 
  deg_results_chr17q_combined$adj.P.Val < 0.05, 
]

# Create regulation status labels for visualization
significant_chr17q_combined$status <- ifelse(
  significant_chr17q_combined$logFC > 1, "OVER_T", 
  ifelse(significant_chr17q_combined$logFC < -1, "DOWN_T", "UNCHANGED")
)

cat("Total significant genes:", nrow(significant_chr17q_combined), "\n")

# Verify key gene presence
cat("BRCA1 in significant genes:", 
    "BRCA1" %in% significant_chr17q_combined$gene_symbol, "\n")

# Extract Chr17-specific significant genes across all contrasts
chr17_deg_all_contrasts <- significant_chr17q_combined[
  significant_chr17q_combined$chromosome == "17", 
]
cat("Chr17 significant genes (all contrasts):", nrow(chr17_deg_all_contrasts), "\n")

# ==============================================================================
# SECTION 14: PRIMARY CONTRAST - GAIN VS DIS
# ==============================================================================

# Extract primary contrast results: GAIN vs DIS
# Scientific rationale: This contrast isolates the effect of Chr17q gain
# while controlling for tumor background (both are tumor samples)
gain_vs_dis_deg <- significant_chr17q_combined[
  significant_chr17q_combined$contrast == "GAIN_vs_DIS", 
]

cat("Significant genes in GAIN vs DIS:", nrow(gain_vs_dis_deg), "\n")

# Focus on Chr17 genes (target chromosome)
chr17_gain_vs_dis_deg <- gain_vs_dis_deg[gain_vs_dis_deg$chromosome == "17", ]

cat("Chr17 significant genes in GAIN vs DIS:", nrow(chr17_gain_vs_dis_deg), "\n")
cat("BRCA1 in Chr17 GAIN vs DIS:", 
    "BRCA1" %in% chr17_gain_vs_dis_deg$gene_symbol, "\n")

# ==============================================================================
# SECTION 15: SECONDARY CONTRAST - GAIN VS CONTROL
# ==============================================================================

# Extract GAIN vs Control comparison results
gain_vs_ctrl_deg <- significant_chr17q_combined[
  significant_chr17q_combined$contrast == "GAIN_vs_Control", 
]

cat("Significant genes in GAIN vs Control:", nrow(gain_vs_ctrl_deg), "\n")

# Chr17-specific genes
chr17_gain_vs_ctrl_deg <- gain_vs_ctrl_deg[gain_vs_ctrl_deg$chromosome == "17", ]

cat("Chr17 significant genes in GAIN vs Control:", nrow(chr17_gain_vs_ctrl_deg), "\n")
cat("BRCA1 in Chr17 GAIN vs Control:", 
    "BRCA1" %in% chr17_gain_vs_ctrl_deg$gene_symbol, "\n")

# ==============================================================================
# SECTION 16: TERTIARY CONTRAST - DIS VS CONTROL
# ==============================================================================

# Extract DIS vs Control comparison results
dis_vs_ctrl_deg <- significant_chr17q_combined[
  significant_chr17q_combined$contrast == "DIS_vs_Control", 
]

cat("Significant genes in DIS vs Control:", nrow(dis_vs_ctrl_deg), "\n")

# Chr17-specific genes
chr17_dis_vs_ctrl_deg <- dis_vs_ctrl_deg[dis_vs_ctrl_deg$chromosome == "17", ]

cat("Chr17 significant genes in DIS vs Control:", nrow(chr17_dis_vs_ctrl_deg), "\n")
cat("BRCA1 in Chr17 DIS vs Control:", 
    "BRCA1" %in% chr17_dis_vs_ctrl_deg$gene_symbol, "\n")

# ==============================================================================
# SECTION 17: VISUALIZATION - VOLCANO PLOT
# ==============================================================================

# Prepare data for volcano plot (GAIN vs Control contrast)
deg_gain_vs_ctrl <- deg_results_chr17q_combined[
  deg_results_chr17q_combined$contrast == "GAIN_vs_Control", 
]

# Filter for Chr17 genes only
deg_chr17q_gain_vs_ctrl <- deg_results_chr17q_combined[
  deg_results_chr17q_combined$contrast == "GAIN_vs_Control" &
  deg_results_chr17q_combined$chromosome == "17",
]

# Create enhanced volcano plot
# Visualization rationale: Volcano plots simultaneously display fold change
# magnitude and statistical significance, ideal for identifying key genes
volcano_plot <- EnhancedVolcano(
  toptable = deg_chr17q_gain_vs_ctrl,            
  lab = deg_chr17q_gain_vs_ctrl$gene_symbol,    
  x = 'logFC',                    
  y = 'adj.P.Val',              
  title = 'DEG Analysis: GAIN-Chr17q vs Control',
  xlab = bquote(~Log[2] ~ "Fold Change"),
  ylab = bquote(~-Log[10] ~ "Adjusted P-value"),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,          
  labSize = 5.0,
  colAlpha = 0.7,
  legendPosition = 'none'
)

# Display plot
print(volcano_plot)

# Save plot to file
# pdf("results/volcano_plot_gain_vs_control.pdf", width = 10, height = 8)
# print(volcano_plot)
# dev.off()

# ==============================================================================
# SECTION 18: RESULTS SUMMARY
# ==============================================================================

cat("\n=== DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY ===\n")
cat("Total genes analyzed:", nrow(deg_results_chr17q_combined)/3, "\n")
cat("Significant genes (all contrasts):", nrow(significant_chr17q_combined), "\n")
cat("Significant genes - GAIN vs DIS:", 
    sum(significant_chr17q_combined$contrast == "GAIN_vs_DIS"), "\n")
cat("Significant genes - GAIN vs Control:", 
    sum(significant_chr17q_combined$contrast == "GAIN_vs_Control"), "\n")
cat("Significant genes - DIS vs Control:", 
    sum(significant_chr17q_combined$contrast == "DIS_vs_Control"), "\n")
cat("Chr17 significant genes (GAIN vs DIS):", nrow(chr17_gain_vs_dis_deg), "\n")

# ==============================================================================
# SECTION 19: DATA EXPORT
# ==============================================================================

# Export results to CSV files for downstream analysis
# write.csv(deg_results_chr17q_combined, 
#           "results/deg_chr17q_all_results.csv", 
#           row.names = FALSE)
# 
# write.csv(significant_chr17q_combined, 
#           "results/deg_chr17q_significant.csv", 
#           row.names = FALSE)
# 
# write.csv(gain_vs_dis_deg, 
#           "results/deg_gain_vs_dis.csv", 
#           row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
