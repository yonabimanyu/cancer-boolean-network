# ==============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS - BRCA1 Expression Stratification
# ==============================================================================
#
# Title: Gene Expression Analysis Stratified by BRCA1 Expression Levels
#
# Description: 
#   This script stratifies Chr17q GAIN tumor samples based on BRCA1 expression
#   levels and performs differential expression analysis. Samples are classified
#   into upBRCA1 (high expression) and downBRCA1 (low expression) groups using
#   statistical thresholds (mean ± 1 SD) derived from the GAIN cohort distribution.
#   Identifies genes whose expression is associated with BRCA1 levels.
#
# Author: [Your Name]
# Date: 2025-01-XX
# Version: 1.0
#
# Input Files:
#   - data/log2_combined_matrix_2.rds      : Log2-transformed expression matrix
#   - data/meta_data_chr17q_exp.rds        : Metadata with Chr17q status
#   - data/gene_info_final_3.5.rds         : Gene annotation reference
#
# Output Files:
#   - results/deg_brca1_combined.csv              : All DEG results
#   - results/deg_brca1_significant.csv           : Significant DEGs
#   - results/deg_upbrca1_vs_downbrca1.csv        : Primary contrast
#   - results/volcano_plot_upbrca1_vs_downbrca1.pdf : Volcano plot
#   - results/brca1_stratification_summary.txt    : Classification summary
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
# SECTION 1: EXPRESSION DATA PREPROCESSING
# ==============================================================================

# Calculate variance for each gene and remove zero-variance genes
# Rationale: Zero-variance genes provide no information for differential
# expression and can cause numerical instability
cat("Calculating gene variance and filtering...\n")
variances_chr17q <- apply(log2_combined_matrix_2, 1, var)
filtered_matrix_chr17q <- log2_combined_matrix_2[variances_chr17q > 0, ]

cat("Genes before filtering:", nrow(log2_combined_matrix_2), "\n")
cat("Genes after filtering:", nrow(filtered_matrix_chr17q), "\n")

# ==============================================================================
# SECTION 2: METADATA PREPROCESSING
# ==============================================================================

# Filter metadata for samples with complete Chr17q status
# Rationale: Analysis focuses on samples with clear genomic classification
cat("Filtering metadata by Chr17q status...\n")
filtered_meta_data_chr17q <- meta_data_chr17q_exp[
  meta_data_chr17q_exp$Chr17q_Status_Final %in% c("Control", "DIS", "GAIN"), 
]
rownames(filtered_meta_data_chr17q) <- filtered_meta_data_chr17q$barcode

# Align expression matrix with filtered metadata
filtered_matrix_chr17q <- filtered_matrix_chr17q[, rownames(filtered_meta_data_chr17q)]

cat("Samples after filtering:", ncol(filtered_matrix_chr17q), "\n")

# ==============================================================================
# SECTION 3: SAMPLE IDENTIFICATION BY CHR17Q STATUS
# ==============================================================================

# Extract GAIN samples
# Rationale: BRCA1 stratification is performed within the GAIN cohort
# to isolate effects specifically related to BRCA1 expression levels
# independent of copy number status
gain_samples <- rownames(filtered_meta_data_chr17q)[
  filtered_meta_data_chr17q$Chr17q_Status_Final == "GAIN"
]
cat("GAIN samples:", length(gain_samples), "\n")

# Extract normal control samples
normal_samples <- rownames(filtered_meta_data_chr17q)[
  filtered_meta_data_chr17q$sample_type == "Solid Tissue Normal"
]
cat("Control samples:", length(normal_samples), "\n")

# ==============================================================================
# SECTION 4: COMBINE GAIN AND CONTROL SAMPLES
# ==============================================================================

# Merge metadata from GAIN and control groups
# Scientific rationale: Including controls allows assessment of whether
# BRCA1-associated expression patterns are tumor-specific or reflect
# normal tissue biology
meta_gain <- filtered_meta_data_chr17q[gain_samples, ]
meta_normal <- filtered_meta_data_chr17q[normal_samples, ]
meta_combined <- rbind(meta_gain, meta_normal)

# Extract expression for combined samples
combined_samples <- rownames(meta_combined)
exp_combined <- filtered_matrix_chr17q[, combined_samples]

cat("Total samples in analysis:", length(combined_samples), "\n")

# ==============================================================================
# SECTION 5: BRCA1 EXPRESSION EXTRACTION
# ==============================================================================

# Convert to data frame and add Ensembl IDs
exp_combined_df <- as.data.frame(exp_combined)
exp_combined_df$ENSEMBL_ID <- rownames(exp_combined_df)

# Extract BRCA1 expression by Ensembl ID
brca1_ensembl_id <- "ENSG00000012048.23"
cat("Extracting BRCA1 expression (", brca1_ensembl_id, ")...\n")

brca1_expression <- exp_combined_df[exp_combined_df$ENSEMBL_ID == brca1_ensembl_id, ]
brca1_expression$ENSEMBL_ID <- NULL
brca1_expression_matrix <- as.numeric(as.matrix(brca1_expression))

# Validation
if (nrow(exp_combined_df[exp_combined_df$ENSEMBL_ID == brca1_ensembl_id, ]) == 0) {
  stop("ERROR: BRCA1 Ensembl ID not found in data!")
}

# ==============================================================================
# SECTION 6: NORMALITY TEST OF BRCA1 EXPRESSION
# ==============================================================================

# Test normality only on GAIN samples
# Methodology rationale: Statistical thresholds (mean ± SD) are most
# appropriate for normally distributed data. Normality testing validates
# the use of parametric statistical methods for threshold determination.
cat("Performing normality test on GAIN samples...\n")

brca1_expression_gain_only <- brca1_expression[, colnames(exp_combined) %in% gain_samples]
brca1_expression_gain_vector <- as.numeric(brca1_expression_gain_only)

# Kolmogorov-Smirnov test for normality
ks_result <- ks.test(
  brca1_expression_gain_vector,
  "pnorm",
  mean = mean(brca1_expression_gain_vector, na.rm = TRUE),
  sd = sd(brca1_expression_gain_vector, na.rm = TRUE)
)

cat("KS test p-value for BRCA1 expression (GAIN samples):", ks_result$p.value, "\n")
if (ks_result$p.value > 0.05) {
  cat("Data follows normal distribution (p > 0.05)\n")
} else {
  cat("Data does not follow normal distribution (p <= 0.05)\n")
}

# ==============================================================================
# SECTION 7: BRCA1 EXPRESSION THRESHOLD DETERMINATION
# ==============================================================================

# Calculate thresholds based on mean and standard deviation
# Methodology: Mean ± 1 SD encompasses approximately 68% of the distribution
# (assuming normality). Samples beyond ±1 SD represent the tails of the
# distribution and exhibit distinctly high or low BRCA1 expression.
# Scientific rationale: This approach balances statistical rigor with
# adequate sample size in each group for downstream analysis.
cat("Calculating BRCA1 expression thresholds...\n")

mean_brca1_exp <- mean(brca1_expression_gain_vector, na.rm = TRUE)
sd_brca1_exp <- sd(brca1_expression_gain_vector, na.rm = TRUE)

# Thresholds: mean ± 1 SD
up_threshold <- mean_brca1_exp + 1.0 * sd_brca1_exp
down_threshold <- mean_brca1_exp - 1.0 * sd_brca1_exp

cat("Mean BRCA1 expression (GAIN samples):", round(mean_brca1_exp, 4), "\n")
cat("Standard deviation:", round(sd_brca1_exp, 4), "\n")
cat("Upper threshold (upBRCA1):", round(up_threshold, 4), "\n")
cat("Lower threshold (downBRCA1):", round(down_threshold, 4), "\n")

# ==============================================================================
# SECTION 8: SAMPLE CLASSIFICATION BY BRCA1 EXPRESSION
# ==============================================================================

# Extract average BRCA1 expression per sample
brca1_expression_avg <- as.numeric(brca1_expression_matrix)
names(brca1_expression_avg) <- colnames(exp_combined)

# Validation: Ensure sample order consistency
stopifnot(identical(rownames(meta_combined), names(brca1_expression_avg)))

# Classify samples based on BRCA1 expression and tissue type
# Classification logic:
# - Control: Normal tissue (biological baseline)
# - upBRCA1: GAIN samples with expression ≥ upper threshold
# - downBRCA1: GAIN samples with expression ≤ lower threshold
# - NA: Samples with intermediate expression (excluded from contrast)
cat("Classifying samples by BRCA1 expression...\n")

meta_combined$BRCA1_status <- ifelse(
  meta_combined$sample_type == "Solid Tissue Normal", 
  "Control",
  ifelse(
    rownames(meta_combined) %in% gain_samples,
    ifelse(
      brca1_expression_avg >= up_threshold, "upBRCA1",
      ifelse(brca1_expression_avg <= down_threshold, "downBRCA1", NA)
    ),
    NA
  )
)

# Convert to ordered factor
meta_combined$BRCA1_status <- factor(
  meta_combined$BRCA1_status,
  levels = c("Control", "downBRCA1", "upBRCA1")
)

# Display classification summary
classification_table <- table(meta_combined$BRCA1_status, useNA = "ifany")
print(classification_table)

# ==============================================================================
# SECTION 9: EXPRESSION MATRIX FILTERING AND CLEANING
# ==============================================================================

# Align expression matrix with classified samples
filtered_matrix_brca1_exp <- exp_combined[, rownames(meta_combined)]

# Remove samples with NA BRCA1 status
# Rationale: Analysis focuses on samples with clear BRCA1 classification
meta_combined_exp_clean <- meta_combined[!is.na(meta_combined$BRCA1_status), ]
filtered_matrix_brca1_exp <- filtered_matrix_brca1_exp[, rownames(meta_combined_exp_clean)]

# Recalculate variance and remove zero-variance genes
# Rationale: Sample filtering may create new zero-variance genes
gene_variances <- apply(filtered_matrix_brca1_exp, 1, var)
filtered_matrix_brca1_exp <- filtered_matrix_brca1_exp[gene_variances > 0, ]

cat("Genes after variance filtering:", nrow(filtered_matrix_brca1_exp), "\n")
cat("Final expression matrix dimensions:", 
    paste(dim(filtered_matrix_brca1_exp), collapse = " x "), "\n")
cat("Clean metadata dimensions:", 
    paste(dim(meta_combined_exp_clean), collapse = " x "), "\n")

# ==============================================================================
# SECTION 10: GENE ANNOTATION - ENSEMBL TO SYMBOL MAPPING
# ==============================================================================

# Clean version numbers from Ensembl IDs
ensembl_ids_brca1variant <- gsub("\\..*", "", rownames(filtered_matrix_brca1_exp))

# Create mapping data frame
matrix_df_brca1variant <- data.frame(
  ensembl_id = ensembl_ids_brca1variant,
  row_index = 1:nrow(filtered_matrix_brca1_exp),
  stringsAsFactors = FALSE
)

# Merge with gene annotation
matrix_gene_info_brca1variant <- merge(
  matrix_df_brca1variant,
  gene_info_final_3.5[, c("ensembl_id", "gene_symbol")],
  by = "ensembl_id",
  all.x = TRUE
)

# Quality control
missing_ensembl <- sum(is.na(matrix_gene_info_brca1variant$ensembl_id) | 
                       matrix_gene_info_brca1variant$ensembl_id == "")
missing_symbol <- sum(is.na(matrix_gene_info_brca1variant$gene_symbol) | 
                      matrix_gene_info_brca1variant$gene_symbol == "")

cat("Missing Ensembl IDs:", missing_ensembl, "\n")
cat("Missing Gene Symbols:", missing_symbol, "\n")

# Identify unannotated genes
ensembl_without_symbol <- matrix_gene_info_brca1variant[
  !is.na(matrix_gene_info_brca1variant$ensembl_id) & 
  matrix_gene_info_brca1variant$ensembl_id != "" & 
  (is.na(matrix_gene_info_brca1variant$gene_symbol) | 
   matrix_gene_info_brca1variant$gene_symbol == ""), 
]

cat("Ensembl IDs without gene symbol:", nrow(ensembl_without_symbol), "\n")

# Filter valid annotations
matrix_gene_info_brca1variant <- matrix_gene_info_brca1variant[
  !(is.na(matrix_gene_info_brca1variant$gene_symbol) | 
    matrix_gene_info_brca1variant$gene_symbol == ""),
]

cat("Genes after annotation filtering:", nrow(matrix_gene_info_brca1variant), "\n")

# ==============================================================================
# SECTION 11: REPRESENTATIVE GENE SELECTION
# ==============================================================================

# Convert to data.table
filtered_matrix_dt_brca1 <- data.table(
  ensembl_id = matrix_gene_info_brca1variant$ensembl_id,
  gene_symbol = matrix_gene_info_brca1variant$gene_symbol,
  row_index = matrix_gene_info_brca1variant$row_index,
  filtered_matrix_brca1_exp[matrix_gene_info_brca1variant$row_index, ]
)

# Identify sample columns
sample_cols_exp_brca1 <- names(filtered_matrix_dt_brca1)[
  !names(filtered_matrix_dt_brca1) %in% c("ensembl_id", "gene_symbol", "row_index")
]

# Select representative Ensembl ID per gene based on highest MAD
# Methodology: Same as previous scripts - MAD-based selection for robustness
representative_genes_brca1_2 <- filtered_matrix_dt_brca1[, {
  if(length(sample_cols_exp_brca1) > 0) {
    gene_vars <- apply(.SD[, sample_cols_exp_brca1, with = FALSE], 1, function(x) {
      mad(as.numeric(x), na.rm = TRUE, constant = 1.4826)
    })
    .SD[which.max(gene_vars)]
  }
}, by = gene_symbol]

cat("Genes before duplicate removal:", nrow(filtered_matrix_dt_brca1), "\n")
cat("Genes after selecting representatives:", nrow(representative_genes_brca1_2), "\n")

# ==============================================================================
# SECTION 12: FINAL EXPRESSION MATRIX PREPARATION
# ==============================================================================

# Create final expression matrix
filtered_matrix_brca1_exp_final <- as.matrix(
  representative_genes_brca1_2[, sample_cols_exp_brca1, with = FALSE]
)
rownames(filtered_matrix_brca1_exp_final) <- representative_genes_brca1_2$ensembl_id

cat("Final expression matrix dimensions:", 
    paste(dim(filtered_matrix_brca1_exp_final), collapse = " x "), "\n")

# ==============================================================================
# SECTION 13: STATISTICAL DESIGN MATRIX
# ==============================================================================

# Create design matrix based on BRCA1 status
design_brca1 <- model.matrix(~ 0 + BRCA1_status, data = meta_combined_exp_clean)
colnames(design_brca1) <- levels(meta_combined_exp_clean$BRCA1_status)

cat("Design matrix columns:\n")
print(colnames(design_brca1))

# ==============================================================================
# SECTION 14: LINEAR MODEL FITTING
# ==============================================================================

# Fit linear model with limma
fit_brca1 <- lmFit(filtered_matrix_brca1_exp_final, design_brca1)

cat("Genes in model:", nrow(fit_brca1$coefficients), "\n")
cat("Coefficients:", paste(colnames(fit_brca1$coefficients), collapse = ", "), "\n")

# ==============================================================================
# SECTION 15: CONTRAST DEFINITIONS
# ==============================================================================

# Define contrasts for BRCA1-stratified comparisons:
# 1. upBRCA1 vs downBRCA1: Direct comparison of expression extremes
#    (isolates BRCA1 expression-associated effects)
# 2. upBRCA1 vs Control: High BRCA1 tumors vs normal
# 3. downBRCA1 vs Control: Low BRCA1 tumors vs normal
contrast_matrix_brca1 <- makeContrasts(
  "upBRCA1_vs_downBRCA1" = upBRCA1 - downBRCA1,   
  "upBRCA1_vs_Control" = upBRCA1 - Control,    
  "downBRCA1_vs_Control" = downBRCA1 - Control, 
  levels = design_brca1
)

cat("Contrast matrix:\n")
print(contrast_matrix_brca1)

# ==============================================================================
# SECTION 16: CONTRAST APPLICATION AND EMPIRICAL BAYES
# ==============================================================================

# Apply contrasts
fit_brca1 <- contrasts.fit(fit_brca1, contrast_matrix_brca1)

# Empirical Bayes moderation
fit_brca1 <- eBayes(fit_brca1)

# ==============================================================================
# SECTION 17: DEG RESULTS EXTRACTION
# ==============================================================================

# Initialize results list
deg_results_brca1_list <- list()

# Extract results for each contrast
for (contrast_name in colnames(contrast_matrix_brca1)) {
  
  cat("Processing contrast:", contrast_name, "\n")
  
  deg_result <- topTable(fit_brca1, coef = contrast_name, 
                         adjust = "fdr", number = Inf)
  
  deg_result$contrast <- contrast_name
  deg_result$ensembl_gene_id <- rownames(deg_result)
  
  deg_results_brca1_list[[contrast_name]] <- deg_result
}

# Combine results
deg_results_brca1_combined <- do.call(rbind, deg_results_brca1_list)

cat("Total DEG results:", nrow(deg_results_brca1_combined), "\n")
cat("BRCA1 present:", "BRCA1" %in% deg_results_brca1_combined$gene_symbol, "\n")

# ==============================================================================
# SECTION 18: GENE ANNOTATION ENRICHMENT
# ==============================================================================

# Merge with complete gene annotation
deg_results_brca1_combined <- merge(
  deg_results_brca1_combined,
  gene_info_final_3.5[, c("ensembl_id", "gene_symbol", "entrez_id", "chromosome")],
  by.x = "ensembl_gene_id",
  by.y = "ensembl_id",
  all.x = TRUE
)

# Quality control
missing_ensembl_final <- sum(is.na(deg_results_brca1_combined$ensembl_gene_id) | 
                             deg_results_brca1_combined$ensembl_gene_id == "")
missing_symbol_final <- sum(is.na(deg_results_brca1_combined$gene_symbol) | 
                            deg_results_brca1_combined$gene_symbol == "")

cat("Final missing Ensembl IDs:", missing_ensembl_final, "\n")
cat("Final missing Gene Symbols:", missing_symbol_final, "\n")

# Create unique identifier
deg_results_brca1_combined$unique_id <- paste(
  deg_results_brca1_combined$gene_symbol,
  deg_results_brca1_combined$contrast,
  sep = "_"
)

# Verify uniqueness
is_unique <- length(unique(deg_results_brca1_combined$unique_id)) == 
             nrow(deg_results_brca1_combined)
cat("All IDs unique:", is_unique, "\n")

rownames(deg_results_brca1_combined) <- deg_results_brca1_combined$unique_id

# ==============================================================================
# SECTION 19: SIGNIFICANT GENE IDENTIFICATION
# ==============================================================================

# Filter significant genes: |logFC| > 1 and FDR < 0.05
significant_brca1_combined <- deg_results_brca1_combined[
  abs(deg_results_brca1_combined$logFC) > 1 & 
  deg_results_brca1_combined$adj.P.Val < 0.05, 
]

# Add regulation status
significant_brca1_combined$status <- ifelse(
  significant_brca1_combined$logFC > 1, "UP_REG", 
  ifelse(significant_brca1_combined$logFC < -1, "DOWN_REG", "UNCHANGED")
)

cat("Total significant genes:", nrow(significant_brca1_combined), "\n")
cat("BRCA1 in significant genes:", 
    "BRCA1" %in% significant_brca1_combined$gene_symbol, "\n")

# Chr17-specific genes
chr17_deg_all_contrasts <- significant_brca1_combined[
  significant_brca1_combined$chromosome == "17", 
]
cat("Chr17 significant genes (all contrasts):", 
    nrow(chr17_deg_all_contrasts), "\n")

# ==============================================================================
# SECTION 20: PRIMARY CONTRAST - upBRCA1 VS downBRCA1
# ==============================================================================

# Primary contrast: upBRCA1 vs downBRCA1
up_down_brca1_deg <- significant_brca1_combined[
  significant_brca1_combined$contrast == "upBRCA1_vs_downBRCA1", 
]

cat("Significant genes in upBRCA1 vs downBRCA1:", nrow(up_down_brca1_deg), "\n")

# Chr17-specific genes
chr17_up_down_brca1_deg <- up_down_brca1_deg[up_down_brca1_deg$chromosome == "17", ]
cat("Chr17 significant genes (upBRCA1 vs downBRCA1):", 
    nrow(chr17_up_down_brca1_deg), "\n")
cat("BRCA1 in Chr17 primary contrast:", 
    "BRCA1" %in% chr17_up_down_brca1_deg$gene_symbol, "\n")

# ==============================================================================
# SECTION 21: SECONDARY CONTRAST - upBRCA1 VS CONTROL
# ==============================================================================

up_ctrl_brca1_deg <- significant_brca1_combined[
  significant_brca1_combined$contrast == "upBRCA1_vs_Control", 
]

cat("Significant genes in upBRCA1 vs Control:", nrow(up_ctrl_brca1_deg), "\n")

chr17_up_ctrl_brca1_deg <- up_ctrl_brca1_deg[up_ctrl_brca1_deg$chromosome == "17", ]
cat("Chr17 significant genes (upBRCA1 vs Control):", 
    nrow(chr17_up_ctrl_brca1_deg), "\n")
cat("BRCA1 in Chr17 upBRCA1 vs Control:", 
    "BRCA1" %in% chr17_up_ctrl_brca1_deg$gene_symbol, "\n")

# ==============================================================================
# SECTION 22: TERTIARY CONTRAST - downBRCA1 VS CONTROL
# ==============================================================================

down_ctrl_brca1_deg <- significant_brca1_combined[
  significant_brca1_combined$contrast == "downBRCA1_vs_Control", 
]

cat("Significant genes in downBRCA1 vs Control:", nrow(down_ctrl_brca1_deg), "\n")

chr17_down_ctrl_brca1_deg <- down_ctrl_brca1_deg[
  down_ctrl_brca1_deg$chromosome == "17", 
]
cat("Chr17 significant genes (downBRCA1 vs Control):", 
    nrow(chr17_down_ctrl_brca1_deg), "\n")
cat("BRCA1 in Chr17 downBRCA1 vs Control:", 
    "BRCA1" %in% chr17_down_ctrl_brca1_deg$gene_symbol, "\n")

# ==============================================================================
# SECTION 23: VISUALIZATION - VOLCANO PLOT
# ==============================================================================

# Prepare data for volcano plot (primary contrast)
deg_up_down_brca1 <- deg_results_brca1_combined[
  deg_results_brca1_combined$contrast == "upBRCA1_vs_downBRCA1", 
]

# Filter Chr17 genes
deg_chr17q_up_down_brca1 <- deg_results_brca1_combined[
  deg_results_brca1_combined$contrast == "upBRCA1_vs_downBRCA1" &
  deg_results_brca1_combined$chromosome == "17",
]

# Generate volcano plot
volcano_plot <- EnhancedVolcano(
  toptable = deg_chr17q_up_down_brca1,            
  lab = deg_chr17q_up_down_brca1$gene_symbol,    
  x = 'logFC',                    
  y = 'adj.P.Val',              
  title = 'DEG Analysis: upBRCA1 vs downBRCA1',
  xlab = bquote(~Log[2] ~ "Fold Change"),
  ylab = bquote(~-Log[10] ~ "Adjusted P-value"),
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,          
  labSize = 5.0,
  colAlpha = 0.7,
  legendPosition = 'none'
)

print(volcano_plot)

# Save plot
# pdf("results/volcano_plot_upbrca1_vs_downbrca1.pdf", width = 10, height = 8)
# print(volcano_plot)
# dev.off()

# ==============================================================================
# SECTION 24: RESULTS SUMMARY
# ==============================================================================

cat("\n=== BRCA1-STRATIFIED EXPRESSION ANALYSIS SUMMARY ===\n")
cat("Total genes analyzed:", nrow(deg_results_brca1_combined)/3, "\n")
cat("Significant genes (all contrasts):", nrow(significant_brca1_combined), "\n")
cat("Significant - upBRCA1 vs downBRCA1:", 
    sum(significant_brca1_combined$contrast == "upBRCA1_vs_downBRCA1"), "\n")
cat("Significant - upBRCA1 vs Control:", 
    sum(significant_brca1_combined$contrast == "upBRCA1_vs_Control"), "\n")
cat("Significant - downBRCA1 vs Control:", 
    sum(significant_brca1_combined$contrast == "downBRCA1_vs_Control"), "\n")
cat("Chr17 significant (upBRCA1 vs downBRCA1):", 
    nrow(chr17_up_down_brca1_deg), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================