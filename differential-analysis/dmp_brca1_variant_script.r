# ==============================================================================
# DIFFERENTIAL METHYLATION ANALYSIS - BRCA1 Expression Stratification
# ==============================================================================
#
# Title: DNA Methylation Analysis Stratified by BRCA1 Expression Levels
#
# Description: 
#   This script stratifies Chr17q GAIN tumor samples based on BRCA1 expression
#   levels and performs differential methylation analysis. Samples are classified
#   into upBRCA1 and downBRCA1 groups using the same thresholds from expression
#   analysis. Identifies CpG probes with methylation changes associated with
#   BRCA1 expression status, revealing epigenetic patterns linked to BRCA1.
#
# Author: [Your Name]
# Date: 2025-01-XX
# Version: 1.0
#
# Input Files:
#   - data/m_met_combined_matrix.rds      : M-value methylation matrix
#   - data/meta_data_chr17q_met.rds       : Methylation metadata with Chr17q status
#   - data/annot_one2one_2.rds            : Probe annotation
#   - data/exp_combined_2.rds             : Expression data (for BRCA1 classification)
#
# Output Files:
#   - results/dmp_brca1_combined.csv              : All DMP results
#   - results/dmp_brca1_significant.csv           : Significant DMPs
#   - results/dmp_upbrca1_vs_downbrca1.csv        : Primary contrast
#   - results/volcano_plot_met_upbrca1_vs_downbrca1.pdf : Volcano plot
#
# Dependencies:
#   - data.table (>= 1.14.0)
#   - dplyr (>= 1.0.0)
#   - tibble (>= 3.0.0)
#   - limma (>= 3.48.0)
#   - EnhancedVolcano (>= 1.10.0)
#
# ==============================================================================

# ==============================================================================
# LIBRARY IMPORTS
# ==============================================================================

library(data.table)
library(dplyr)
library(tibble)
library(limma)
library(EnhancedVolcano)

# ==============================================================================
# SECTION 1: METHYLATION DATA PREPROCESSING
# ==============================================================================

# Calculate variance per probe and remove zero-variance probes
variances_chr17q_met <- apply(m_met_combined_matrix, 1, var)
filtered_matrix_chr17q_met <- m_met_combined_matrix[variances_chr17q_met > 0, ]

cat("Probes before filtering:", nrow(m_met_combined_matrix), "\n")
cat("Probes after filtering:", nrow(filtered_matrix_chr17q_met), "\n")

# ==============================================================================
# SECTION 2: METADATA PREPROCESSING
# ==============================================================================

# Filter metadata for complete Chr17q status
filtered_meta_data_chr17q_met <- meta_data_chr17q_met[
  meta_data_chr17q_met$Chr17q_Status_Final %in% c("Control", "DIS", "GAIN"), 
]
rownames(filtered_meta_data_chr17q_met) <- filtered_meta_data_chr17q_met$barcode

# Align methylation matrix with metadata
filtered_matrix_chr17q_met <- filtered_matrix_chr17q_met[
  , rownames(filtered_meta_data_chr17q_met)
]

cat("Samples after filtering:", ncol(filtered_matrix_chr17q_met), "\n")

# ==============================================================================
# SECTION 3: SAMPLE IDENTIFICATION BY STATUS
# ==============================================================================

# Extract GAIN samples for stratification
gain_samples_met <- rownames(filtered_meta_data_chr17q_met)[
  filtered_meta_data_chr17q_met$Chr17q_Status_Final == "GAIN"
]
cat("GAIN samples:", length(gain_samples_met), "\n")

# Extract normal control samples
normal_samples_met <- rownames(filtered_meta_data_chr17q_met)[
  filtered_meta_data_chr17q_met$sample_type == "Solid Tissue Normal"
]
cat("Control samples:", length(normal_samples_met), "\n")

# ==============================================================================
# SECTION 4: COMBINE GAIN AND CONTROL SAMPLES
# ==============================================================================

# Merge metadata
meta_gain_met <- filtered_meta_data_chr17q_met[gain_samples_met, ]
meta_normal_met <- filtered_meta_data_chr17q_met[normal_samples_met, ]
meta_combined_met <- rbind(meta_gain_met, meta_normal_met)

# Extract methylation for combined samples
combined_samples_met <- rownames(meta_combined_met)
met_combined <- filtered_matrix_chr17q_met[, combined_samples_met]

cat("Total samples in analysis:", length(combined_samples_met), "\n")

# ==============================================================================
# SECTION 5: BRCA1 METHYLATION EXTRACTION
# ==============================================================================

# Convert methylation matrix to data frame with probe IDs
met_matrix_combined_df <- as.data.frame(met_combined) %>%
  tibble::rownames_to_column(var = "probeID")

# Join with probe annotation to get gene information
met_matrix_with_gene <- met_matrix_combined_df %>%
  left_join(annot_one2one_2 %>% select(probeID, gene_name), by = "probeID")

# Filter BRCA1-specific probes
# Scientific rationale: Multiple CpG probes map to BRCA1 gene regions,
# capturing methylation patterns across promoter, gene body, and regulatory regions
brca1_met_matrix <- met_matrix_with_gene %>%
  filter(gene_name == "BRCA1")

cat("BRCA1 methylation probes found:", nrow(brca1_met_matrix), "\n")

# Extract numeric methylation data
brca1_met_matrix_numeric <- brca1_met_matrix %>%
  select(-gene_name, -probeID) %>%
  as.matrix()

# Calculate average methylation across all BRCA1 probes per sample
# Methodology: Averaging across multiple probes provides a robust measure
# of overall BRCA1 methylation status, reducing single-probe noise
brca1_methylation_avg <- colMeans(brca1_met_matrix_numeric, na.rm = TRUE)

cat("Extracting BRCA1 methylation from", nrow(brca1_met_matrix), "probe(s)...\n")

# ==============================================================================
# SECTION 6: NORMALITY TEST OF BRCA1 METHYLATION
# ==============================================================================

# Test normality on GAIN samples only
cat("Performing normality test on GAIN samples...\n")

brca1_methylation_gain_only <- brca1_methylation_avg[
  names(brca1_methylation_avg) %in% gain_samples_met
]
brca1_methylation_gain_vector <- as.numeric(brca1_methylation_gain_only)

# Kolmogorov-Smirnov test
ks_result_met <- ks.test(
  brca1_methylation_gain_vector,
  "pnorm",
  mean = mean(brca1_methylation_gain_vector, na.rm = TRUE),
  sd = sd(brca1_methylation_gain_vector, na.rm = TRUE)
)

cat("KS test p-value for BRCA1 methylation (GAIN samples):", 
    ks_result_met$p.value, "\n")
if (ks_result_met$p.value > 0.05) {
  cat("Data follows normal distribution (p > 0.05)\n")
} else {
  cat("Data does not follow normal distribution (p <= 0.05)\n")
}

# ==============================================================================
# SECTION 7: BRCA1 METHYLATION THRESHOLD DETERMINATION
# ==============================================================================

# Calculate thresholds: mean ± 1 SD
# Methodology note: These thresholds are based on METHYLATION levels,
# not expression. High methylation (hyperBRCA1) typically correlates with
# low expression, and vice versa, due to methylation's gene-silencing role.
cat("Calculating BRCA1 methylation thresholds...\n")

mean_brca1_methyl <- mean(brca1_methylation_gain_vector, na.rm = TRUE)
sd_brca1_methyl <- sd(brca1_methylation_gain_vector, na.rm = TRUE)

up_threshold_met <- mean_brca1_methyl + 1.0 * sd_brca1_methyl
down_threshold_met <- mean_brca1_methyl - 1.0 * sd_brca1_methyl

cat("Mean BRCA1 methylation (GAIN samples):", round(mean_brca1_methyl, 4), "\n")
cat("Standard deviation:", round(sd_brca1_methyl, 4), "\n")
cat("Upper threshold (hyperBRCA1):", round(up_threshold_met, 4), "\n")
cat("Lower threshold (hypoBRCA1):", round(down_threshold_met, 4), "\n")

# ==============================================================================
# SECTION 8: SAMPLE CLASSIFICATION BY BRCA1 METHYLATION
# ==============================================================================

# Validate sample order consistency
stopifnot(identical(rownames(meta_combined_met), names(brca1_methylation_avg)))

# Classify samples based on BRCA1 methylation
# Classification note: Using same labels (upBRCA1/downBRCA1) as expression
# analysis for consistency, though biological interpretation differs:
# - upBRCA1 methylation = hypermethylated BRCA1 region
# - downBRCA1 methylation = hypomethylated BRCA1 region
cat("Classifying samples by BRCA1 methylation...\n")

meta_combined_met$BRCA1_methylation_status <- ifelse(
  meta_combined_met$sample_type == "Solid Tissue Normal", 
  "Control",
  ifelse(
    rownames(meta_combined_met) %in% gain_samples_met,
    ifelse(
      brca1_methylation_avg >= up_threshold_met, "upBRCA1",
      ifelse(brca1_methylation_avg <= down_threshold_met, "downBRCA1", NA)
    ),
    NA
  )
)

# Convert to factor
meta_combined_met$BRCA1_methylation_status <- factor(
  meta_combined_met$BRCA1_methylation_status,
  levels = c("Control", "downBRCA1", "upBRCA1")
)

# Display classification
classification_summary <- table(meta_combined_met$BRCA1_methylation_status, useNA = "ifany")
print(classification_summary)

# ==============================================================================
# SECTION 9: METHYLATION MATRIX FILTERING
# ==============================================================================

# Align with classified samples
filtered_matrix_brca1_met <- met_combined[, rownames(meta_combined_met)]

# Remove samples with NA classification
meta_combined_met_clean <- meta_combined_met[
  !is.na(meta_combined_met$BRCA1_methylation_status), 
]
filtered_matrix_brca1_met <- filtered_matrix_brca1_met[
  , rownames(meta_combined_met_clean)
]

cat("Final methylation matrix dimensions:", 
    paste(dim(filtered_matrix_brca1_met), collapse = " x "), "\n")
cat("Clean metadata dimensions:", 
    paste(dim(meta_combined_met_clean), collapse = " x "), "\n")

# ==============================================================================
# SECTION 10: PROBE ANNOTATION MAPPING
# ==============================================================================

# Extract probe IDs
probe_ids <- rownames(filtered_matrix_brca1_met)

# Create mapping data frame
matrix_met_df <- data.frame(
  probeID = probe_ids,
  row_index = 1:nrow(filtered_matrix_brca1_met),
  stringsAsFactors = FALSE
)

# Filter invalid probe IDs
matrix_met_df <- matrix_met_df[
  !(is.na(matrix_met_df$probeID) | matrix_met_df$probeID == ""), 
]

cat("Probes after ID validation:", nrow(matrix_met_df), "\n")

# Merge with annotation
matrix_met_gene_info <- merge(
  matrix_met_df,
  annot_one2one_2[, c("probeID", "gene_name")],
  by = "probeID",
  all.x = TRUE
)

# Filter probes with valid gene annotation
matrix_met_gene_info <- matrix_met_gene_info[
  !(is.na(matrix_met_gene_info$gene_name) | matrix_met_gene_info$gene_name == ""), 
]

cat("Probes with valid annotation:", nrow(matrix_met_gene_info), "\n")

# ==============================================================================
# SECTION 11: REPRESENTATIVE PROBE SELECTION
# ==============================================================================

# Convert to data.table
filtered_matrix_met_dt <- data.table(
  probeID = matrix_met_gene_info$probeID,
  gene_name = matrix_met_gene_info$gene_name,
  row_index = matrix_met_gene_info$row_index,
  filtered_matrix_brca1_met[matrix_met_gene_info$row_index, ]
)

# Identify sample columns
sample_cols <- names(filtered_matrix_met_dt)[
  !names(filtered_matrix_met_dt) %in% c("probeID", "gene_name", "row_index")
]

# Select representative probe per gene based on MAD
representative_probes_2 <- filtered_matrix_met_dt[, {
  if (length(sample_cols) > 0) {
    probe_vars <- apply(.SD[, sample_cols, with = FALSE], 1, function(x) {
      mad(as.numeric(x), na.rm = TRUE, constant = 1.4826)
    })
    .SD[which.max(probe_vars)]
  }
}, by = gene_name]

cat("Probes before duplicate removal:", nrow(filtered_matrix_met_dt), "\n")
cat("Probes after selecting representatives:", nrow(representative_probes_2), "\n")

# ==============================================================================
# SECTION 12: FINAL METHYLATION MATRIX CONSTRUCTION
# ==============================================================================

# Create final matrix
filtered_matrix_brca1_met_final <- as.matrix(
  representative_probes_2[, sample_cols, with = FALSE]
)
rownames(filtered_matrix_brca1_met_final) <- representative_probes_2$probeID

cat("Final methylation matrix dimensions:", 
    paste(dim(filtered_matrix_brca1_met_final), collapse = " x "), "\n")

# ==============================================================================
# SECTION 13: STATISTICAL DESIGN MATRIX
# ==============================================================================

# Create design matrix
design_brca1_met <- model.matrix(
  ~ 0 + BRCA1_methylation_status, 
  data = meta_combined_met_clean
)
colnames(design_brca1_met) <- levels(meta_combined_met_clean$BRCA1_methylation_status)

cat("Design matrix columns:\n")
print(colnames(design_brca1_met))

# ==============================================================================
# SECTION 14: LINEAR MODEL FITTING
# ==============================================================================

# Fit linear model
fit_brca1_met <- lmFit(filtered_matrix_brca1_met_final, design_brca1_met)

cat("Probes in model:", nrow(fit_brca1_met$coefficients), "\n")
cat("Coefficients:", colnames(fit_brca1_met$coefficients), "\n")

# ==============================================================================
# SECTION 15: CONTRAST DEFINITIONS
# ==============================================================================

# Define contrasts for BRCA1-stratified methylation comparisons
# Scientific interpretation:
# 1. upBRCA1 vs downBRCA1: Compares hyper- vs hypo-methylated BRCA1 regions
# 2. upBRCA1 vs Control: Hypermethylated tumors vs normal baseline
# 3. downBRCA1 vs Control: Hypomethylated tumors vs normal baseline
contrast_matrix_brca1_met <- makeContrasts(
  "upBRCA1_vs_downBRCA1" = upBRCA1 - downBRCA1,   
  "upBRCA1_vs_Control" = upBRCA1 - Control,       
  "downBRCA1_vs_Control" = downBRCA1 - Control, 
  levels = design_brca1_met
)

cat("Contrast matrix:\n")
print(contrast_matrix_brca1_met)

# ==============================================================================
# SECTION 16: CONTRAST APPLICATION AND EMPIRICAL BAYES
# ==============================================================================

# Apply contrasts
fit_brca1_met <- contrasts.fit(fit_brca1_met, contrast_matrix_brca1_met)

# Empirical Bayes moderation
fit_brca1_met <- eBayes(fit_brca1_met)

# ==============================================================================
# SECTION 17: DMP RESULTS EXTRACTION
# ==============================================================================

# Initialize results list
dmp_results_brca1_list <- list()

# Extract results for each contrast
for (contrast_name in colnames(contrast_matrix_brca1_met)) {
  
  dmp_result <- topTable(
    fit_brca1_met, 
    coef = contrast_name, 
    adjust = "fdr", 
    number = Inf
  )
  
  dmp_result$contrast <- contrast_name
  dmp_result$probe_id <- rownames(dmp_result)
  
  dmp_results_brca1_list[[contrast_name]] <- dmp_result
}

# Combine results
dmp_results_brca1_combined <- do.call(rbind, dmp_results_brca1_list)

cat("Total DMP results:", nrow(dmp_results_brca1_combined), "\n")

# ==============================================================================
# SECTION 18: PROBE ANNOTATION ENRICHMENT
# ==============================================================================

# Merge with complete annotation
dmp_results_brca1_combined <- merge(
  dmp_results_brca1_combined, 
  annot_one2one_2[, c("probeID", "gene_name", "chromosome_name")],
  by.x = "probe_id", 
  by.y = "probeID", 
  all.x = TRUE
)

# Quality control
cat("Probes without gene name:", 
    sum(is.na(dmp_results_brca1_combined$gene_name) | 
        dmp_results_brca1_combined$gene_name == ""), "\n")
cat("Probes without chromosome:", 
    sum(is.na(dmp_results_brca1_combined$chromosome_name) | 
        dmp_results_brca1_combined$chromosome_name == ""), "\n")

# ==============================================================================
# SECTION 19: UNIQUENESS VERIFICATION
# ==============================================================================

# Check for multiple probes per gene
multiple_probes_check <- dmp_results_brca1_combined %>%
  distinct(gene_name, probe_id) %>%
  dplyr::count(gene_name) %>%
  filter(n > 1)

cat("Genes with multiple probes:", nrow(multiple_probes_check), "\n")

# Check for multiple gene annotations
multiple_gene_check <- table(
  multiple = grepl(";|,|\\|", dmp_results_brca1_combined$gene_name)
)
cat("Probe distribution by gene annotation:\n")
print(multiple_gene_check)

# Create unique identifier
dmp_results_brca1_combined$unique_id <- paste(
  dmp_results_brca1_combined$gene_name,
  dmp_results_brca1_combined$contrast,
  sep = "_"
)

# Verify uniqueness
unique_check <- length(unique(dmp_results_brca1_combined$unique_id)) == 
                nrow(dmp_results_brca1_combined)
cat("Unique ID verification:", ifelse(unique_check, "PASSED", "FAILED"), "\n")

rownames(dmp_results_brca1_combined) <- dmp_results_brca1_combined$unique_id

cat("BRCA1 present in results:", 
    "BRCA1" %in% dmp_results_brca1_combined$gene_name, "\n")

# ==============================================================================
# SECTION 20: SIGNIFICANT PROBE FILTERING
# ==============================================================================

# Filter significant probes: |logFC| > 0.2 and FDR < 0.05
# Methodology: Threshold reflects methylation-specific effect sizes
significant_brca1_met_combined <- dmp_results_brca1_combined[
  abs(dmp_results_brca1_combined$logFC) > 0.2 & 
  dmp_results_brca1_combined$adj.P.Val < 0.05, 
]

# Create methylation status labels
significant_brca1_met_combined$status <- ifelse(
  significant_brca1_met_combined$logFC > 0.2, "HYPERMETH", 
  ifelse(significant_brca1_met_combined$logFC < -0.2, "HYPOMETH", "UNCHANGED")
)

cat("Total significant probes:", nrow(significant_brca1_met_combined), "\n")
cat("BRCA1 in significant probes:", 
    "BRCA1" %in% significant_brca1_met_combined$gene_name, "\n")

# Chr17-specific probes
chr17_dmp_all_contrasts <- significant_brca1_met_combined[
  significant_brca1_met_combined$chromosome_name == "17", 
]
cat("Chr17 significant probes (all contrasts):", 
    nrow(chr17_dmp_all_contrasts), "\n")

# ==============================================================================
# SECTION 21: PRIMARY CONTRAST - upBRCA1 VS downBRCA1
# ==============================================================================

# Primary contrast
upbrca1_vs_downbrca1_dmp <- significant_brca1_met_combined[
  significant_brca1_met_combined$contrast == "upBRCA1_vs_downBRCA1", 
]

cat("Significant probes in upBRCA1 vs downBRCA1:", 
    nrow(upbrca1_vs_downbrca1_dmp), "\n")

# Chr17-specific
chr17_upbrca1_vs_downbrca1_dmp <- upbrca1_vs_downbrca1_dmp[
  upbrca1_vs_downbrca1_dmp$chromosome_name == "17", 
]

cat("Chr17 significant probes (upBRCA1 vs downBRCA1):", 
    nrow(chr17_upbrca1_vs_downbrca1_dmp), "\n")
cat("BRCA1 in Chr17 primary contrast:", 
    "BRCA1" %in% chr17_upbrca1_vs_downbrca1_dmp$gene_name, "\n")

# ==============================================================================
# SECTION 22: SECONDARY CONTRAST - upBRCA1 VS CONTROL
# ==============================================================================

upbrca1_vs_ctrl_dmp <- significant_brca1_met_combined[
  significant_brca1_met_combined$contrast == "upBRCA1_vs_Control", 
]

cat("Significant probes in upBRCA1 vs Control:", 
    nrow(upbrca1_vs_ctrl_dmp), "\n")

chr17_upbrca1_vs_ctrl_dmp <- upbrca1_vs_ctrl_dmp[
  upbrca1_vs_ctrl_dmp$chromosome_name == "17", 
]

cat("Chr17 significant probes (upBRCA1 vs Control):", 
    nrow(chr17_upbrca1_vs_ctrl_dmp), "\n")
cat("BRCA1 in Chr17 upBRCA1 vs Control:", 
    "BRCA1" %in% chr17_upbrca1_vs_ctrl_dmp$gene_name, "\n")

# ==============================================================================
# SECTION 23: TERTIARY CONTRAST - downBRCA1 VS CONTROL
# ==============================================================================

downbrca1_vs_ctrl_dmp <- significant_brca1_met_combined[
  significant_brca1_met_combined$contrast == "downBRCA1_vs_Control", 
]

cat("Significant probes in downBRCA1 vs Control:", 
    nrow(downbrca1_vs_ctrl_dmp), "\n")

chr17_downbrca1_vs_ctrl_dmp <- downbrca1_vs_ctrl_dmp[
  downbrca1_vs_ctrl_dmp$chromosome_name == "17", 
]

cat("Chr17 significant probes (downBRCA1 vs Control):", 
    nrow(chr17_downbrca1_vs_ctrl_dmp), "\n")

# ==============================================================================
# SECTION 24: VISUALIZATION - VOLCANO PLOT
# ==============================================================================

# Prepare data for volcano plot
dmp_upbrca1_vs_downbrca1 <- dmp_results_brca1_combined[
  dmp_results_brca1_combined$contrast == "upBRCA1_vs_downBRCA1", 
]

# Filter Chr17 probes
dmp_chr17q_up_down_brca1 <- dmp_results_brca1_combined[
  dmp_results_brca1_combined$contrast == "upBRCA1_vs_downBRCA1" &
  dmp_results_brca1_combined$chromosome_name == "17",
]

# Create volcano plot
volcano_plot_brca1_met <- EnhancedVolcano(
  toptable = dmp_chr17q_up_down_brca1,            
  lab = dmp_chr17q_up_down_brca1$gene_name,    
  x = 'logFC',                    
  y = 'adj.P.Val',              
  title = 'DMP Analysis: upBRCA1 vs downBRCA1',
  xlab = bquote(~Log[2] ~ "Methylation Change"),
  ylab = bquote(~-Log[10] ~ "Adjusted P-value"),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  pointSize = 2.0,          
  labSize = 5.0,
  colAlpha = 0.7,
  legendPosition = 'none'
)

print(volcano_plot_brca1_met)

# Save plot
# pdf("results/volcano_plot_met_upbrca1_vs_downbrca1.pdf", width = 10, height = 8)
# print(volcano_plot_brca1_met)
# dev.off()

# ==============================================================================
# SECTION 25: RESULTS SUMMARY
# ==============================================================================

cat("\n=== BRCA1-STRATIFIED METHYLATION ANALYSIS SUMMARY ===\n")
cat("Total probes analyzed:", nrow(dmp_results_brca1_combined)/3, "\n")
cat("Significant probes (all contrasts):", nrow(significant_brca1_met_combined), "\n")
cat("Significant - upBRCA1 vs downBRCA1:", 
    sum(significant_brca1_met_combined$contrast == "upBRCA1_vs_downBRCA1"), "\n")
cat("Significant - upBRCA1 vs Control:", 
    sum(significant_brca1_met_combined$contrast == "upBRCA1_vs_Control"), "\n")
cat("Significant - downBRCA1 vs Control:", 
    sum(significant_brca1_met_combined$contrast == "downBRCA1_vs_Control"), "\n")
cat("Chr17 significant (upBRCA1 vs downBRCA1):", 
    nrow(chr17_upbrca1_vs_downbrca1_dmp), "\n")

# Analyze hyper/hypomethylation ratio
if(nrow(upbrca1_vs_downbrca1_dmp) > 0) {
  hyper_count <- sum(upbrca1_vs_downbrca1_dmp$status == "HYPERMETH")
  hypo_count <- sum(upbrca1_vs_downbrca1_dmp$status == "HYPOMETH")
  
  cat("\nHypermethylated probes (upBRCA1 vs downBRCA1):", hyper_count, "\n")
  cat("Hypomethylated probes (upBRCA1 vs downBRCA1):", hypo_count, "\n")
  cat("Hyper:Hypo ratio:", round(hyper_count/hypo_count, 2), "\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================