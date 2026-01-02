# ==============================================================================
# DIFFERENTIAL METHYLATION PROBE (DMP) ANALYSIS - Chr17q Copy Number Status
# ==============================================================================
#
# Title: Differential Methylation Analysis Based on Chr17q Copy Number Status
#
# Description: 
#   This script performs differential methylation analysis comparing tumor samples
#   with Chr17q copy number alterations (GAIN vs Disomic) against normal controls.
#   Uses M-values from Illumina 450K/EPIC arrays and limma for statistical testing.
#   Identifies CpG probes with significant methylation changes associated with
#   Chr17q copy number alterations in breast cancer.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Input Files:
#   - data/combined_barcodes_MET_sum.rds  : Combined methylation data (SummarizedExperiment)
#   - data/copy_number_prepared.rds       : Copy number status annotations
#   - data/annot_one2one_2.rds            : Probe annotation (450K/EPIC)
#   - data/m_met_combined_matrix.rds      : M-value matrix
#
# Output Files:
#   - results/dmp_chr17q_combined.csv          : All DMP results (3 contrasts)
#   - results/dmp_chr17q_significant.csv       : Significant DMPs
#   - results/dmp_gain_vs_dis.csv              : Primary contrast (GAIN vs DIS)
#   - results/volcano_plot_met_gain_vs_ctrl.pdf: Volcano plot visualization
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
# SECTION 1: METADATA INTEGRATION FOR METHYLATION DATA
# ==============================================================================

# Extract metadata from methylation SummarizedExperiment object
meta_data_met <- colData(combined_barcodes_MET_sum)
meta_data_met_df <- as.data.frame(meta_data_met)

# Create truncated sample IDs (15 characters) for matching
# Rationale: Standardize sample identifiers for integration with copy number data
meta_data_met_df$Sample_15 <- substr(meta_data_met_df$sample, 1, 15)

# Join with copy number annotations
# Methodology: Integrate genomic copy number status with epigenomic methylation data
meta_data_met_df <- meta_data_met_df %>%
  left_join(copy_number_prepared[, c("SAMPLE_ID_15", "Chr17q_Status")], 
            by = c("Sample_15" = "SAMPLE_ID_15"))

# Assign "Control" for normal tissue, keep GAIN/DIS for tumors
# Scientific rationale: Normal tissue represents unmutated epigenetic baseline
meta_data_met_df$Chr17q_Status_Final <- ifelse(
  meta_data_met_df$sample_type == "Solid Tissue Normal", 
  "Control",
  meta_data_met_df$Chr17q_Status
)

# Convert to ordered factor
meta_data_met_df$Chr17q_Status_Final <- factor(
  meta_data_met_df$Chr17q_Status_Final, 
  levels = c("Control", "DIS", "GAIN")
)

# Convert back to DFrame
meta_data_chr17q_met <- as(meta_data_met_df, "DFrame")

# Display sample distribution
cat("Methylation sample distribution by Chr17q Status:\n")
print(table(meta_data_chr17q_met$Chr17q_Status_Final))

# ==============================================================================
# SECTION 2: PROBE FILTERING BY VARIANCE
# ==============================================================================

# Calculate variance per probe on M-value data
# Methodology rationale: M-values are more suitable for statistical analysis
# than beta-values due to better homoscedasticity. M-value = log2(Beta/(1-Beta))
# provides better performance in differential methylation testing.
variances_chr17q_met <- apply(m_met_combined_matrix, 1, var)

# Remove probes with zero variance
# Rationale: Zero-variance probes provide no discriminatory power and can
# cause numerical issues in linear modeling
filtered_matrix_chr17q_met <- m_met_combined_matrix[variances_chr17q_met > 0, ]

cat("Probes before variance filtering:", nrow(m_met_combined_matrix), "\n")
cat("Probes after variance filtering:", nrow(filtered_matrix_chr17q_met), "\n")

# Filter metadata for samples with complete Chr17q status
filtered_meta_data_chr17q_met <- meta_data_chr17q_met[
  meta_data_chr17q_met$Chr17q_Status_Final %in% c("Control", "DIS", "GAIN"), 
]
rownames(filtered_meta_data_chr17q_met) <- filtered_meta_data_chr17q_met$barcode  

# Align methylation matrix with filtered metadata
filtered_matrix_chr17q_met <- filtered_matrix_chr17q_met[
  , rownames(filtered_meta_data_chr17q_met)
]

cat("Methylation matrix dimensions:", 
    paste(dim(filtered_matrix_chr17q_met), collapse = " x "), "\n")
cat("Metadata dimensions:", 
    paste(dim(filtered_meta_data_chr17q_met), collapse = " x "), "\n")
cat("Sample distribution (after filtering):\n")
print(table(filtered_meta_data_chr17q_met$Chr17q_Status_Final))

# ==============================================================================
# SECTION 3: PROBE ANNOTATION - MAPPING TO GENE NAMES
# ==============================================================================

# Extract probe IDs from methylation matrix
probe_ids <- rownames(filtered_matrix_chr17q_met)

# Create data frame for annotation merging
matrix_met_df <- data.frame(
  probeID = probe_ids,
  row_index = 1:nrow(filtered_matrix_chr17q_met),
  stringsAsFactors = FALSE
)

# Filter out invalid probe IDs
# Quality control: Remove probes with NA or empty IDs
matrix_met_df <- matrix_met_df[
  !(is.na(matrix_met_df$probeID) | matrix_met_df$probeID == ""), 
]

cat("Probes after ID validation:", nrow(matrix_met_df), "\n")

# Merge with probe annotation to get gene names
# Methodology: annot_one2one_2 contains curated 1-to-1 probe-to-gene mappings
# from Illumina array annotation, ensuring unambiguous gene assignment
matrix_met_gene_info <- merge(
  matrix_met_df,
  annot_one2one_2[, c("probeID", "gene_name")],
  by = "probeID",
  all.x = TRUE
)

# Filter probes without valid gene annotation
# Rationale: Unannotated probes cannot be biologically interpreted
matrix_met_gene_info <- matrix_met_gene_info[
  !(is.na(matrix_met_gene_info$gene_name) | matrix_met_gene_info$gene_name == ""), 
]

cat("Probes with valid gene annotation:", nrow(matrix_met_gene_info), "\n")

# ==============================================================================
# SECTION 4: REPRESENTATIVE PROBE SELECTION FOR DUPLICATE GENES
# ==============================================================================

# Convert to data.table for efficient processing
filtered_matrix_met_dt <- data.table(
  probeID = matrix_met_gene_info$probeID,
  gene_name = matrix_met_gene_info$gene_name,
  row_index = matrix_met_gene_info$row_index,
  filtered_matrix_chr17q_met[matrix_met_gene_info$row_index, ]
)

# Identify sample columns
sample_cols <- names(filtered_matrix_met_dt)[
  !names(filtered_matrix_met_dt) %in% c("probeID", "gene_name", "row_index")
]

# Select one representative probe per gene based on highest MAD
# Methodology: Multiple CpG probes often map to the same gene. We select
# the most variable probe using Median Absolute Deviation (MAD).
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
# SECTION 5: FINAL METHYLATION MATRIX CONSTRUCTION
# ==============================================================================

# Construct final methylation matrix with representative probes
filtered_matrix_chr17q_met_final <- as.matrix(
  representative_probes_2[, sample_cols, with = FALSE]
)
rownames(filtered_matrix_chr17q_met_final) <- representative_probes_2$probeID

cat("Final methylation matrix dimensions:", 
    paste(dim(filtered_matrix_chr17q_met_final), collapse = " x "), "\n")

# ==============================================================================
# SECTION 6: STATISTICAL DESIGN - LINEAR MODEL SETUP
# ==============================================================================

# Create design matrix without intercept
# Methodology: Allows direct comparison between all groups
design_chr17q_met <- model.matrix(
  ~ 0 + Chr17q_Status_Final, 
  data = filtered_meta_data_chr17q_met
)
colnames(design_chr17q_met) <- levels(filtered_meta_data_chr17q_met$Chr17q_Status_Final)

cat("Design matrix columns:\n")
print(colnames(design_chr17q_met))

# ==============================================================================
# SECTION 7: LINEAR MODEL FITTING WITH LIMMA
# ==============================================================================

# Fit linear model for each probe
# Methodology: Limma's empirical Bayes framework is well-suited for
# methylation array data, providing stable variance estimates
fit_chr17q_met <- lmFit(filtered_matrix_chr17q_met_final, design_chr17q_met)

cat("Probes in fitted model:", nrow(fit_chr17q_met$coefficients), "\n")
cat("Model coefficients:", colnames(fit_chr17q_met$coefficients), "\n")

# ==============================================================================
# SECTION 8: CONTRAST DEFINITIONS FOR METHYLATION COMPARISONS
# ==============================================================================

# Define three key contrasts:
# 1. GAIN vs DIS: Isolates effect of Chr17q gain on methylation
# 2. GAIN vs Control: Tumor with gain vs normal methylation baseline
# 3. DIS vs Control: Disomic tumor vs normal methylation
# Scientific rationale: These contrasts separate copy number-specific effects
# from general cancer-associated methylation changes
contrast_matrix_chr17q_met <- makeContrasts(
  "GAIN_vs_DIS" = GAIN - DIS,      
  "GAIN_vs_Control" = GAIN - Control,
  "DIS_vs_Control" = DIS - Control, 
  levels = design_chr17q_met
)

cat("Contrast matrix:\n")
print(contrast_matrix_chr17q_met)

# ==============================================================================
# SECTION 9: CONTRAST APPLICATION AND EMPIRICAL BAYES
# ==============================================================================

# Apply contrast matrix
fit_chr17q_met <- contrasts.fit(fit_chr17q_met, contrast_matrix_chr17q_met)

# Empirical Bayes moderation
# Methodology: Borrows information across probes to stabilize variance estimates
fit_chr17q_met <- eBayes(fit_chr17q_met)

# ==============================================================================
# SECTION 10: DIFFERENTIAL METHYLATION RESULTS EXTRACTION
# ==============================================================================

# Initialize results list
dmp_results_chr17q_list <- list()

# Extract DMP results for each contrast
for (contrast_name in colnames(contrast_matrix_chr17q_met)) {
  # Extract with FDR correction
  dmp_result <- topTable(
    fit_chr17q_met, 
    coef = contrast_name, 
    adjust = "fdr", 
    number = Inf
  )
  
  # Add contrast and probe identifiers
  dmp_result$contrast <- contrast_name
  dmp_result$probe_id <- rownames(dmp_result)
  
  dmp_results_chr17q_list[[contrast_name]] <- dmp_result
}

# Combine all results
dmp_results_chr17q_combined <- do.call(rbind, dmp_results_chr17q_list)

cat("Total DMP results:", nrow(dmp_results_chr17q_combined), "\n")
cat("BRCA1 present:", "BRCA1" %in% dmp_results_chr17q_combined$gene_name, "\n")

# ==============================================================================
# SECTION 11: PROBE ANNOTATION ENRICHMENT
# ==============================================================================

# Merge with complete probe annotation
dmp_results_chr17q_combined <- merge(
  dmp_results_chr17q_combined, 
  annot_one2one_2[, c("probeID", "gene_name", "chromosome_name")],
  by.x = "probe_id", 
  by.y = "probeID", 
  all.x = TRUE
)

# Quality control checks
cat("Probes without gene name:", 
    sum(is.na(dmp_results_chr17q_combined$gene_name) | 
        dmp_results_chr17q_combined$gene_name == ""), "\n")
cat("Probes without chromosome:", 
    sum(is.na(dmp_results_chr17q_combined$chromosome_name) | 
        dmp_results_chr17q_combined$chromosome_name == ""), "\n")

# ==============================================================================
# SECTION 12: UNIQUENESS VERIFICATION
# ==============================================================================

# Check for genes with multiple probes (should be absent after filtering)
multiple_probes_check <- dmp_results_chr17q_combined %>%
  distinct(gene_name, probe_id) %>%
  dplyr::count(gene_name) %>%
  filter(n > 1)

cat("Genes with multiple probes:", nrow(multiple_probes_check), "\n")

# Check for probes with multiple gene annotations
multiple_gene_check <- table(
  multiple = grepl(";|,|\\|", dmp_results_chr17q_combined$gene_name)
)
cat("Probe distribution by gene annotation:\n")
print(multiple_gene_check)

# Create unique identifier
dmp_results_chr17q_combined$unique_id <- paste(
  dmp_results_chr17q_combined$gene_name,
  dmp_results_chr17q_combined$contrast,
  sep = "_"
)

# Verify uniqueness
unique_check <- length(unique(dmp_results_chr17q_combined$unique_id)) == 
                nrow(dmp_results_chr17q_combined)
cat("Unique ID verification:", ifelse(unique_check, "PASSED", "FAILED"), "\n")

rownames(dmp_results_chr17q_combined) <- dmp_results_chr17q_combined$unique_id

# ==============================================================================
# SECTION 13: SIGNIFICANT PROBE FILTERING
# ==============================================================================

# Define significance criteria: |logFC| > 0.2 and FDR < 0.05
# Methodology rationale: Methylation changes are typically smaller in magnitude
# than gene expression changes. A threshold of 0.2 on M-value scale corresponds
# to approximately 15% difference in beta-value, representing biologically
# meaningful methylation change. Lower threshold than expression data reflects
# the more subtle nature of epigenetic regulation.
significant_chr17q_met_combined <- dmp_results_chr17q_combined[
  abs(dmp_results_chr17q_combined$logFC) > 0.2 & 
  dmp_results_chr17q_combined$adj.P.Val < 0.05, 
]

# Create methylation status labels
significant_chr17q_met_combined$status <- ifelse(
  significant_chr17q_met_combined$logFC > 0.2, "HYPERMETH", 
  ifelse(significant_chr17q_met_combined$logFC < -0.2, "HYPOMETH", "UNCHANGED")
)

cat("Total significant probes:", nrow(significant_chr17q_met_combined), "\n")
cat("BRCA1 in significant probes:", 
    "BRCA1" %in% significant_chr17q_met_combined$gene_name, "\n")

# Extract Chr17-specific probes
chr17_dmp_all_contrasts <- significant_chr17q_met_combined[
  significant_chr17q_met_combined$chromosome_name == "17", 
]
cat("Chr17 significant probes (all contrasts):", 
    nrow(chr17_dmp_all_contrasts), "\n")

# ==============================================================================
# SECTION 14: PRIMARY CONTRAST - GAIN VS DIS
# ==============================================================================

# Extract primary contrast
gain_vs_dis_dmp <- significant_chr17q_met_combined[
  significant_chr17q_met_combined$contrast == "GAIN_vs_DIS", 
]

cat("Significant probes in GAIN vs DIS:", nrow(gain_vs_dis_dmp), "\n")

# Chr17-specific probes
chr17_gain_vs_dis_dmp <- gain_vs_dis_dmp[
  gain_vs_dis_dmp$chromosome_name == "17", 
]

cat("Chr17 significant probes in GAIN vs DIS:", 
    nrow(chr17_gain_vs_dis_dmp), "\n")

# ==============================================================================
# SECTION 15: SECONDARY CONTRAST - GAIN VS CONTROL
# ==============================================================================

gain_vs_ctrl_dmp <- significant_chr17q_met_combined[
  significant_chr17q_met_combined$contrast == "GAIN_vs_Control", 
]

cat("Significant probes in GAIN vs Control:", nrow(gain_vs_ctrl_dmp), "\n")

chr17_gain_vs_ctrl_dmp <- gain_vs_ctrl_dmp[
  gain_vs_ctrl_dmp$chromosome_name == "17", 
]

cat("Chr17 significant probes in GAIN vs Control:", 
    nrow(chr17_gain_vs_ctrl_dmp), "\n")

# ==============================================================================
# SECTION 16: TERTIARY CONTRAST - DIS VS CONTROL
# ==============================================================================

dis_vs_ctrl_dmp <- significant_chr17q_met_combined[
  significant_chr17q_met_combined$contrast == "DIS_vs_Control", 
]

cat("Significant probes in DIS vs Control:", nrow(dis_vs_ctrl_dmp), "\n")

chr17_dis_vs_ctrl_dmp <- dis_vs_ctrl_dmp[
  dis_vs_ctrl_dmp$chromosome_name == "17", 
]

cat("Chr17 significant probes in DIS vs Control:", 
    nrow(chr17_dis_vs_ctrl_dmp), "\n")

# ==============================================================================
# SECTION 17: VISUALIZATION - VOLCANO PLOT
# ==============================================================================

# Prepare data for volcano plot (GAIN vs Control)
dmp_gain_vs_control <- dmp_results_chr17q_combined[
  dmp_results_chr17q_combined$contrast == "GAIN_vs_Control", 
]

# Filter for Chr17 probes
dmp_chr17q_gain_vs_ctrl <- dmp_results_chr17q_combined[
  dmp_results_chr17q_combined$contrast == "GAIN_vs_Control" &
  dmp_results_chr17q_combined$chromosome_name == "17",
]

# Create volcano plot
# Visualization rationale: Lower FCcutoff (0.2) reflects smaller magnitude
# of methylation changes compared to expression changes
volcano_plot_met <- EnhancedVolcano(
  toptable = dmp_chr17q_gain_vs_ctrl,            
  lab = dmp_chr17q_gain_vs_ctrl$gene_name,    
  x = 'logFC',                    
  y = 'adj.P.Val',              
  title = 'DMP Analysis: GAIN-Chr17q vs Control',
  xlab = bquote(~Log[2] ~ "Methylation Change"),
  ylab = bquote(~-Log[10] ~ "Adjusted P-value"),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  pointSize = 2.0,          
  labSize = 5.0,
  colAlpha = 0.7,
  legendPosition = 'none'
)

print(volcano_plot_met)

# Save plot
# pdf("results/volcano_plot_met_gain_vs_control.pdf", width = 10, height = 8)
# print(volcano_plot_met)
# dev.off()

# ==============================================================================
# SECTION 18: GENOMIC REGION ANALYSIS (OPTIONAL)
# ==============================================================================

# Analyze probe distribution by genomic region (if annotation available)
if("UCSC_RefGene_Group" %in% names(significant_chr17q_met_combined) &&
   !all(is.na(significant_chr17q_met_combined$UCSC_RefGene_Group))) {
  
  region_distribution <- significant_chr17q_met_combined %>%
    filter(contrast == "GAIN_vs_DIS" & !is.na(UCSC_RefGene_Group)) %>%
    dplyr::count(UCSC_RefGene_Group, status) %>%
    arrange(desc(n))
  
  cat("\nProbe distribution by genomic region (GAIN vs DIS):\n")
  print(region_distribution)
}

# ==============================================================================
# SECTION 19: RESULTS SUMMARY
# ==============================================================================

cat("\n=== DIFFERENTIAL METHYLATION ANALYSIS SUMMARY ===\n")
cat("Total probes analyzed:", nrow(dmp_results_chr17q_combined)/3, "\n")
cat("Significant probes (all contrasts):", nrow(significant_chr17q_met_combined), "\n")
cat("Significant probes - GAIN vs DIS:", 
    sum(significant_chr17q_met_combined$contrast == "GAIN_vs_DIS"), "\n")
cat("Significant probes - GAIN vs Control:", 
    sum(significant_chr17q_met_combined$contrast == "GAIN_vs_Control"), "\n")
cat("Significant probes - DIS vs Control:", 
    sum(significant_chr17q_met_combined$contrast == "DIS_vs_Control"), "\n")
cat("Chr17 significant probes (GAIN vs DIS):", nrow(chr17_gain_vs_dis_dmp), "\n")

# Analyze hyper/hypomethylation ratio for primary contrast
if(nrow(gain_vs_dis_dmp) > 0) {
  hyper_count <- sum(gain_vs_dis_dmp$status == "HYPERMETH")
  hypo_count <- sum(gain_vs_dis_dmp$status == "HYPOMETH")
  
  cat("\nHypermethylated probes (GAIN vs DIS):", hyper_count, "\n")
  cat("Hypomethylated probes (GAIN vs DIS):", hypo_count, "\n")
  cat("Hyper:Hypo ratio:", round(hyper_count/hypo_count, 2), "\n")
}

# ==============================================================================
# SECTION 20: DATA EXPORT
# ==============================================================================

# Export results
# write.csv(dmp_results_chr17q_combined, 
#           "results/dmp_chr17q_all_results.csv", 
#           row.names = FALSE)
# 
# write.csv(significant_chr17q_met_combined, 
#           "results/dmp_chr17q_significant.csv", 
#           row.names = FALSE)
# 
# write.csv(gain_vs_dis_dmp, 
#           "results/dmp_gain_vs_dis.csv", 
#           row.names = FALSE)

cat("\n=== ANALYSIS COMPLETE ===\n")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
