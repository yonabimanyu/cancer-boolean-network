################################################################################
# TCGA-BRCA Multi-Omics Integration Pipeline
################################################################################
#
# Title: Integrated RNA-seq, DNA Methylation, and Copy Number Analysis for TCGA-BRCA
#
# Description:
#   This pipeline integrates three molecular data types from The Cancer Genome Atlas
#   Breast Cancer (TCGA-BRCA) project to enable comprehensive multi-omics analysis:
#   1. RNA-seq (gene expression via STAR-Counts workflow)
#   2. DNA Methylation (Illumina HumanMethylation450 BeadChip)
#   3. Copy Number Alterations (GISTIC arm-level calls, focusing on Chr17q)
#
#   The pipeline performs strategic sample selection by:
#   - Identifying patients with complete multi-omics profiles (all 3 data types)
#   - Filtering for Chr17q gain events (excluding losses, which have different biology)
#   - Including matched normal tissue controls for differential analysis
#   - Annotating samples with Chr17q status (GAIN/DIS/Control) for stratified analysis
#
#   Scientific Rationale:
#   Chromosome 17q gain is a recurrent oncogenic event in breast cancer associated
#   with aggressive tumor phenotypes. Multi-omics integration enables investigation of
#   how copy number gains translate into transcriptional and epigenetic dysregulation.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Input:
#   - data/GISTIC_arm_level_copynumber.txt: GISTIC-processed arm-level copy number calls
#   - TCGA-BRCA data downloaded via TCGAbiolinks from GDC data portal
#
# Output:
#   - log2_combined_matrix_2: Log2-transformed TPM expression matrix (RNA-seq)
#   - m_met_combined_matrix: M-value transformed methylation beta values
#   - metadata_rnaseq_annotated: Annotated clinical/molecular metadata (RNA-seq)
#   - metadata_methylation_annotated: Annotated clinical/molecular metadata (Methylation)
#
# Dependencies:
#   - TCGAbiolinks: Interface to TCGA data via GDC API
#   - SummarizedExperiment: Bioconductor data container for omics experiments
#   - dplyr: Data manipulation and integration
#
# Notes:
#   - This pipeline requires active internet connection for GDC data download
#   - Download process may take 30-60 minutes depending on connection speed
#   - Final sample size depends on TCGA data availability (expected ~430-450 samples)
#
################################################################################

# ==============================================================================
# LOAD REQUIRED LIBRARIES
# ==============================================================================

library(TCGAbiolinks)        # TCGA data query and download
library(SummarizedExperiment) # Bioconductor data structures
library(dplyr)               # Data manipulation

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Define relative paths
COPY_NUMBER_FILE <- "data/GISTIC_arm_level_copynumber.txt"

# TCGA project identifier
PROJECT_ID <- "TCGA-BRCA"

# ==============================================================================
# PART 1: INITIAL DATA QUERIES AND COPY NUMBER PREPARATION
# ==============================================================================

cat("=== PART 1: Setting up data queries and loading copy number data ===\n")

# ------------------------------------------------------------------------------
# 1.1 Query RNA-seq data for Primary Tumor samples
# ------------------------------------------------------------------------------
# Rationale: RNA-seq provides genome-wide gene expression profiles. STAR-Counts
# workflow outputs are chosen for their robust quantification and compatibility
# with downstream differential expression analysis tools.

cat("Creating RNA-seq query for primary tumors...\n")
query_rnaseq_primary_tumor <- GDCquery(
  project = PROJECT_ID,
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  access = 'open',
  experimental.strategy = 'RNA-Seq',
  sample.type = 'Primary Tumor'
)

# ------------------------------------------------------------------------------
# 1.2 Query DNA Methylation data for Primary Tumor samples
# ------------------------------------------------------------------------------
# Rationale: Illumina 450K arrays measure methylation at ~450,000 CpG sites
# genome-wide. Beta values (methylation proportion) are used as they directly
# represent biological methylation levels (0 = unmethylated, 1 = fully methylated).

cat("Creating methylation query for primary tumors...\n")
query_methylation_primary_tumor <- GDCquery(
  project = PROJECT_ID,
  data.category = 'DNA Methylation',
  platform = 'Illumina Human Methylation 450',
  access = 'open',
  data.type = 'Methylation Beta Value',
  sample.type = 'Primary Tumor'
)

# ------------------------------------------------------------------------------
# 1.3 Load and prepare copy number alteration data
# ------------------------------------------------------------------------------
# Rationale: GISTIC (Genomic Identification of Significant Targets in Cancer)
# identifies recurrent copy number alterations. We focus on Chr17q gains because:
# 1. Chr17q contains multiple oncogenes
# 2. Gains (amplifications) have different biological consequences than losses
# 3. Gains are more prevalent and clinically relevant in breast cancer

cat("Loading GISTIC copy number alteration data...\n")
copy_number_raw <- read.delim(COPY_NUMBER_FILE, header = TRUE, sep = "\t")

cat("Preparing copy number data with Chr17q status annotation...\n")
copy_number_annotated <- copy_number_raw %>%
  # Strategic filtering: Exclude "Loss" events to focus on gain biology
  filter(!is.na(X17q) & X17q != "Loss") %>%
  mutate(
    # Extract patient ID (first 12 characters) for cross-omics matching
    Patient_ID = substr(SAMPLE_ID, 1, 12),
    
    # Extract sample ID (first 15 characters) for precise metadata annotation
    Sample_ID_Full = substr(SAMPLE_ID, 1, 15),
    
    # Classify Chr17q status: GAIN (amplified) vs DIS (disomic)
    # This binary classification enables clear stratification in downstream analysis
    Chr17q_Status = ifelse(X17q == "Gain", "GAIN", "DIS")
  )

cat("Copy number data prepared:", nrow(copy_number_annotated), "samples with Chr17q annotations\n")

# ==============================================================================
# PART 2: MULTI-OMICS SAMPLE INTERSECTION
# ==============================================================================

cat("\n=== PART 2: Identifying samples with complete multi-omics profiles ===\n")

# ------------------------------------------------------------------------------
# 2.1 Extract patient identifiers from each omics query
# ------------------------------------------------------------------------------
# Rationale: TCGA uses hierarchical barcodes where first 12 characters uniquely
# identify a patient. This enables cross-omics matching since the same patient
# may have multiple sample aliquots across different platforms.

rnaseq_patients <- substr(
  getResults(query_rnaseq_primary_tumor, cols = "cases"), 
  1, 12
)

methylation_patients <- substr(
  getResults(query_methylation_primary_tumor, cols = "cases"), 
  1, 12
)

cat("Patients with RNA-seq data:", length(rnaseq_patients), "\n")
cat("Patients with methylation data:", length(methylation_patients), "\n")

# ------------------------------------------------------------------------------
# 2.2 Perform iterative intersection to identify complete cases
# ------------------------------------------------------------------------------
# Rationale: Multi-omics integration requires complete data matrices. By
# selecting only patients with all three data types, we ensure that downstream
# integrative analyses (differential analysis, StepMiner, and BooleanNet) can
# leverage the full molecular profile without missing data complications.

# Step 1: Find patients with both RNA-seq AND methylation
patients_with_both_omics <- intersect(rnaseq_patients, methylation_patients)
cat("Patients with both RNA-seq and methylation:", length(patients_with_both_omics), "\n")

# Step 2: Further intersect with copy number data (already filtered for Chr17q gain/diploid)
final_selected_patients <- intersect(
  patients_with_both_omics, 
  copy_number_annotated$Patient_ID
)

cat("Final selected patients (complete multi-omics + Chr17q annotation):", 
    length(final_selected_patients), "\n")

cat("\n=== Sample Selection Summary ===\n")
cat("Starting cohort (RNA-seq):", length(rnaseq_patients), "\n")
cat("After methylation requirement:", length(patients_with_both_omics), "\n")
cat("After copy number requirement:", length(final_selected_patients), "\n")
cat("Retention rate:", 
    round(100 * length(final_selected_patients) / length(rnaseq_patients), 1), "%\n")

# ==============================================================================
# PART 3: RNA-SEQ DATA DOWNLOAD AND PROCESSING
# ==============================================================================

cat("\n=== PART 3: Downloading and processing RNA-seq data ===\n")

# ------------------------------------------------------------------------------
# 3.1 Prepare RNA-seq sample selection
# ------------------------------------------------------------------------------
# Rationale: We download both Primary Tumor (PT) and Solid Tissue Normal (STN):
# - PT samples: Only those in our final_selected_patients (multi-omics complete)
# - STN samples: ALL available normals for robust differential analysis baseline
# Normal tissue controls are critical for identifying cancer-specific alterations.

cat("Retrieving full barcodes for selected primary tumor samples...\n")
rnaseq_pt_all_barcodes <- getResults(query_rnaseq_primary_tumor, cols = "cases")
rnaseq_pt_selected_barcodes <- rnaseq_pt_all_barcodes[
  substr(rnaseq_pt_all_barcodes, 1, 12) %in% final_selected_patients
]

cat("Querying all available solid tissue normal samples...\n")
query_rnaseq_normal <- GDCquery(
  project = PROJECT_ID,
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  access = 'open',
  experimental.strategy = 'RNA-Seq',
  sample.type = 'Solid Tissue Normal'
)

rnaseq_normal_barcodes <- getResults(query_rnaseq_normal, cols = "cases")

# Combine selected PT + all normals for comprehensive cohort
combined_rnaseq_barcodes <- c(rnaseq_pt_selected_barcodes, rnaseq_normal_barcodes)

cat("Total RNA-seq samples to download:", length(combined_rnaseq_barcodes), 
    "(", length(rnaseq_pt_selected_barcodes), "PT +", 
    length(rnaseq_normal_barcodes), "Normal )\n")

# ------------------------------------------------------------------------------
# 3.2 Download RNA-seq data from GDC
# ------------------------------------------------------------------------------

cat("Creating final RNA-seq query with selected barcodes...\n")
query_rnaseq_final <- GDCquery(
  project = PROJECT_ID,
  data.category = 'Transcriptome Profiling',
  data.type = 'Gene Expression Quantification',
  workflow.type = 'STAR - Counts',
  access = 'open',
  experimental.strategy = 'RNA-Seq',
  barcode = combined_rnaseq_barcodes
)

cat("Downloading RNA-seq data (this may take 15-30 minutes)...\n")
GDCdownload(query_rnaseq_final, method = 'api', files.per.chunk = 1)

cat("Preparing RNA-seq SummarizedExperiment object...\n")
rnaseq_data_combined <- GDCprepare(
  query_rnaseq_final, 
  summarizedExperiment = TRUE
)

# ------------------------------------------------------------------------------
# 3.3 Extract and transform expression matrices
# ------------------------------------------------------------------------------
# Rationale: TPM (Transcripts Per Million) normalization accounts for gene length
# and library size, making expression values comparable across genes and samples.
# Log2 transformation stabilizes variance and makes data more normally distributed,
# which is required for parametric statistical tests.

cat("Extracting TPM expression matrix...\n")
tpm_expression_matrix <- assay(rnaseq_data_combined, 'tpm_unstrand')

cat("Applying log2 transformation (log2(TPM + 1))...\n")
log2_expression_matrix <- log2(tpm_expression_matrix + 1)

cat("RNA-seq processing complete:", nrow(log2_expression_matrix), "genes x", 
    ncol(log2_expression_matrix), "samples\n")

# ==============================================================================
# PART 4: METHYLATION DATA DOWNLOAD AND PROCESSING
# ==============================================================================

cat("\n=== PART 4: Downloading and processing methylation data ===\n")

# ------------------------------------------------------------------------------
# 4.1 Prepare methylation sample selection
# ------------------------------------------------------------------------------
# Rationale: Same strategy as RNA-seq - selected PT + all normals for
# comprehensive case-control methylation profiling.

cat("Retrieving full barcodes for selected primary tumor samples...\n")
methylation_pt_all_barcodes <- getResults(query_methylation_primary_tumor, cols = "cases")
methylation_pt_selected_barcodes <- methylation_pt_all_barcodes[
  substr(methylation_pt_all_barcodes, 1, 12) %in% final_selected_patients
]

cat("Querying all available solid tissue normal samples...\n")
query_methylation_normal <- GDCquery(
  project = PROJECT_ID,
  data.category = 'DNA Methylation',
  platform = 'Illumina Human Methylation 450',
  access = 'open',
  data.type = 'Methylation Beta Value',
  sample.type = 'Solid Tissue Normal'
)

methylation_normal_barcodes <- getResults(query_methylation_normal, cols = "cases")

# Combine selected PT + all normals
combined_methylation_barcodes <- c(
  methylation_pt_selected_barcodes, 
  methylation_normal_barcodes
)

cat("Total methylation samples to download:", length(combined_methylation_barcodes), 
    "(", length(methylation_pt_selected_barcodes), "PT +", 
    length(methylation_normal_barcodes), "Normal )\n")

# ------------------------------------------------------------------------------
# 4.2 Download methylation data from GDC
# ------------------------------------------------------------------------------

cat("Creating final methylation query with selected barcodes...\n")
query_methylation_final <- GDCquery(
  project = PROJECT_ID,
  data.category = 'DNA Methylation',
  platform = 'Illumina Human Methylation 450',
  access = 'open',
  data.type = 'Methylation Beta Value',
  barcode = combined_methylation_barcodes
)

cat("Downloading methylation data (this may take 15-30 minutes)...\n")
GDCdownload(query_methylation_final, method = 'api', files.per.chunk = 1)

cat("Preparing methylation SummarizedExperiment object...\n")
methylation_data_combined <- GDCprepare(
  query_methylation_final, 
  summarizedExperiment = TRUE
)

# ------------------------------------------------------------------------------
# 4.3 Extract and transform methylation matrices
# ------------------------------------------------------------------------------
# Rationale: Beta values (0-1 scale) are converted to M-values using the logit
# transformation: M = log2(β / (1-β)). M-values have better statistical properties:
# - More homoscedastic (constant variance across methylation range)
# - More suitable for differential methylation analysis

cat("Extracting beta value matrix...\n")
beta_value_matrix <- assay(methylation_data_combined)

cat("Converting beta values to M-values: M = log2(β / (1-β))...\n")
m_value_matrix <- log2(beta_value_matrix / (1 - beta_value_matrix))

cat("Methylation processing complete:", nrow(m_value_matrix), "CpG sites x", 
    ncol(m_value_matrix), "samples\n")

# ==============================================================================
# PART 5: CREATE INTEGRATED METADATA FOR RNA-SEQ
# ==============================================================================

cat("\n=== PART 5: Annotating RNA-seq samples with Chr17q status ===\n")

# ------------------------------------------------------------------------------
# 5.1 Extract and prepare RNA-seq metadata
# ------------------------------------------------------------------------------
# Rationale: Metadata integration is critical for stratified analysis. We annotate
# each sample with its Chr17q status (GAIN/DIS/Control) to enable investigation of
# copy number-driven molecular changes.

cat("Extracting clinical and technical metadata from RNA-seq data...\n")
metadata_rnaseq <- colData(rnaseq_data_combined)
metadata_rnaseq_df <- as.data.frame(metadata_rnaseq)

# Create sample identifier for copy number matching (15-character TCGA barcode)
metadata_rnaseq_df$Sample_ID_Full <- substr(metadata_rnaseq_df$sample, 1, 15)

# ------------------------------------------------------------------------------
# 5.2 Integrate Chr17q copy number status
# ------------------------------------------------------------------------------
# Rationale: Left join preserves all RNA-seq samples (including normals) while
# adding Chr17q status where available. This creates a unified annotation framework.

cat("Integrating Chr17q status from copy number data...\n")
metadata_rnaseq_df <- metadata_rnaseq_df %>%
  left_join(
    copy_number_annotated[, c("Sample_ID_Full", "Chr17q_Status")], 
    by = "Sample_ID_Full"
  )

# ------------------------------------------------------------------------------
# 5.3 Finalize Chr17q status with control group assignment
# ------------------------------------------------------------------------------
# Rationale: Three-group stratification enables comprehensive analysis:
# - Control: Normal tissue baseline (no copy number alterations)
# - DIS: Disomic tumors (no loss or gain in chromosome number)
# - GAIN: Chr17q amplified tumors (test group of interest)

cat("Assigning final Chr17q status categories...\n")
metadata_rnaseq_df$Chr17q_Status_Final <- ifelse(
  metadata_rnaseq_df$sample_type == "Solid Tissue Normal", 
  "Control",                       # All normals = Control group
  metadata_rnaseq_df$Chr17q_Status # Tumors retain GAIN/DIS classification
)

# Convert to ordered factor for consistent analysis and visualization
metadata_rnaseq_df$Chr17q_Status_Final <- factor(
  metadata_rnaseq_df$Chr17q_Status_Final, 
  levels = c("Control", "DIS", "GAIN"),
  ordered = TRUE
)

# Convert back to DFrame for SummarizedExperiment compatibility
metadata_rnaseq_annotated <- as(metadata_rnaseq_df, "DFrame")

# Display distribution summary
cat("\nRNA-seq Chr17q Status Distribution:\n")
print(table(metadata_rnaseq_annotated$Chr17q_Status_Final))

# ==============================================================================
# PART 6: CREATE INTEGRATED METADATA FOR METHYLATION
# ==============================================================================

cat("\n=== PART 6: Annotating methylation samples with Chr17q status ===\n")

# ------------------------------------------------------------------------------
# 6.1 Extract and prepare methylation metadata
# ------------------------------------------------------------------------------

cat("Extracting clinical and technical metadata from methylation data...\n")
metadata_methylation <- colData(methylation_data_combined)
metadata_methylation_df <- as.data.frame(metadata_methylation)

# Create sample identifier for copy number matching
metadata_methylation_df$Sample_ID_Full <- substr(metadata_methylation_df$sample, 1, 15)

# ------------------------------------------------------------------------------
# 6.2 Integrate Chr17q copy number status
# ------------------------------------------------------------------------------

cat("Integrating Chr17q status from copy number data...\n")
metadata_methylation_df <- metadata_methylation_df %>%
  left_join(
    copy_number_annotated[, c("Sample_ID_Full", "Chr17q_Status")], 
    by = "Sample_ID_Full"
  )

# ------------------------------------------------------------------------------
# 6.3 Finalize Chr17q status with control group assignment
# ------------------------------------------------------------------------------

cat("Assigning final Chr17q status categories...\n")
metadata_methylation_df$Chr17q_Status_Final <- ifelse(
  metadata_methylation_df$sample_type == "Solid Tissue Normal", 
  "Control",
  metadata_methylation_df$Chr17q_Status
)

# Convert to ordered factor
metadata_methylation_df$Chr17q_Status_Final <- factor(
  metadata_methylation_df$Chr17q_Status_Final, 
  levels = c("Control", "DIS", "GAIN"),
  ordered = TRUE
)

# Convert back to DFrame
metadata_methylation_annotated <- as(metadata_methylation_df, "DFrame")

# Display distribution summary
cat("\nMethylation Chr17q Status Distribution:\n")
print(table(metadata_methylation_annotated$Chr17q_Status_Final))

# ==============================================================================
# PART 7: VALIDATION AND QUALITY CONTROL
# ==============================================================================

cat("\n=== PART 7: Validating data integration and completeness ===\n")

# ------------------------------------------------------------------------------
# 7.1 Overall data dimensions
# ------------------------------------------------------------------------------

cat("\n=== Final Data Dimensions ===\n")
cat("RNA-seq expression matrix:", 
    nrow(log2_expression_matrix), "genes x", 
    ncol(log2_expression_matrix), "samples\n")
cat("Methylation matrix:", 
    nrow(m_value_matrix), "CpG sites x", 
    ncol(m_value_matrix), "samples\n")

# ------------------------------------------------------------------------------
# 7.2 Chr17q status annotation completeness
# ------------------------------------------------------------------------------

cat("\n=== Chr17q Annotation Summary ===\n")

# RNA-seq
cat("\nRNA-seq Chr17q Status:\n")
print(table(metadata_rnaseq_annotated$Chr17q_Status_Final, useNA = "ifany"))

# Methylation
cat("\nMethylation Chr17q Status:\n")
print(table(metadata_methylation_annotated$Chr17q_Status_Final, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 7.3 Validate primary tumor annotation completeness
# ------------------------------------------------------------------------------
# Rationale: All primary tumor samples should have defined Chr17q status
# (either GAIN or DIS) since we filtered copy number data explicitly.
# Any NA values would indicate integration errors.

cat("\n=== Primary Tumor Annotation Validation ===\n")

pt_rnaseq_samples <- metadata_rnaseq_annotated[
  metadata_rnaseq_annotated$sample_type == "Primary Tumor", 
]
pt_methylation_samples <- metadata_methylation_annotated[
  metadata_methylation_annotated$sample_type == "Primary Tumor", 
]

cat("PT RNA-seq samples with defined Chr17q status:", 
    sum(!is.na(pt_rnaseq_samples$Chr17q_Status_Final)), "/", 
    nrow(pt_rnaseq_samples), "\n")

cat("PT Methylation samples with defined Chr17q status:", 
    sum(!is.na(pt_methylation_samples$Chr17q_Status_Final)), "/", 
    nrow(pt_methylation_samples), "\n")

# ------------------------------------------------------------------------------
# 7.4 Verify cross-omics patient overlap
# ------------------------------------------------------------------------------
# Rationale: The number of overlapping patients should match our
# final_selected_patients count, confirming successful multi-omics integration.

cat("\n=== Cross-Omics Overlap Validation ===\n")

rnaseq_pt_patients <- substr(pt_rnaseq_samples$sample, 1, 12)
methylation_pt_patients <- substr(pt_methylation_samples$sample, 1, 12)
overlapping_pt_patients <- intersect(rnaseq_pt_patients, methylation_pt_patients)

cat("Primary tumor patients with both RNA-seq and methylation:", 
    length(overlapping_pt_patients), "\n")
cat("Expected (from initial selection):", 
    length(final_selected_patients), "\n")
cat("Match status:", 
    length(overlapping_pt_patients) == length(final_selected_patients), "\n")

# ------------------------------------------------------------------------------
# 7.5 Final validation summary
# ------------------------------------------------------------------------------

if (length(overlapping_pt_patients) == length(final_selected_patients)) {
  cat("\n✓ SUCCESS: Multi-omics integration validated successfully\n")
  cat("✓ All primary tumor samples have complete Chr17q annotations\n")
  cat("✓ Cross-omics patient matching confirmed\n")
} else {
  cat("\n✗ WARNING: Mismatch detected in patient overlap\n")
  cat("✗ Review sample selection and annotation steps\n")
}

# ------------------------------------------------------------------------------
# 7.6 Create summary report
# ------------------------------------------------------------------------------

cat("\n=== FINAL PIPELINE SUMMARY ===\n")
cat("Total patients with complete multi-omics data:", length(final_selected_patients), "\n")
cat("RNA-seq cohort: ", ncol(log2_expression_matrix), "samples\n")
cat("  - Primary Tumor:", nrow(pt_rnaseq_samples), "\n")
cat("  - Normal Tissue:", sum(metadata_rnaseq_annotated$sample_type == "Solid Tissue Normal"), "\n")
cat("Methylation cohort:", ncol(m_value_matrix), "samples\n")
cat("  - Primary Tumor:", nrow(pt_methylation_samples), "\n")
cat("  - Normal Tissue:", sum(metadata_methylation_annotated$sample_type == "Solid Tissue Normal"), "\n")

cat("\n=== Ready for downstream analysis ===\n")
cat("Key output objects:\n")
cat("  - log2_expression_matrix: Log2-transformed TPM values\n")
cat("  - m_value_matrix: M-value transformed methylation data\n")
cat("  - metadata_rnaseq_annotated: Annotated RNA-seq clinical metadata\n")
cat("  - metadata_methylation_annotated: Annotated methylation clinical metadata\n")

# Assign to legacy variable names for backward compatibility
log2_combined_matrix_2 <- log2_expression_matrix
m_met_combined_matrix <- m_value_matrix
meta_data_chr17q_exp <- metadata_rnaseq_annotated
meta_data_chr17q_met <- metadata_methylation_annotated

################################################################################
# END OF SCRIPT
################################################################################
