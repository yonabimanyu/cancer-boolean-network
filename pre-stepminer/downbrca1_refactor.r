# ===============================================================================
# PROCESS downBRCA1 SAMPLES: TRANSCRIPTOMICS AND EPIGENOMICS
# ===============================================================================
#
# Title: Data Processing Pipeline for downBRCA1 Samples
# Description: Extracts and processes chromosome 17 gene expression and DNA
#              methylation data from samples classified as downBRCA1 (downregulated
#              BRCA1). Implements probe/gene selection based on highest variance
#              to handle multiple isoforms and probes mapping to single genes.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Input Files:
#   - meta_combined: Metadata with BRCA1_status classification
#   - exp_combined: Combined gene expression matrix (genes x samples)
#   - gene_info_final_3.5: Gene annotation reference (ensembl_id, gene_symbol, chromosome)
#   - meta_combined_met: Methylation metadata with BRCA1_methylation_status
#   - met_combined: Combined DNA methylation matrix (probes x samples)
#   - annot_one2one_2: Probe annotation reference (probeID, gene_name, chromosome_name)
#
# Output Files:
#   - results/exp_downbrca1-expr.txt: Processed gene expression data (Chr17 only)
#   - results/met_downbrca1-expr.txt: Processed DNA methylation data (Chr17 only)
#
# Dependencies: data.table, dplyr
#
# ===============================================================================

# ===============================================================================
# LOAD REQUIRED LIBRARIES
# ===============================================================================

library(data.table)
library(dplyr)

# ===============================================================================
# SETUP OUTPUT DIRECTORY
# ===============================================================================

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# ===============================================================================
# PART 1: GENE EXPRESSION DATA PROCESSING (TRANSCRIPTOMICS)
# ===============================================================================

cat("=== PROCESSING downBRCA1 GENE EXPRESSION DATA ===\n")

# -------------------------------------------------------------------------------
# 1.1 Extract downBRCA1 Samples
# -------------------------------------------------------------------------------
# Scientific rationale: Sample stratification based on BRCA1 expression status
# allows investigation of BRCA1-specific regulatory mechanisms and downstream
# effects on chromosome 17 genes. downBRCA1 samples represent the complementary
# regulatory state to upBRCA1, enabling comparative analysis
# -------------------------------------------------------------------------------

cat("1. Extracting downBRCA1 samples from metadata...\n")

downbrca1_sample_ids <- rownames(meta_combined)[
  meta_combined$BRCA1_status == "downBRCA1" & 
  !is.na(meta_combined$BRCA1_status)
]

cat("   Found", length(downbrca1_sample_ids), "downBRCA1 samples\n")

# Subset expression matrix to downBRCA1 samples only
exp_matrix_downbrca1 <- exp_combined[, downbrca1_sample_ids]

# -------------------------------------------------------------------------------
# 1.2 Filter Non-Variable Genes
# -------------------------------------------------------------------------------
# Scientific rationale: Genes with zero variance across samples provide no
# biological information and can cause numerical issues in downstream analyses.
# Removing them improves computational efficiency without loss of information
# -------------------------------------------------------------------------------

gene_variance_downbrca1 <- apply(exp_matrix_downbrca1, 1, var)
exp_matrix_downbrca1 <- exp_matrix_downbrca1[gene_variance_downbrca1 > 0, ]

cat("   Genes with variance > 0:", sum(gene_variance_downbrca1 > 0), "\n")
cat("   Matrix dimensions:", dim(exp_matrix_downbrca1)[1], "genes x", 
    dim(exp_matrix_downbrca1)[2], "samples\n")

# -------------------------------------------------------------------------------
# 1.3 Map Ensembl IDs to Gene Symbols
# -------------------------------------------------------------------------------
# Scientific rationale: Ensembl IDs with version numbers (e.g., ENSG.23) need
# standardization for proper mapping. Gene symbols provide human-readable
# identifiers for biological interpretation
# -------------------------------------------------------------------------------

cat("2. Mapping Ensembl IDs to gene symbols...\n")

# Remove version suffixes from Ensembl IDs (ENSG00000012048.23 -> ENSG00000012048)
ensembl_ids_clean <- gsub("\\..*", "", rownames(exp_matrix_downbrca1))

# Create mapping dataframe with row indices for later subsetting
ensembl_mapping_df <- data.frame(
  ensembl_id = ensembl_ids_clean,
  row_index = seq_len(nrow(exp_matrix_downbrca1)),
  stringsAsFactors = FALSE
)

# Merge with gene annotation reference
gene_annotation_downbrca1 <- merge(
  ensembl_mapping_df,
  gene_info_final_3.5[, c("ensembl_id", "gene_symbol", "chromosome")],
  by = "ensembl_id",
  all.x = TRUE
)

# -------------------------------------------------------------------------------
# 1.4 Quality Control
# -------------------------------------------------------------------------------

cat("3. Performing quality control...\n")

missing_ensembl <- sum(is.na(gene_annotation_downbrca1$ensembl_id) | 
                       gene_annotation_downbrca1$ensembl_id == "")
missing_symbol <- sum(is.na(gene_annotation_downbrca1$gene_symbol) | 
                      gene_annotation_downbrca1$gene_symbol == "")

cat("   Missing Ensembl IDs:", missing_ensembl, "\n")
cat("   Missing Gene Symbols:", missing_symbol, "\n")

unmapped_genes <- gene_annotation_downbrca1[
  !is.na(gene_annotation_downbrca1$ensembl_id) & 
  gene_annotation_downbrca1$ensembl_id != "" & 
  (is.na(gene_annotation_downbrca1$gene_symbol) | 
   gene_annotation_downbrca1$gene_symbol == ""), 
]
cat("   Ensembl IDs without gene symbols:", nrow(unmapped_genes), "\n")

# -------------------------------------------------------------------------------
# 1.5 Filter for Chromosome 17 Genes
# -------------------------------------------------------------------------------
# Scientific rationale: Analysis focuses on chromosome 17 genes to investigate
# cis-regulatory effects and local genomic alterations associated with BRCA1
# dysregulation. This reduces dimensionality while maintaining biological focus
# -------------------------------------------------------------------------------

cat("4. Filtering for chromosome 17 genes...\n")

# Remove genes without valid gene symbols
gene_annotation_downbrca1 <- gene_annotation_downbrca1[
  !(is.na(gene_annotation_downbrca1$gene_symbol) | 
    gene_annotation_downbrca1$gene_symbol == ""),
]

cat("   Genes with valid symbols:", nrow(gene_annotation_downbrca1), "\n")

# Filter for chromosome 17 only
gene_annotation_chr17_downbrca1 <- gene_annotation_downbrca1[
  !is.na(gene_annotation_downbrca1$chromosome) & 
  gene_annotation_downbrca1$chromosome == "17",
]

cat("   Chromosome 17 genes:", nrow(gene_annotation_chr17_downbrca1), "\n")

# Create data.table with expression values
expr_data_chr17_downbrca1 <- data.table(
  ensembl_id = gene_annotation_chr17_downbrca1$ensembl_id,
  gene_symbol = gene_annotation_chr17_downbrca1$gene_symbol,
  chromosome = gene_annotation_chr17_downbrca1$chromosome,
  row_index = gene_annotation_chr17_downbrca1$row_index,
  exp_matrix_downbrca1[gene_annotation_chr17_downbrca1$row_index, ]
)

# -------------------------------------------------------------------------------
# 1.6 Select Representative Gene Isoforms
# -------------------------------------------------------------------------------
# Scientific rationale: Multiple Ensembl IDs can map to the same gene symbol
# (alternative isoforms, transcripts). We select the most variable isoform
# using MAD (Median Absolute Deviation) as it is robust to outliers compared
# to variance. This ensures biological signal is preserved while avoiding
# redundancy in downstream analyses
# -------------------------------------------------------------------------------

cat("5. Selecting representative gene isoforms (highest MAD)...\n")

sample_columns_expr <- names(expr_data_chr17_downbrca1)[
  !names(expr_data_chr17_downbrca1) %in% c("ensembl_id", "gene_symbol", 
                                            "chromosome", "row_index")
]

representative_genes_downbrca1 <- expr_data_chr17_downbrca1[, {
  if (length(sample_columns_expr) > 0) {
    # Calculate MAD for each isoform
    isoform_mad_values <- apply(.SD[, sample_columns_expr, with = FALSE], 
                                1, function(x) {
      mad(as.numeric(x), na.rm = TRUE, constant = 1.4826) 
    })
    # Select isoform with highest MAD
    .SD[which.max(isoform_mad_values)]
  }
}, by = gene_symbol]

cat("   Representative genes selected:", nrow(representative_genes_downbrca1), "\n")

# -------------------------------------------------------------------------------
# 1.7 Finalize and Export Expression Data
# -------------------------------------------------------------------------------

cat("6. Finalizing and exporting expression data...\n")

# Reorder columns for readability
setcolorder(representative_genes_downbrca1, 
            c("ensembl_id", "gene_symbol", 
              setdiff(colnames(representative_genes_downbrca1), 
                      c("ensembl_id", "gene_symbol", "chromosome", "row_index"))))

# Remove unnecessary columns
representative_genes_downbrca1[, c("row_index", "chromosome") := NULL]

# Rename columns to standard format
representative_genes_downbrca1 <- representative_genes_downbrca1 %>%
  rename(
    ProbeID = ensembl_id,
    Name = gene_symbol
  )

# Export to file
write.table(representative_genes_downbrca1, 
            file = "results/exp_downbrca1-expr.txt", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

cat("   Expression data saved to: results/exp_downbrca1-expr.txt\n")
cat("=== GENE EXPRESSION PROCESSING COMPLETE ===\n\n")

# ===============================================================================
# PART 2: DNA METHYLATION DATA PROCESSING (EPIGENOMICS)
# ===============================================================================

cat("=== PROCESSING downBRCA1 DNA METHYLATION DATA ===\n")

# -------------------------------------------------------------------------------
# 2.1 Extract downBRCA1 Samples
# -------------------------------------------------------------------------------
# Scientific rationale: DNA methylation patterns are analyzed separately from
# expression as they represent distinct regulatory layers. Coordinated analysis
# of both omics layers reveals epigenetic mechanisms underlying gene regulation
# -------------------------------------------------------------------------------

cat("1. Extracting downBRCA1 samples from methylation metadata...\n")

downbrca1_met_sample_ids <- rownames(meta_combined_met)[
  meta_combined_met$BRCA1_methylation_status == "downBRCA1" & 
  !is.na(meta_combined_met$BRCA1_methylation_status)
]

cat("   Found", length(downbrca1_met_sample_ids), "downBRCA1 methylation samples\n")

# Subset methylation matrix
met_matrix_downbrca1 <- met_combined[, downbrca1_met_sample_ids]

cat("   Matrix dimensions:", dim(met_matrix_downbrca1)[1], "probes x", 
    dim(met_matrix_downbrca1)[2], "samples\n")

# -------------------------------------------------------------------------------
# 2.2 Map Probes to Genes
# -------------------------------------------------------------------------------

cat("2. Mapping methylation probes to genes...\n")

probe_ids_downbrca1 <- rownames(met_matrix_downbrca1)

probe_mapping_df <- data.frame(
  probeID = probe_ids_downbrca1,
  row_index = seq_len(nrow(met_matrix_downbrca1)),
  stringsAsFactors = FALSE
)

# Remove invalid probe IDs
initial_probe_count <- nrow(probe_mapping_df)
probe_mapping_df <- probe_mapping_df[
  !(is.na(probe_mapping_df$probeID) | probe_mapping_df$probeID == ""), 
]
cat("   Valid probe IDs:", nrow(probe_mapping_df), "of", initial_probe_count, "\n")

# -------------------------------------------------------------------------------
# 2.3 Annotate Probes with Gene Information
# -------------------------------------------------------------------------------

cat("3. Annotating probes with gene information...\n")

probe_annotation_downbrca1 <- merge(
  probe_mapping_df,
  annot_one2one_2[, c("probeID", "gene_name", "chromosome_name")],
  by = "probeID",
  all.x = TRUE
)

# Filter for valid gene annotations
probe_annotation_downbrca1 <- probe_annotation_downbrca1[
  !(is.na(probe_annotation_downbrca1$gene_name) | 
    probe_annotation_downbrca1$gene_name == ""), 
]

cat("   Probes with valid gene annotations:", nrow(probe_annotation_downbrca1), "\n")

# -------------------------------------------------------------------------------
# 2.4 Filter for Chromosome 17 Probes
# -------------------------------------------------------------------------------
# Scientific rationale: Consistent with expression analysis, methylation
# analysis focuses on chromosome 17 to enable integrative multi-omics
# investigation of local regulatory mechanisms
# -------------------------------------------------------------------------------

probe_annotation_chr17_downbrca1 <- probe_annotation_downbrca1[
  !is.na(probe_annotation_downbrca1$chromosome_name) & 
  probe_annotation_downbrca1$chromosome_name == "17",
]

cat("   Chromosome 17 probes:", nrow(probe_annotation_chr17_downbrca1), "\n")

# -------------------------------------------------------------------------------
# 2.5 Create Methylation Data Table
# -------------------------------------------------------------------------------

cat("4. Creating methylation data table...\n")

met_data_chr17_downbrca1 <- data.table(
  probeID = probe_annotation_chr17_downbrca1$probeID,
  gene_name = probe_annotation_chr17_downbrca1$gene_name,
  chromosome_name = probe_annotation_chr17_downbrca1$chromosome_name,
  row_index = probe_annotation_chr17_downbrca1$row_index,
  met_matrix_downbrca1[probe_annotation_chr17_downbrca1$row_index, ]
)

# -------------------------------------------------------------------------------
# 2.6 Select Representative Probes per Gene
# -------------------------------------------------------------------------------
# Scientific rationale: Multiple probes often target different CpG sites
# within the same gene. Selecting the probe with highest variance captures
# the most dynamic methylation changes and reduces data redundancy while
# maintaining biological signal
# -------------------------------------------------------------------------------

cat("5. Selecting representative probes (highest variance)...\n")

sample_columns_met <- names(met_data_chr17_downbrca1)[
  !names(met_data_chr17_downbrca1) %in% c("probeID", "gene_name", 
                                           "chromosome_name", "row_index")
]

representative_probes_downbrca1 <- met_data_chr17_downbrca1[, {
  if(length(sample_columns_met) > 0) {
    # Calculate variance for each probe targeting the same gene
    probe_variance_values <- apply(.SD[, sample_columns_met, with = FALSE], 
                                   1, function(x) {
      var(as.numeric(x), na.rm = TRUE)
    })
    # Select probe with highest variance
    .SD[which.max(probe_variance_values)]
  }
}, by = gene_name]

cat("   Representative probes selected:", nrow(representative_probes_downbrca1), "\n")

# -------------------------------------------------------------------------------
# 2.7 Finalize and Export Methylation Data
# -------------------------------------------------------------------------------

cat("6. Finalizing and exporting methylation data...\n")

# Reorder columns
setcolorder(representative_probes_downbrca1, 
            c("probeID", "gene_name", 
              setdiff(colnames(representative_probes_downbrca1), 
                      c("probeID", "gene_name", "chromosome_name", "row_index"))))

# Standardize column names
setnames(representative_probes_downbrca1, "gene_name", "gene_symbol")

# Remove unnecessary columns
representative_probes_downbrca1[, c("row_index", "chromosome_name") := NULL]

# Rename to standard format
representative_probes_downbrca1 <- representative_probes_downbrca1 %>%
  rename(
    ProbeID = probeID,
    Name = gene_symbol
  )

# Export to file
write.table(representative_probes_downbrca1, 
            file = "results/met_downbrca1-expr.txt", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

cat("   Methylation data saved to: results/met_downbrca1-expr.txt\n")
cat("=== DNA METHYLATION PROCESSING COMPLETE ===\n\n")

# ===============================================================================
# END OF SCRIPT
# ===============================================================================
