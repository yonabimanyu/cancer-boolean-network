################################################################################
# Gene and Probe Annotation Pipeline for Illumina 450K Array
################################################################################
#
# Title: Comprehensive Gene Annotation and Coordinate Liftover Pipeline
#
# Description: 
#   This script creates a unified gene annotation database by integrating 
#   multiple sources (Ensembl biomaRt, org.Hs.eg.db) and performs coordinate
#   liftover from hg19 to hg38 for Illumina HumanMethylation450 BeadChip probes.
#   The pipeline generates multiple gene-centric annotation tables grouped by
#   Entrez ID, Ensembl ID, and gene symbol, along with probe-to-gene mappings.
#
# Author: [Your Name]
# Date: 2026-01-01
# Version: 1.0
#
# Input:
#   - data/HM450.hg38.manifest.gencode.v36.tsv.gz (Illumina 450K manifest, hg38)
#   - data/hg19ToHg38.over2.chain (UCSC chain file for liftover)
#
# Output:
#   - Gene annotation tables: gene_info_final_2.6, gene_info_final_3.5, gene_info_final_4
#   - Probe annotation table: probeInfo
#
# Dependencies:
#   - biomaRt: Access to Ensembl gene annotation database
#   - org.Hs.eg.db: Human genome annotation from Bioconductor
#   - minfi: Illumina methylation array analysis
#   - rtracklayer: Genomic coordinate liftover
#   - dplyr, tidyr, purrr, stringr: Data manipulation
#   - data.table: Fast data processing
#   - IlluminaHumanMethylation450kanno.ilmn12.hg19: Array annotation package
#
################################################################################

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(data.table)
library(minfi)
library(rtracklayer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Define relative paths (assumes data files are in 'data/' subdirectory)
MANIFEST_HG38_PATH <- "data/HM450.hg38.manifest.gencode.v36.tsv.gz"
CHAIN_FILE_PATH <- "data/hg19ToHg38.over2.chain"

# Mirror selection for Ensembl (use US East mirror for stability)
ENSEMBL_MIRROR <- "useast"

# ==============================================================================
# PART 1: RETRIEVE GENE ANNOTATIONS FROM MULTIPLE SOURCES
# ==============================================================================

cat("=== PART 1: Retrieving gene annotations ===\n")

# ------------------------------------------------------------------------------
# 1.1 Query Ensembl biomaRt (Source 1: Basic gene information)
# ------------------------------------------------------------------------------
# Rationale: biomaRt provides the most current Ensembl gene IDs and symbols
# with chromosome locations. This serves as the primary annotation source.

cat("Connecting to Ensembl biomaRt...\n")
ensembl <- useEnsembl(
  biomart = "ensembl", 
  dataset = "hsapiens_gene_ensembl", 
  mirror = ENSEMBL_MIRROR
)

# Retrieve core gene attributes: Ensembl ID, gene symbol, Entrez ID, chromosome
attributes_basic <- c(
  "ensembl_gene_id", 
  "external_gene_name",
  "entrezgene_id",
  "chromosome_name"
)

cat("Querying basic gene information from biomaRt...\n")
gene_info_biomart_basic <- getBM(
  attributes = attributes_basic,
  mart = ensembl
)

# Standardize column names for downstream integration
gene_info_biomart_basic_clean <- data.frame(
  ensembl_gene_id = gene_info_biomart_basic$ensembl_gene_id,
  gene_name = gene_info_biomart_basic$external_gene_name,
  entrezgene_id = gene_info_biomart_basic$entrezgene_id,
  chromosome_name = gene_info_biomart_basic$chromosome_name,
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 1.2 Query Ensembl biomaRt (Source 2: Extended gene information)
# ------------------------------------------------------------------------------
# Rationale: Second query includes Entrez accessions which may provide 
# additional gene symbols not captured in external_gene_name field.

attributes_extended <- c(
  "ensembl_gene_id",
  "external_gene_name", 
  "entrezgene_accession",
  "entrezgene_id",
  "chromosome_name"
)

cat("Querying extended gene information from biomaRt...\n")
gene_info_biomart_extended <- getBM(
  attributes = attributes_extended,
  mart = ensembl
)

# ------------------------------------------------------------------------------
# 1.3 Query org.Hs.eg.db (Source 3: Bioconductor annotation)
# ------------------------------------------------------------------------------
# Rationale: org.Hs.eg.db provides curated Entrez-centric annotations including
# full gene descriptions and alternative gene-chromosome mappings that may
# supplement biomaRt data.

cat("Querying org.Hs.eg.db annotation database...\n")

# Extract all Entrez Gene IDs as keys
entrez_ids_all <- keys(org.Hs.eg.db, keytype = "ENTREZID")

# Map Entrez IDs to gene symbols
symbol_mapping <- as.list(org.Hs.egSYMBOL)
gene_symbols <- sapply(symbol_mapping[entrez_ids_all], `[`, 1)

# Map Entrez IDs to full gene names (descriptions)
genename_mapping <- as.list(org.Hs.egGENENAME)
gene_descriptions <- sapply(genename_mapping[entrez_ids_all], `[`, 1)

# Map Entrez IDs to chromosome locations
chr_mapping <- as.list(org.Hs.egCHR)
chr_locations <- sapply(chr_mapping[entrez_ids_all], `[`, 1)

# Map Entrez IDs to Ensembl IDs
ensembl_mapping <- as.list(org.Hs.egENSEMBL)
ensembl_ids <- sapply(ensembl_mapping[entrez_ids_all], `[`, 1)

# Construct comprehensive annotation data frame
gene_info_orgdb <- data.frame(
  ENTREZID = entrez_ids_all,
  SYMBOL = gene_symbols,
  GENENAME = gene_descriptions,
  CHR = chr_locations,
  ENSEMBL = ensembl_ids,
  stringsAsFactors = FALSE
)

# ==============================================================================
# PART 2: INTEGRATE GENE ANNOTATIONS WITH PRIORITIZATION
# ==============================================================================

cat("\n=== PART 2: Integrating annotations from multiple sources ===\n")

# ------------------------------------------------------------------------------
# 2.1 Standardize format across all sources
# ------------------------------------------------------------------------------
# Rationale: Different sources use different column names and data types.
# Standardization enables seamless integration and prioritization.

# Priority 1: biomaRt basic (most current, primary source)
source1_standardized <- gene_info_biomart_basic_clean %>%
  transmute(
    entrez_id = as.character(entrezgene_id),
    ensembl_id = ensembl_gene_id,
    gene_symbol = gene_name,
    gene_description = NA_character_,
    chromosome = chromosome_name,
    priority = 1
  )

# Priority 2: biomaRt extended (supplementary gene symbols from accessions)
source2_standardized <- gene_info_biomart_extended %>%
  transmute(
    entrez_id = as.character(entrezgene_id),
    ensembl_id = ensembl_gene_id,
    gene_symbol = entrezgene_accession,
    gene_description = NA_character_,
    chromosome = chromosome_name,
    priority = 2
  )

# Priority 3: org.Hs.eg.db (curated descriptions, alternative mappings)
source3_standardized <- gene_info_orgdb %>%
  transmute(
    entrez_id = as.character(ENTREZID),
    ensembl_id = ENSEMBL,
    gene_symbol = SYMBOL,
    gene_description = GENENAME,
    chromosome = CHR,
    priority = 3
  )

# ------------------------------------------------------------------------------
# 2.2 Combine all sources with priority ordering
# ------------------------------------------------------------------------------
# Rationale: By combining sources and maintaining priority, we ensure that
# more reliable data (biomaRt) takes precedence while still capturing
# comprehensive annotations from all sources.

cat("Merging annotation sources with prioritization...\n")
gene_annotations_combined <- bind_rows(
  source1_standardized,
  source2_standardized,
  source3_standardized
) %>%
  arrange(priority) %>%
  distinct()

cat("Combined annotation records:", nrow(gene_annotations_combined), "\n")

# ==============================================================================
# PART 3: CREATE GENE-CENTRIC ANNOTATION TABLES
# ==============================================================================

cat("\n=== PART 3: Creating gene-centric annotation tables ===\n")

# ------------------------------------------------------------------------------
# 3.1 Version 2.6: Entrez ID-centric annotation
# ------------------------------------------------------------------------------
# Rationale: Entrez Gene IDs are stable, widely-used identifiers in genomics.
# This version is optimal for analyses requiring Entrez-based gene references
# (e.g., pathway analysis, cross-database integration).

cat("Creating Entrez ID-centric annotation (v2.6)...\n")
gene_annotation_by_entrez <- gene_annotations_combined %>%
  group_by(entrez_id) %>%
  summarise(
    ensembl_id = first(ensembl_id[!is.na(ensembl_id) & ensembl_id != ""]),
    gene_symbol = first(gene_symbol[!is.na(gene_symbol) & gene_symbol != ""]),
    gene_description = first(gene_description[!is.na(gene_description) & gene_description != ""]),
    # Prioritize numeric chromosomes (1-22, X, Y) over scaffolds/patches
    chromosome = {
      chr_values <- chromosome[!is.na(chromosome) & chromosome != ""]
      chr_numeric <- chr_values[grepl("^[0-9]+$", chr_values)]
      if (length(chr_numeric) > 0) {
        first(chr_numeric)
      } else {
        first(chr_values)
      }
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(entrez_id)) %>%
  select(entrez_id, ensembl_id, gene_symbol, gene_description, chromosome)

# Record chromosome 17 count for QC
chr17_count_entrez <- sum(gene_annotation_by_entrez$chromosome == "17", na.rm = TRUE)
cat("Entrez-centric annotation created:", nrow(gene_annotation_by_entrez), "genes\n")
cat("Chromosome 17 genes:", chr17_count_entrez, "\n")

# ------------------------------------------------------------------------------
# 3.2 Version 3.5: Ensembl ID-centric annotation
# ------------------------------------------------------------------------------
# Rationale: Ensembl IDs provide the most current gene models and are essential
# for genomic coordinate-based analyses. This version is optimal for genomic
# region analyses and Ensembl-based workflows.

cat("Creating Ensembl ID-centric annotation (v3.5)...\n")
gene_annotation_by_ensembl <- gene_annotations_combined %>%
  filter(!is.na(ensembl_id) & str_trim(ensembl_id) != "") %>%
  group_by(ensembl_id) %>%
  summarise(
    gene_symbol = first(gene_symbol[!is.na(gene_symbol) & gene_symbol != ""]),
    entrez_id = first(entrez_id[!is.na(entrez_id) & entrez_id != ""]),
    gene_description = first(gene_description[!is.na(gene_description) & gene_description != ""]),
    chromosome = {
      chr_values <- chromosome[!is.na(chromosome) & chromosome != ""]
      chr_numeric <- chr_values[grepl("^[0-9]+$", chr_values)]
      if (length(chr_numeric) > 0) {
        first(chr_numeric)
      } else {
        first(chr_values)
      }
    },
    .groups = "drop"
  ) %>%
  select(ensembl_id, gene_symbol, entrez_id, gene_description, chromosome)

chr17_count_ensembl <- sum(gene_annotation_by_ensembl$chromosome == "17", na.rm = TRUE)
cat("Ensembl-centric annotation created:", nrow(gene_annotation_by_ensembl), "genes\n")
cat("Chromosome 17 genes:", chr17_count_ensembl, "\n")

# ------------------------------------------------------------------------------
# 3.3 Version 4: Gene symbol-centric annotation
# ------------------------------------------------------------------------------
# Rationale: Gene symbols are human-readable and most commonly used in
# publications. This version is optimal for literature-based analyses and
# result interpretation where symbol-based lookups are required.

cat("Creating gene symbol-centric annotation (v4)...\n")
gene_annotation_by_symbol <- gene_annotations_combined %>%
  filter(!is.na(gene_symbol) & str_trim(gene_symbol) != "") %>%
  group_by(gene_symbol) %>%
  summarise(
    ensembl_id = first(na.omit(ensembl_id)),
    entrez_id = first(na.omit(entrez_id)),
    gene_description = first(na.omit(gene_description)),
    chromosome = {
      chr_values <- na.omit(chromosome)
      chr_numeric <- chr_values[grepl("^[0-9]+$", chr_values)]
      if (length(chr_numeric) > 0) {
        first(chr_numeric)
      } else {
        first(chr_values)
      }
    },
    .groups = "drop"
  ) %>%
  select(gene_symbol, ensembl_id, entrez_id, gene_description, chromosome)

chr17_count_symbol <- sum(gene_annotation_by_symbol$chromosome == "17", na.rm = TRUE)
cat("Symbol-centric annotation created:", nrow(gene_annotation_by_symbol), "genes\n")
cat("Chromosome 17 genes:", chr17_count_symbol, "\n")

# ==============================================================================
# PART 4: COORDINATE LIFTOVER FROM HG19 TO HG38
# ==============================================================================

cat("\n=== PART 4: Performing coordinate liftover (hg19 -> hg38) ===\n")

# ------------------------------------------------------------------------------
# 4.1 Load hg19 annotation and hg38 manifest
# ------------------------------------------------------------------------------
# Rationale: Illumina's original 450K annotation uses hg19 coordinates.
# Many current analyses require hg38 coordinates. Liftover ensures
# compatibility with modern reference genomes.

cat("Loading hg38 manifest data...\n")
manifest_hg38 <- read.delim(MANIFEST_HG38_PATH)

cat("Loading hg19 annotation data...\n")
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_hg19 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annotation_hg19_df <- as.data.frame(annotation_hg19)

# ------------------------------------------------------------------------------
# 4.2 Perform coordinate liftover using UCSC chain file
# ------------------------------------------------------------------------------
# Rationale: UCSC chain files provide validated coordinate mappings between
# genome assemblies. The liftOver process accounts for assembly differences
# including insertions, deletions, and rearrangements.

cat("Loading UCSC chain file for liftover...\n")
chain_hg19_to_hg38 <- import.chain(CHAIN_FILE_PATH)

cat("Creating GRanges object from hg19 CpG coordinates...\n")
cpg_granges_hg19 <- GRanges(
  seqnames = annotation_hg19$chr,
  ranges = IRanges(start = annotation_hg19$pos, width = 1),
  strand = annotation_hg19$strand,
  name = rownames(annotation_hg19)
)

cat("Performing coordinate liftover to hg38...\n")
cpg_granges_hg38 <- liftOver(cpg_granges_hg19, chain_hg19_to_hg38)
cpg_granges_hg38 <- unlist(cpg_granges_hg38)

# ------------------------------------------------------------------------------
# 4.3 Merge hg38 coordinates with hg19 annotations
# ------------------------------------------------------------------------------

cat("Merging hg38 coordinates with hg19 annotations...\n")
cpg_coordinates_hg38 <- as.data.frame(cpg_granges_hg38)
cpg_coordinates_hg38$CpG <- cpg_granges_hg38$name
annotation_hg19_df$CpG <- rownames(annotation_hg19_df)

annotation_hg38_merged <- cpg_coordinates_hg38 %>%
  left_join(annotation_hg19_df, by = "CpG")

cat("Merged annotation contains", nrow(annotation_hg38_merged), "CpG sites\n")

# ==============================================================================
# PART 5: CREATE ONE-TO-ONE PROBE-TO-GENE MAPPING
# ==============================================================================

cat("\n=== PART 5: Creating one-to-one probe-to-gene mapping ===\n")

# ------------------------------------------------------------------------------
# 5.1 Define nearest gene selection function
# ------------------------------------------------------------------------------
# Rationale: Many probes map to multiple genes. To create unambiguous
# probe-gene relationships, we select the nearest gene based on distance
# to transcription start site (TSS). This is biologically justified as
# CpG methylation typically affects the nearest gene's regulation.

pick_nearest_gene_to_probe <- function(gene_string, distance_string) {
  # Parse semicolon-separated gene names and distances
  genes <- unlist(strsplit(gene_string, ";", fixed = TRUE))
  distances <- suppressWarnings(
    as.numeric(unlist(strsplit(distance_string, ";", fixed = TRUE)))
  )
  
  # Validate input
  if (length(genes) == 0 || length(distances) == 0) {
    return(NULL)
  }
  
  # Strategy: Keep all genes regardless of distance (no arbitrary cutoff)
  # This maximizes probe coverage, especially for intergenic regions
  keep_indices <- seq_along(genes)
  
  if (length(keep_indices) == 0) {
    return(NULL)
  }
  
  genes <- genes[keep_indices]
  distances <- distances[keep_indices]
  
  # Find gene(s) with minimum absolute distance to TSS
  min_abs_distance <- min(abs(distances))
  candidate_indices <- which(abs(distances) == min_abs_distance)
  
  # Tie-breaking rule: If multiple genes equidistant, prefer upstream genes
  # (negative distance) as they are more likely to be regulatory targets
  if (length(candidate_indices) > 1) {
    upstream_indices <- candidate_indices[distances[candidate_indices] < 0]
    if (length(upstream_indices) > 0) {
      selected_index <- upstream_indices[1]
    } else {
      selected_index <- candidate_indices[1]
    }
  } else {
    selected_index <- candidate_indices
  }
  
  data.frame(
    best_gene = genes[selected_index], 
    best_distToTSS = distances[selected_index]
  )
}

# ------------------------------------------------------------------------------
# 5.2 Apply gene selection to all probes
# ------------------------------------------------------------------------------

cat("Applying nearest gene selection to all probes...\n")
probe_to_gene_mapping <- manifest_hg38 %>%
  mutate(
    gene_selection = pmap(
      list(geneNames, distToTSS),
      ~ pick_nearest_gene_to_probe(..1, ..2)
    )
  ) %>%
  filter(!map_lgl(gene_selection, is.null)) %>%
  unnest(cols = c(gene_selection)) %>%
  select(probeID, CpG_chrm, best_gene, best_distToTSS)

# Standardize column names and chromosome notation (remove 'chr' prefix)
probe_to_gene_mapping <- probe_to_gene_mapping %>%
  rename(
    chromosome_name = CpG_chrm,
    gene_name = best_gene
  ) %>%
  mutate(chromosome_name = gsub("^chr", "", chromosome_name))

cat("One-to-one probe-gene mapping created:", nrow(probe_to_gene_mapping), "probes\n")

# ==============================================================================
# PART 6: BUILD COMPREHENSIVE PROBE INFORMATION TABLE
# ==============================================================================

cat("\n=== PART 6: Building comprehensive probe annotation table ===\n")

# ------------------------------------------------------------------------------
# 6.1 Integrate probe-gene mapping with UCSC RefGene annotations
# ------------------------------------------------------------------------------
# Rationale: Combining our nearest-gene mapping with UCSC RefGene group
# annotations (TSS, 5'UTR, Body, 3'UTR, etc.) provides complete functional
# context for each probe's genomic location.

comprehensive_probe_info <- bind_rows(
  # Source 1: Our one-to-one gene mapping
  probe_to_gene_mapping %>%
    select(probeID, gene_name, chromosome_name),
  
  # Source 2: UCSC RefGene functional region annotations
  annotation_hg38_merged %>%
    rename(probeID = CpG) %>%
    select(probeID, UCSC_RefGene_Group)
) %>%
  group_by(probeID) %>%
  summarise(
    gene_name = first(na.omit(gene_name)),
    UCSC_RefGene_Group = first(na.omit(UCSC_RefGene_Group)),
    chromosome_name = first(na.omit(chromosome_name)),
    .groups = "drop"
  )

cat("Comprehensive probe annotation created:", nrow(comprehensive_probe_info), "probes\n")

# Verify one-to-one relationship (each probe should appear exactly once)
is_one_to_one <- nrow(comprehensive_probe_info) == n_distinct(comprehensive_probe_info$probeID)
cat("One-to-one probe mapping verified:", is_one_to_one, "\n")

# ==============================================================================
# FINAL OUTPUT SUMMARY
# ==============================================================================

cat("\n=== PIPELINE COMPLETE ===\n")
cat("\nGenerated annotation tables:\n")
cat("1. gene_annotation_by_entrez (gene_info_final_2.6):", nrow(gene_annotation_by_entrez), "genes\n")
cat("2. gene_annotation_by_ensembl (gene_info_final_3.5):", nrow(gene_annotation_by_ensembl), "genes\n")
cat("3. gene_annotation_by_symbol (gene_info_final_4):", nrow(gene_annotation_by_symbol), "genes\n")
cat("4. comprehensive_probe_info (probeInfo):", nrow(comprehensive_probe_info), "probes\n")

# Assign to legacy variable names for backward compatibility
gene_info_final_2.6 <- gene_annotation_by_entrez
gene_info_final_3.5 <- gene_annotation_by_ensembl
gene_info_final_4 <- gene_annotation_by_symbol
probeInfo <- comprehensive_probe_info

cat("\nAll annotation tables are ready for downstream analysis.\n")

################################################################################
# END OF SCRIPT
################################################################################
