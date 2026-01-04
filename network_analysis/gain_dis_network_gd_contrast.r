# ==============================================================================
# BIOLOGICAL NETWORK ANALYSIS: CHROMOSOME 17 COPY NUMBER VARIATION
# ==============================================================================
# Title: Network Construction for GAIN vs DIS (disomic) Contrast
# Description: Constructs gene regulatory networks from expression and methylation
#              data, focusing on chromosome 17 genes showing significant changes
#              between GAIN (amplification) and DIS (disomic) conditions.
#              This analysis identifies molecular changes specifically associated
#              with chromosome 17 copy number gain, independent of BRCA1 status.
#              Uses stricter thresholds (2.3 for expression, 2.4 for methylation)
#              to focus on high-confidence relationships.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Biological Context:
#   Chromosome 17 copy number gain (GAIN) is common in breast cancer and can
#   affect numerous genes including BRCA1, HER2/ERBB2, and others. This analysis
#   compares samples with chr17 gain against disomic (DIS) samples to
#   identify dosage-sensitive genes and regulatory networks affected by
#   chromosomal amplification.
#
# Input Files:
#   - data/expression/GAIN/all_relations_gainexp_2.3.txt
#   - data/expression/DIS/all_relations_disexp_2.3.txt
#   - data/methylation/GAIN/all_relations_gainmet_2.4.txt
#   - data/methylation/DIS/all_relations_dismet_2.4.txt
#   - data/reference/gene_info_final_3.5.rds (Ensembl to gene symbol mapping)
#   - data/reference/gene_info_final_4.rds (Gene annotations)
#   - data/reference/probeInfo.rds (Methylation probe annotations)
#   - data/deg/gain_vs_dis_deg.rds (Differential expression results)
#   - data/dmp/gain_vs_dis_dmp.rds (Differential methylation results)
#
# Output: Network nodes and edges data frames for Cytoscape
#
# Dependencies: dplyr, stringr, purrr
# ==============================================================================

# === LOAD REQUIRED LIBRARIES ===
library(dplyr)
library(stringr)
library(purrr)

# === LOAD REFERENCE DATA ===
# Gene annotation files for mapping IDs and adding biological context
gene_info_final_3.5 <- readRDS("data/reference/gene_info_final_3.5.rds")
gene_info_final_4 <- readRDS("data/reference/gene_info_final_4.rds")
probeInfo <- readRDS("data/reference/probeInfo.rds")

# Differential analysis results for status annotation
gain_vs_dis_deg <- readRDS("data/deg/gain_vs_dis_deg.rds")
gain_vs_dis_dmp <- readRDS("data/dmp/gain_vs_dis_dmp.rds")


# ==============================================================================
# PART 1: EXPRESSION NETWORK - GAIN Samples (Chr17 Focus)
# ==============================================================================
# Rationale: Focus on chromosome 17 where copy number variation occurs. Use
# stricter threshold (2.3) to identify high-confidence regulatory relationships.
# Only include genes showing significant differential expression (|logFC| > 1)
# to capture dosage-sensitive genes responding to chromosomal amplification.

# --- Load Expression Relations Data ---
exp_gain_chr17_sthr2.3 <- read.table(
  "data/expression/GAIN/all_relations_gainexp_2.3.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Genes ---
# Create master list of all genes involved in regulatory relationships
unique_genes_gain_df <- data.frame(
  ensembl_gene_id = unique(c(
    exp_gain_chr17_sthr2.3$source_ensembl_id,
    exp_gain_chr17_sthr2.3$target_ensembl_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ensembl_gene_id))

# --- Node Construction ---
# Build comprehensive node table with multi-layer annotations:
# 1. Map Ensembl IDs to gene symbols and chromosomes
# 2. Add gene descriptions and Entrez IDs for external database linking
# 3. Annotate expression status (OVER_T/DOWN_T/UNCHANGED) based on DEG results
#
# Status nomenclature difference from BRCA1 analysis:
# - OVER_T (overexpressed in tumor/GAIN): logFC > 1
# - DOWN_T (underexpressed in tumor/GAIN): logFC < -1
# - UNCHANGED: |logFC| <= 1
nodes_exp_gain_chr17_sthr2.3 <- unique_genes_gain_df %>%
  # Layer 1: Basic gene identity (Ensembl -> Symbol + Chromosome)
  left_join(gene_info_final_3.5, by = c("ensembl_gene_id" = "ensembl_id")) %>%
  
  # Layer 2: Detailed annotations (Symbol -> Description + Entrez)
  left_join(
    gene_info_final_4 %>%
      select(gene_symbol, gene_description, entrez_id) %>%
      rename(
        gene_symbol_lookup = gene_symbol,
        gene_description_lookup = gene_description,
        entrez_id_lookup = entrez_id
      ),
    by = c("gene_symbol" = "gene_symbol_lookup")
  ) %>%
  
  # Layer 3: Expression status classification
  # Using GAIN vs DIS contrast to identify dosage-sensitive genes
  left_join(
    gain_vs_dis_deg %>%
      mutate(
        expression_status = case_when(
          logFC > 1  ~ "OVER_T",
          logFC < -1 ~ "DOWN_T",
          TRUE       ~ NA_character_
        )
      ) %>%
      select(gene_symbol, expression_status),
    by = "gene_symbol"
  ) %>%
  
  # Finalize node attributes
  mutate(
    id         = gene_symbol,
    symbol     = gene_description_lookup,
    entrez_id  = entrez_id_lookup,
    ensembl_id = ensembl_gene_id,
    chromosome = chromosome,
    status     = if_else(is.na(expression_status), "UNCHANGED", expression_status),
    source     = "TCGA-RNAseq",
    type       = "gene",
    layer      = "expression"
  ) %>%
  filter(!is.na(id), id != "") %>%
  select(id, symbol, entrez_id, ensembl_id, type, layer,
         chromosome, status, source)

# Quality check: Verify chromosome and status distributions
message("Expression GAIN - Chromosome distribution:")
print(table(nodes_exp_gain_chr17_sthr2.3$chromosome))
message("Expression GAIN - Status distribution:")
print(table(nodes_exp_gain_chr17_sthr2.3$status))

# Add layer-specific suffix to prevent ID collisions across omics layers
nodes_exp_gain_chr17_sthr2.3 <- nodes_exp_gain_chr17_sthr2.3 %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
# Build regulatory relationship table with stricter threshold (2.3)
rels_exp_gain_chr17_sthr2.3 <- exp_gain_chr17_sthr2.3 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "expression"
  ) %>%
  distinct(from, to, .keep_all = TRUE)  # Remove duplicate edges

# Add layer suffix to match node IDs
rels_exp_gain_chr17_sthr2.3 <- rels_exp_gain_chr17_sthr2.3 %>%
  mutate(
    from = paste0(from, "_expression"),
    to   = paste0(to, "_expression")
  )

# --- Filter for Biologically Relevant Edges ---
# Rationale: Focus only on edges between genes showing significant expression
# changes (OVER_T or DOWN_T). This removes noise from unchanged genes and
# highlights the core regulatory network responding to chr17 copy number gain.
# Note: Very low count suggests most genes don't meet strict threshold in both
# directions, indicating sparse high-confidence network
rels_exp_gain_GD_sthr2.3 <- rels_exp_gain_chr17_sthr2.3 %>%
  left_join(
    nodes_exp_gain_chr17_sthr2.3 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_gain_chr17_sthr2.3 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  # Both nodes must be differentially expressed
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
# Rationale: "Equivalent" and "Opposite" relationships are symmetric (A-B = B-A).
# Keep only one direction to avoid redundancy in undirected network visualization.
rels_exp_gain_GD_equiv_opp <- rels_exp_gain_GD_sthr2.3 %>%
  filter(relation_type %in% c("Equivalent", "Opposite")) %>%
  mutate(
    from_base = str_remove(from, "_.*$"),
    to_base   = str_remove(to, "_.*$"),
    # Create canonical pair key (sorted alphabetically)
    pair_key  = map2_chr(from_base, to_base, ~ paste(sort(c(.x, .y)), collapse = "_"))
  ) %>%
  group_by(pair_key) %>%
  slice(1) %>%  # Keep first occurrence only
  ungroup() %>%
  select(from, to, relation_type, layer)

# Keep all other (directional) relations
rels_exp_gain_GD_others <- rels_exp_gain_GD_sthr2.3 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

# Merge deduplicated edges back
rels_exp_gain_GD_sthr2.3 <- bind_rows(
  rels_exp_gain_GD_others,
  rels_exp_gain_GD_equiv_opp
)

# Summary statistics
message("Expression GAIN GD - Relation type distribution:")
print(table(rels_exp_gain_GD_sthr2.3$relation_type))

# --- Network Integrity Validation ---
# Check for orphaned nodes (nodes without any edges)
orphaned_nodes <- setdiff(
  nodes_exp_gain_chr17_sthr2.3$id,
  c(rels_exp_gain_GD_sthr2.3$from, rels_exp_gain_GD_sthr2.3$to)
)
message("Expression GAIN - Orphaned nodes: ", length(orphaned_nodes))

# Check for edges referencing non-existent nodes
missing_from_nodes <- setdiff(
  rels_exp_gain_GD_sthr2.3$from,
  nodes_exp_gain_chr17_sthr2.3$id
)
missing_to_nodes <- setdiff(
  rels_exp_gain_GD_sthr2.3$to,
  nodes_exp_gain_chr17_sthr2.3$id
)
message("Expression GAIN - Missing FROM nodes: ", length(missing_from_nodes))
message("Expression GAIN - Missing TO nodes: ", length(missing_to_nodes))


# ==============================================================================
# PART 2: EXPRESSION NETWORK - DIS Samples (Chr17 Focus)
# ==============================================================================
# Rationale: Analyze disomic (DIS) samples to establish baseline
# regulatory network. Comparison with GAIN network identifies changes
# specifically due to chromosomal amplification.

# --- Load Expression Relations Data ---
exp_dis_chr17_sthr2.3 <- read.table(
  "data/expression/DIS/all_relations_disexp_2.3.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Genes ---
unique_genes_dis_df <- data.frame(
  ensembl_gene_id = unique(c(
    exp_dis_chr17_sthr2.3$source_ensembl_id,
    exp_dis_chr17_sthr2.3$target_ensembl_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ensembl_gene_id))

# --- Node Construction ---
# Same process as GAIN, but using DIS-specific gene list
nodes_exp_dis_chr17_sthr2.3 <- unique_genes_dis_df %>%
  left_join(gene_info_final_3.5, by = c("ensembl_gene_id" = "ensembl_id")) %>%
  left_join(
    gene_info_final_4 %>%
      select(gene_symbol, gene_description, entrez_id) %>%
      rename(
        gene_symbol_lookup = gene_symbol,
        gene_description_lookup = gene_description,
        entrez_id_lookup = entrez_id
      ),
    by = c("gene_symbol" = "gene_symbol_lookup")
  ) %>%
  # Use same GAIN vs DIS contrast for status classification
  left_join(
    gain_vs_dis_deg %>%
      mutate(
        expression_status = case_when(
          logFC > 1  ~ "OVER_T",
          logFC < -1 ~ "DOWN_T",
          TRUE       ~ NA_character_
        )
      ) %>%
      select(gene_symbol, expression_status),
    by = "gene_symbol"
  ) %>%
  mutate(
    id         = gene_symbol,
    symbol     = gene_description_lookup,
    entrez_id  = entrez_id_lookup,
    ensembl_id = ensembl_gene_id,
    chromosome = chromosome,
    status     = if_else(is.na(expression_status), "UNCHANGED", expression_status),
    source     = "TCGA-RNAseq",
    type       = "gene",
    layer      = "expression"
  ) %>%
  filter(!is.na(id), id != "") %>%
  select(id, symbol, entrez_id, ensembl_id, type, layer,
         chromosome, status, source)

message("Expression DIS - Chromosome distribution:")
print(table(nodes_exp_dis_chr17_sthr2.3$chromosome))
message("Expression DIS - Status distribution:")
print(table(nodes_exp_dis_chr17_sthr2.3$status))

nodes_exp_dis_chr17_sthr2.3 <- nodes_exp_dis_chr17_sthr2.3 %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
rels_exp_dis_chr17_sthr2.3 <- exp_dis_chr17_sthr2.3 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "expression"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_exp_dis_chr17_sthr2.3 <- rels_exp_dis_chr17_sthr2.3 %>%
  mutate(
    from = paste0(from, "_expression"),
    to   = paste0(to, "_expression")
  )

# --- Filter for Biologically Relevant Edges ---
# Note: Extremely sparse network suggests DIS samples have very few genes
# meeting both the strict threshold and differential expression criteria
rels_exp_dis_GD_sthr2.3 <- rels_exp_dis_chr17_sthr2.3 %>%
  left_join(
    nodes_exp_dis_chr17_sthr2.3 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_dis_chr17_sthr2.3 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_exp_dis_GD_equiv_opp <- rels_exp_dis_GD_sthr2.3 %>%
  filter(relation_type %in% c("Equivalent", "Opposite")) %>%
  mutate(
    from_base = str_remove(from, "_.*$"),
    to_base   = str_remove(to, "_.*$"),
    pair_key  = map2_chr(from_base, to_base, ~ paste(sort(c(.x, .y)), collapse = "_"))
  ) %>%
  group_by(pair_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(from, to, relation_type, layer)

rels_exp_dis_GD_others <- rels_exp_dis_GD_sthr2.3 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_exp_dis_GD_sthr2.3 <- bind_rows(
  rels_exp_dis_GD_others,
  rels_exp_dis_GD_equiv_opp
)

message("Expression DIS GD - Relation type distribution:")
print(table(rels_exp_dis_GD_sthr2.3$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_exp_dis_chr17_sthr2.3$id,
  c(rels_exp_dis_GD_sthr2.3$from, rels_exp_dis_GD_sthr2.3$to)
)
missing_from_nodes <- setdiff(
  rels_exp_dis_GD_sthr2.3$from,
  nodes_exp_dis_chr17_sthr2.3$id
)
missing_to_nodes <- setdiff(
  rels_exp_dis_GD_sthr2.3$to,
  nodes_exp_dis_chr17_sthr2.3$id
)
message("Expression DIS - Orphaned nodes: ", length(orphaned_nodes))
message("Expression DIS - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 3: METHYLATION NETWORK - GAIN Samples (Chr17 Focus)
# ==============================================================================
# Rationale: Methylation patterns can be affected by copy number alterations.
# Focus on CpG sites showing significant differential methylation (|logFC| > 0.2)
# in GAIN vs DIS contrast. Uses stricter threshold (2.4) for high-confidence
# methylation relationships.

# --- Load Methylation Relations Data ---
met_gain_chr17_sthr2.4 <- read.table(
  "data/methylation/GAIN/all_relations_gainmet_2.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Probes ---
# Methylation uses probe IDs rather than gene IDs directly
unique_probes_met_df <- data.frame(
  probeID = unique(c(
    met_gain_chr17_sthr2.4$source_probe_id,
    met_gain_chr17_sthr2.4$target_probe_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(probeID))

# --- Node Construction ---
# Multi-step annotation process for methylation probes
nodes_met_gain_chr17_sthr2.4 <- unique_probes_met_df %>%
  # Layer 1: Probe annotations (probe -> gene + location)
  left_join(probeInfo, by = "probeID") %>%
  
  # Layer 2: Methylation status classification
  # HYPERMETH: logFC > 0.2 (increased methylation in GAIN)
  # HYPOMETH: logFC < -0.2 (decreased methylation in GAIN)
  left_join(
    gain_vs_dis_dmp %>%
      mutate(
        methyl_status = case_when(
          logFC > 0.2  ~ "HYPERMETH",
          logFC < -0.2 ~ "HYPOMETH",
          TRUE         ~ NA_character_
        )
      ) %>%
      select(gene_name, methyl_status),
    by = "gene_name"
  ) %>%
  
  # Layer 3: Gene annotations
  left_join(
    gene_info_final_4 %>%
      select(gene_symbol, gene_description, entrez_id) %>%
      rename(
        gene_name = gene_symbol,
        symbol    = gene_description,
        entrez_id = entrez_id
      ),
    by = "gene_name"
  ) %>%
  
  mutate(
    id = gene_name,
    chromosome = if_else(
      is.na(chromosome_name) | chromosome_name == "",
      NA_character_,
      chromosome_name
    ),
    cpg_location = if_else(
      is.na(UCSC_RefGene_Group) | UCSC_RefGene_Group == "",
      NA_character_,
      UCSC_RefGene_Group
    ),
    probe_id = probeID,
    status   = if_else(is.na(methyl_status), "UNCHANGED", methyl_status),
    source   = "TCGA-Methylation",
    type     = "gene",
    layer    = "methylation"
  ) %>%
  filter(!is.na(id), id != "") %>%
  select(id, symbol, entrez_id, probe_id, type, layer,
         chromosome, cpg_location, status, source)

message("Methylation GAIN - Chromosome distribution:")
print(table(nodes_met_gain_chr17_sthr2.4$chromosome))
message("Methylation GAIN - Status distribution:")
print(table(nodes_met_gain_chr17_sthr2.4$status))

nodes_met_gain_chr17_sthr2.4 <- nodes_met_gain_chr17_sthr2.4 %>%
  mutate(id = paste0(id, "_methylation"))

# --- Edge Construction ---
rels_met_gain_chr17_sthr2.4 <- met_gain_chr17_sthr2.4 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "methylation"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_met_gain_chr17_sthr2.4 <- rels_met_gain_chr17_sthr2.4 %>%
  mutate(
    from = paste0(from, "_methylation"),
    to   = paste0(to, "_methylation")
  )

# --- Filter for Biologically Relevant Edges ---
rels_met_gain_GD_sthr2.4 <- rels_met_gain_chr17_sthr2.4 %>%
  left_join(
    nodes_met_gain_chr17_sthr2.4 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_gain_chr17_sthr2.4 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_gain_GD_equiv_opp <- rels_met_gain_GD_sthr2.4 %>%
  filter(relation_type %in% c("Equivalent", "Opposite")) %>%
  mutate(
    from_base = str_remove(from, "_.*$"),
    to_base   = str_remove(to, "_.*$"),
    pair_key  = map2_chr(from_base, to_base, ~ paste(sort(c(.x, .y)), collapse = "_"))
  ) %>%
  group_by(pair_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(from, to, relation_type, layer)

rels_met_gain_GD_others <- rels_met_gain_GD_sthr2.4 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_gain_GD_sthr2.4 <- bind_rows(
  rels_met_gain_GD_others,
  rels_met_gain_GD_equiv_opp
)

message("Methylation GAIN GD - Relation type distribution:")
print(table(rels_met_gain_GD_sthr2.4$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_met_gain_chr17_sthr2.4$id,
  c(rels_met_gain_GD_sthr2.4$from, rels_met_gain_GD_sthr2.4$to)
)
missing_from_nodes <- setdiff(
  rels_met_gain_GD_sthr2.4$from,
  nodes_met_gain_chr17_sthr2.4$id
)
missing_to_nodes <- setdiff(
  rels_met_gain_GD_sthr2.4$to,
  nodes_met_gain_chr17_sthr2.4$id
)
message("Methylation GAIN - Orphaned nodes: ", length(orphaned_nodes))
message("Methylation GAIN - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 4: METHYLATION NETWORK - DIS Samples (Chr17 Focus)
# ==============================================================================

# --- Load Methylation Relations Data ---
met_dis_chr17_sthr2.4 <- read.table(
  "data/methylation/DIS/all_relations_dismet_2.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Probes ---
unique_probes_met_dis_df <- data.frame(
  probeID = unique(c(
    met_dis_chr17_sthr2.4$source_probe_id,
    met_dis_chr17_sthr2.4$target_probe_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(probeID))

# --- Node Construction ---
nodes_met_dis_chr17_sthr2.4 <- unique_probes_met_dis_df %>%
  left_join(probeInfo, by = "probeID") %>%
  left_join(
    gain_vs_dis_dmp %>%
      mutate(
        methyl_status = case_when(
          logFC > 0.2  ~ "HYPERMETH",
          logFC < -0.2 ~ "HYPOMETH",
          TRUE         ~ NA_character_
        )
      ) %>%
      select(gene_name, methyl_status),
    by = "gene_name"
  ) %>%
  left_join(
    gene_info_final_4 %>%
      select(gene_symbol, gene_description, entrez_id) %>%
      rename(
        gene_name = gene_symbol,
        symbol    = gene_description,
        entrez_id = entrez_id
      ),
    by = "gene_name"
  ) %>%
  mutate(
    id = gene_name,
    chromosome = if_else(
      is.na(chromosome_name) | chromosome_name == "",
      NA_character_,
      chromosome_name
    ),
    cpg_location = if_else(
      is.na(UCSC_RefGene_Group) | UCSC_RefGene_Group == "",
      NA_character_,
      UCSC_RefGene_Group
    ),
    probe_id = probeID,
    status   = if_else(is.na(methyl_status), "UNCHANGED", methyl_status),
    source   = "TCGA-Methylation",
    type     = "gene",
    layer    = "methylation"
  ) %>%
  filter(!is.na(id), id != "") %>%
  select(id, symbol, entrez_id, probe_id, type, layer,
         chromosome, cpg_location, status, source)

message("Methylation DIS - Chromosome distribution:")
print(table(nodes_met_dis_chr17_sthr2.4$chromosome))
message("Methylation DIS - Status distribution:")
print(table(nodes_met_dis_chr17_sthr2.4$status))

nodes_met_dis_chr17_sthr2.4 <- nodes_met_dis_chr17_sthr2.4 %>%
  mutate(id = paste0(id, "_methylation"))

# --- Edge Construction ---
rels_met_dis_chr17_sthr2.4 <- met_dis_chr17_sthr2.4 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "methylation"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_met_dis_chr17_sthr2.4 <- rels_met_dis_chr17_sthr2.4 %>%
  mutate(
    from = paste0(from, "_methylation"),
    to   = paste0(to, "_methylation")
  )

# --- Filter for Biologically Relevant Edges ---
rels_met_dis_GD_sthr2.4 <- rels_met_dis_chr17_sthr2.4 %>%
  left_join(
    nodes_met_dis_chr17_sthr2.4 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_dis_chr17_sthr2.4 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_dis_GD_equiv_opp <- rels_met_dis_GD_sthr2.4 %>%
  filter(relation_type %in% c("Equivalent", "Opposite")) %>%
  mutate(
    from_base = str_remove(from, "_.*$"),
    to_base   = str_remove(to, "_.*$"),
    pair_key  = map2_chr(from_base, to_base, ~ paste(sort(c(.x, .y)), collapse = "_"))
  ) %>%
  group_by(pair_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(from, to, relation_type, layer)

rels_met_dis_GD_others <- rels_met_dis_GD_sthr2.4 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_dis_GD_sthr2.4 <- bind_rows(
  rels_met_dis_GD_others,
  rels_met_dis_GD_equiv_opp
)

message("Methylation DIS GD - Relation type distribution:")
print(table(rels_met_dis_GD_sthr2.4$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_met_dis_chr17_sthr2.4$id,
  c(rels_met_dis_GD_sthr2.4$from, rels_met_dis_GD_sthr2.4$to)
)
missing_from_nodes <- setdiff(
  rels_met_dis_GD_sthr2.4$from,
  nodes_met_dis_chr17_sthr2.4$id
)
missing_to_nodes <- setdiff(
  rels_met_dis_GD_sthr2.4$to,
  nodes_met_dis_chr17_sthr2.4$id
)
message("Methylation DIS - Orphaned nodes: ", length(orphaned_nodes))
message("Methylation DIS - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
# Output objects available in workspace:
# Expression networks:
# - nodes_exp_gain_chr17_sthr2.3, rels_exp_gain_GD_sthr2.3 (GAIN samples)
# - nodes_exp_dis_chr17_sthr2.3, rels_exp_dis_GD_sthr2.3 (DIS samples)
#
# Methylation networks:
# - nodes_met_gain_chr17_sthr2.4, rels_met_gain_GD_sthr2.4 (GAIN samples)
# - nodes_met_dis_chr17_sthr2.4, rels_met_dis_GD_sthr2.4 (DIS samples)
#
# Unique gene/probe lists for downstream analysis:
# - unique_genes_gain_df, unique_genes_dis_df
# - unique_probes_met_df, unique_probes_met_dis_df
# ==============================================================================
