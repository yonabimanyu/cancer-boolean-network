# ==============================================================================
# BIOLOGICAL NETWORK ANALYSIS: BRCA1 EXPRESSION AND METHYLATION
# ==============================================================================
# Title: Network Construction for BRCA1 Upregulated vs Downregulated Contrast
# Description: Constructs gene regulatory networks from expression and methylation
#              data, focusing on chromosome 17 genes that show significant changes
#              between upBRCA1 and downBRCA1 conditions. Filters for biologically
#              relevant relationships (UP_REG/DOWN_REG for expression, 
#              HYPERMETH/HYPOMETH for methylation) and removes redundant edges.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Input Files:
#   - data/expression/upBRCA1/all_relations_upexp_1.4.txt
#   - data/expression/downBRCA1/all_relations_downexp_1.4.txt
#   - data/methylation/upBRCA1/all_relations_upmet_1.4.txt
#   - data/methylation/downBRCA1/all_relations_downmet_1.4.txt
#   - data/reference/gene_info_final_3.5.rds (Ensembl to gene symbol mapping)
#   - data/reference/gene_info_final_4.rds (Gene annotations)
#   - data/reference/probeInfo.rds (Methylation probe annotations)
#   - data/deg/up_down_brca1_deg.rds (Differential expression results)
#   - data/dmp/upbrca1_vs_downbrca1_dmp.rds (Differential methylation results)
#
# Output: Network nodes and edges data frames for Cytoscape/igraph
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
up_down_brca1_deg <- readRDS("data/deg/up_down_brca1_deg.rds")
upbrca1_vs_downbrca1_dmp <- readRDS("data/dmp/upbrca1_vs_downbrca1_dmp.rds")


# ==============================================================================
# PART 1: EXPRESSION NETWORK - upBRCA1 (Chr17 Focus)
# ==============================================================================
# Rationale: Focus on chromosome 17 where BRCA1 is located. Only include genes
# showing significant differential expression (|logFC| > 1) to capture 
# biologically meaningful regulatory relationships.

# --- Load Expression Relations Data ---
exp_upBRCA1_chr17_sthr1.40 <- read.table(
  "data/expression/upBRCA1/all_relations_upexp_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Genes ---
# Create master list of all genes involved in regulatory relationships
# Expected count: 702 unique genes (as of 2025-11-19)
unique_genes_up_df <- data.frame(
  ensembl_gene_id = unique(c(
    exp_upBRCA1_chr17_sthr1.40$source_ensembl_id,
    exp_upBRCA1_chr17_sthr1.40$target_ensembl_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ensembl_gene_id))

# --- Node Construction ---
# Build comprehensive node table with multi-layer annotations:
# 1. Map Ensembl IDs to gene symbols and chromosomes
# 2. Add gene descriptions and Entrez IDs for external database linking
# 3. Annotate expression status (UP_REG/DOWN_REG/UNCHANGED) based on DEG results
nodes_exp_upBRCA1_chr17_sthr1.40 <- unique_genes_up_df %>%
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
  # UP_REG: logFC > 1 (2-fold upregulation)
  # DOWN_REG: logFC < -1 (2-fold downregulation)
  # UNCHANGED: -1 <= logFC <= 1 (not differentially expressed)
  left_join(
    up_down_brca1_deg %>%
      mutate(
        expression_status = case_when(
          logFC > 1  ~ "UP_REG",
          logFC < -1 ~ "DOWN_REG",
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
message("Expression upBRCA1 - Chromosome distribution:")
print(table(nodes_exp_upBRCA1_chr17_sthr1.40$chromosome))
message("Expression upBRCA1 - Status distribution:")
print(table(nodes_exp_upBRCA1_chr17_sthr1.40$status))

# Add layer-specific suffix to prevent ID collisions across omics layers
nodes_exp_upBRCA1_chr17_sthr1.40 <- nodes_exp_upBRCA1_chr17_sthr1.40 %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
# Build regulatory relationship table
rels_exp_upBRCA1_chr17_sthr1.40 <- exp_upBRCA1_chr17_sthr1.40 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "expression"
  ) %>%
  distinct(from, to, .keep_all = TRUE)  # Remove duplicate edges

# Add layer suffix to match node IDs
rels_exp_upBRCA1_chr17_sthr1.40 <- rels_exp_upBRCA1_chr17_sthr1.40 %>%
  mutate(
    from = paste0(from, "_expression"),
    to   = paste0(to, "_expression")
  )

# --- Filter for Biologically Relevant Edges ---
# Rationale: Focus only on edges between genes showing significant expression
# changes (UP_REG or DOWN_REG). This removes noise from unchanged genes and
# highlights the core regulatory network responding to BRCA1 status.
rels_exp_upbrca1_UD_sthr1.40 <- rels_exp_upBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_exp_upBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_upBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  # Both nodes must be differentially expressed
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
# Rationale: "Equivalent" and "Opposite" relationships are symmetric (A-B = B-A).
# Keep only one direction to avoid redundancy in undirected network visualization.
# Other relation types (regulatory, pathway) are directional and kept as-is.
rels_exp_upbrca1_UD_equiv_opp <- rels_exp_upbrca1_UD_sthr1.40 %>%
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
rels_exp_upbrca1_UD_others <- rels_exp_upbrca1_UD_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

# Merge deduplicated edges back
rels_exp_upbrca1_UD_sthr1.40 <- bind_rows(
  rels_exp_upbrca1_UD_others,
  rels_exp_upbrca1_UD_equiv_opp
)

# Summary statistics
message("Expression upBRCA1 UD - Relation type distribution:")
print(table(rels_exp_upbrca1_UD_sthr1.40$relation_type))

# --- Network Integrity Validation ---
# Check for orphaned nodes (nodes without any edges)
orphaned_nodes <- setdiff(
  nodes_exp_upBRCA1_chr17_sthr1.40$id,
  c(rels_exp_upbrca1_UD_sthr1.40$from, rels_exp_upbrca1_UD_sthr1.40$to)
)
message("Expression upBRCA1 - Orphaned nodes: ", length(orphaned_nodes))

# Check for edges referencing non-existent nodes
missing_from_nodes <- setdiff(
  rels_exp_upbrca1_UD_sthr1.40$from,
  nodes_exp_upBRCA1_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_exp_upbrca1_UD_sthr1.40$to,
  nodes_exp_upBRCA1_chr17_sthr1.40$id
)
message("Expression upBRCA1 - Missing FROM nodes: ", length(missing_from_nodes))
message("Expression upBRCA1 - Missing TO nodes: ", length(missing_to_nodes))


# ==============================================================================
# PART 2: EXPRESSION NETWORK - downBRCA1 (Chr17 Focus)
# ==============================================================================

# --- Load Expression Relations Data ---
exp_downBRCA1_chr17_sthr1.40 <- read.table(
  "data/expression/downBRCA1/all_relations_downexp_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Genes ---
unique_genes_down_df <- data.frame(
  ensembl_gene_id = unique(c(
    exp_downBRCA1_chr17_sthr1.40$source_ensembl_id,
    exp_downBRCA1_chr17_sthr1.40$target_ensembl_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(ensembl_gene_id))

# --- Node Construction ---
nodes_exp_downBRCA1_chr17_sthr1.40 <- unique_genes_down_df %>%
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
  left_join(
    up_down_brca1_deg %>%
      mutate(
        expression_status = case_when(
          logFC > 1  ~ "UP_REG",
          logFC < -1 ~ "DOWN_REG",
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

message("Expression downBRCA1 - Chromosome distribution:")
print(table(nodes_exp_downBRCA1_chr17_sthr1.40$chromosome))
message("Expression downBRCA1 - Status distribution:")
print(table(nodes_exp_downBRCA1_chr17_sthr1.40$status))

nodes_exp_downBRCA1_chr17_sthr1.40 <- nodes_exp_downBRCA1_chr17_sthr1.40 %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
rels_exp_downBRCA1_chr17_sthr1.40 <- exp_downBRCA1_chr17_sthr1.40 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "expression"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_exp_downBRCA1_chr17_sthr1.40 <- rels_exp_downBRCA1_chr17_sthr1.40 %>%
  mutate(
    from = paste0(from, "_expression"),
    to   = paste0(to, "_expression")
  )

# --- Filter for Biologically Relevant Edges ---
rels_exp_downBRCA1_UD_sthr1.40 <- rels_exp_downBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_exp_downBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_downBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_exp_UD_down_equiv_opp <- rels_exp_downBRCA1_UD_sthr1.40 %>%
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

rels_exp_UD_down_others <- rels_exp_downBRCA1_UD_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_exp_downBRCA1_UD_sthr1.40 <- bind_rows(
  rels_exp_UD_down_others,
  rels_exp_UD_down_equiv_opp
)

message("Expression downBRCA1 UD - Relation type distribution:")
print(table(rels_exp_downBRCA1_UD_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_exp_downBRCA1_chr17_sthr1.40$id,
  c(rels_exp_downBRCA1_UD_sthr1.40$from, rels_exp_downBRCA1_UD_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_exp_downBRCA1_UD_sthr1.40$from,
  nodes_exp_downBRCA1_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_exp_downBRCA1_UD_sthr1.40$to,
  nodes_exp_downBRCA1_chr17_sthr1.40$id
)
message("Expression downBRCA1 - Orphaned nodes: ", length(orphaned_nodes))
message("Expression downBRCA1 - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 3: METHYLATION NETWORK - upBRCA1 (Chr17 Focus)
# ==============================================================================
# Rationale: Methylation patterns often correlate with gene expression changes.
# Focus on CpG sites showing significant differential methylation (|logFC| > 0.2)
# as these are more likely to have functional impact on gene regulation.

# --- Load Methylation Relations Data ---
met_upBRCA1_chr17_sthr1.40 <- read.table(
  "data/methylation/upBRCA1/all_relations_upmet_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Probes ---
# Methylation uses probe IDs rather than gene IDs directly
unique_probes_up_df <- data.frame(
  probeID = unique(c(
    met_upBRCA1_chr17_sthr1.40$source_probe_id,
    met_upBRCA1_chr17_sthr1.40$target_probe_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(probeID))

# --- Node Construction ---
# Multi-step annotation process:
# 1. Map probe IDs to gene names and genomic locations
# 2. Add methylation status (HYPERMETH/HYPOMETH) from DMP results
# 3. Add gene descriptions and Entrez IDs
# CpG location context (TSS200, 5'UTR, Body, etc.) is important for
# understanding regulatory mechanism
nodes_met_upBRCA1_chr17_sthr1.40 <- unique_probes_up_df %>%
  # Layer 1: Probe annotations (probe -> gene + location)
  left_join(probeInfo, by = "probeID") %>%
  
  # Layer 2: Methylation status classification
  # HYPERMETH: logFC > 0.2 (increased methylation)
  # HYPOMETH: logFC < -0.2 (decreased methylation)
  # Threshold of 0.2 chosen based on biological significance
  left_join(
    upbrca1_vs_downbrca1_dmp %>%
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

message("Methylation upBRCA1 - Chromosome distribution:")
print(table(nodes_met_upBRCA1_chr17_sthr1.40$chromosome))
message("Methylation upBRCA1 - Status distribution:")
print(table(nodes_met_upBRCA1_chr17_sthr1.40$status))

nodes_met_upBRCA1_chr17_sthr1.40 <- nodes_met_upBRCA1_chr17_sthr1.40 %>%
  mutate(id = paste0(id, "_methylation"))

# --- Edge Construction ---
rels_met_upBRCA1_chr17_sthr1.40 <- met_upBRCA1_chr17_sthr1.40 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "methylation"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_met_upBRCA1_chr17_sthr1.40 <- rels_met_upBRCA1_chr17_sthr1.40 %>%
  mutate(
    from = paste0(from, "_methylation"),
    to   = paste0(to, "_methylation")
  )

# --- Filter for Biologically Relevant Edges ---
rels_met_upbrca1_UD_sthr1.40 <- rels_met_upBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_met_upBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_upBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_up_UD_equiv_opp <- rels_met_upbrca1_UD_sthr1.40 %>%
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

rels_met_up_UD_others <- rels_met_upbrca1_UD_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_upbrca1_UD_sthr1.40 <- bind_rows(
  rels_met_up_UD_others,
  rels_met_up_UD_equiv_opp
)

message("Methylation upBRCA1 UD - Relation type distribution:")
print(table(rels_met_upbrca1_UD_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_met_upBRCA1_chr17_sthr1.40$id,
  c(rels_met_upbrca1_UD_sthr1.40$from, rels_met_upbrca1_UD_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_met_upbrca1_UD_sthr1.40$from,
  nodes_met_upBRCA1_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_met_upbrca1_UD_sthr1.40$to,
  nodes_met_upBRCA1_chr17_sthr1.40$id
)
message("Methylation upBRCA1 - Orphaned nodes: ", length(orphaned_nodes))
message("Methylation upBRCA1 - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 4: METHYLATION NETWORK - downBRCA1 (Chr17 Focus)
# ==============================================================================

# --- Load Methylation Relations Data ---
met_downBRCA1_chr17_sthr1.40 <- read.table(
  "data/methylation/downBRCA1/all_relations_downmet_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)

# --- Extract Unique Probes ---
unique_probes_down_df <- data.frame(
  probeID = unique(c(
    met_downBRCA1_chr17_sthr1.40$source_probe_id,
    met_downBRCA1_chr17_sthr1.40$target_probe_id
  )),
  stringsAsFactors = FALSE
) %>%
  filter(!is.na(probeID))

# --- Node Construction ---
nodes_met_downBRCA1_chr17_sthr1.40 <- unique_probes_down_df %>%
  left_join(probeInfo, by = "probeID") %>%
  left_join(
    upbrca1_vs_downbrca1_dmp %>%
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

message("Methylation downBRCA1 - Chromosome distribution:")
print(table(nodes_met_downBRCA1_chr17_sthr1.40$chromosome))
message("Methylation downBRCA1 - Status distribution:")
print(table(nodes_met_downBRCA1_chr17_sthr1.40$status))

nodes_met_downBRCA1_chr17_sthr1.40 <- nodes_met_downBRCA1_chr17_sthr1.40 %>%
  mutate(id = paste0(id, "_methylation"))

# --- Edge Construction ---
rels_met_downBRCA1_chr17_sthr1.40 <- met_downBRCA1_chr17_sthr1.40 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "methylation"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_met_downBRCA1_chr17_sthr1.40 <- rels_met_downBRCA1_chr17_sthr1.40 %>%
  mutate(
    from = paste0(from, "_methylation"),
    to   = paste0(to, "_methylation")
  )

# --- Filter for Biologically Relevant Edges ---
rels_met_downbrca1_UD_sthr1.40 <- rels_met_downBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_met_downBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_downBRCA1_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_UD_down_equiv_opp <- rels_met_downbrca1_UD_sthr1.40 %>%
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

rels_met_UD_down_others <- rels_met_downbrca1_UD_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_downbrca1_UD_sthr1.40 <- bind_rows(
  rels_met_UD_down_others,
  rels_met_UD_down_equiv_opp
)

message("Methylation downBRCA1 UD - Relation type distribution:")
print(table(rels_met_downbrca1_UD_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_met_downBRCA1_chr17_sthr1.40$id,
  c(rels_met_downbrca1_UD_sthr1.40$from, rels_met_downbrca1_UD_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_met_downbrca1_UD_sthr1.40$from,
  nodes_met_downBRCA1_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_met_downbrca1_UD_sthr1.40$to,
  nodes_met_downBRCA1_chr17_sthr1.40$id
)
message("Methylation downBRCA1 - Orphaned nodes: ", length(orphaned_nodes))
message("Methylation downBRCA1 - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
# Output objects available in workspace:
# - nodes_exp_upBRCA1_chr17_sthr1.40, rels_exp_upbrca1_UD_sthr1.40
# - nodes_exp_downBRCA1_chr17_sthr1.40, rels_exp_downBRCA1_UD_sthr1.40
# - nodes_met_upBRCA1_chr17_sthr1.40, rels_met_upbrca1_UD_sthr1.40
# - nodes_met_downBRCA1_chr17_sthr1.40, rels_met_downbrca1_UD_sthr1.40
# ==============================================================================
