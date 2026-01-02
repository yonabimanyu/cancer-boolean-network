# ==============================================================================
# BIOLOGICAL NETWORK ANALYSIS: BRCA1 vs CONTROL CONTRASTS
# ==============================================================================
# Title: Network Construction for BRCA1 Conditions vs Control Comparison
# Description: Constructs gene regulatory networks comparing upBRCA1 vs Control
#              (UC) and downBRCA1 vs Control (DC). Focuses on chromosome 17 genes
#              showing significant changes relative to control samples. This
#              analysis complements the UD (upBRCA1 vs downBRCA1) contrast by
#              establishing baseline differences from normal tissue.
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
#   - data/reference/gene_info_final_3.5.rds
#   - data/reference/gene_info_final_4.rds
#   - data/reference/probeInfo.rds
#   - data/deg/up_ctrl_brca1_deg.rds (upBRCA1 vs Control DEG)
#   - data/deg/down_ctrl_brca1_deg.rds (downBRCA1 vs Control DEG)
#   - data/dmp/upbrca1_vs_ctrl_dmp.rds (upBRCA1 vs Control DMP)
#   - data/dmp/downbrca1_vs_ctrl_dmp.rds (downBRCA1 vs Control DMP)
#
# Prerequisites: Run the UD contrast script first to generate unique gene/probe
#                lists (unique_genes_up_df, unique_genes_down_df,
#                unique_probes_up_df, unique_probes_down_df)
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
gene_info_final_3.5 <- readRDS("data/reference/gene_info_final_3.5.rds")
gene_info_final_4 <- readRDS("data/reference/gene_info_final_4.rds")
probeInfo <- readRDS("data/reference/probeInfo.rds")

# Control contrast differential analysis results
up_ctrl_brca1_deg <- readRDS("data/deg/up_ctrl_brca1_deg.rds")
down_ctrl_brca1_deg <- readRDS("data/deg/down_ctrl_brca1_deg.rds")
upbrca1_vs_ctrl_dmp <- readRDS("data/dmp/upbrca1_vs_ctrl_dmp.rds")
downbrca1_vs_ctrl_dmp <- readRDS("data/dmp/downbrca1_vs_ctrl_dmp.rds")

# === LOAD UNIQUE GENE/PROBE LISTS ===
# These should be generated from the UD contrast script
# unique_genes_up_df <- readRDS("data/processed/unique_genes_up.rds")
# unique_genes_down_df <- readRDS("data/processed/unique_genes_down.rds")
# unique_probes_up_df <- readRDS("data/processed/unique_probes_up.rds")
# unique_probes_down_df <- readRDS("data/processed/unique_probes_down.rds")

# Load relation tables (same as UD contrast, will be filtered differently)
exp_upBRCA1_chr17_sthr1.40 <- read.table(
  "data/expression/upBRCA1/all_relations_upexp_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
exp_downBRCA1_chr17_sthr1.40 <- read.table(
  "data/expression/downBRCA1/all_relations_downexp_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
met_upBRCA1_chr17_sthr1.40 <- read.table(
  "data/methylation/upBRCA1/all_relations_upmet_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
met_downBRCA1_chr17_sthr1.40 <- read.table(
  "data/methylation/downBRCA1/all_relations_downmet_1.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)


# ==============================================================================
# PART 1: EXPRESSION NETWORK - upBRCA1 vs Control (UC)
# ==============================================================================
# Rationale: Compare upBRCA1 samples against control to identify genes that are
# specifically dysregulated when BRCA1 is upregulated. This establishes which
# changes are BRCA1-dependent vs general tumor characteristics.

# --- Node Construction ---
# Same gene list as UD contrast, but annotated with UC differential expression
nodes_exp_UC_chr17_sthr1.40 <- unique_genes_up_df %>%
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
  # Key difference: Use upBRCA1 vs Control DEG results instead of UD
  left_join(
    up_ctrl_brca1_deg %>%
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

message("Expression UC - Chromosome distribution:")
print(table(nodes_exp_UC_chr17_sthr1.40$chromosome))
message("Expression UC - Status distribution:")
print(table(nodes_exp_UC_chr17_sthr1.40$status))

nodes_exp_UC_chr17_sthr1.40 <- nodes_exp_UC_chr17_sthr1.40 %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
# Use same relations as UD, will be filtered with UC status
rels_exp_upBRCA1_chr17_sthr1.40 <- exp_upBRCA1_chr17_sthr1.40 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "expression"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_exp_upBRCA1_chr17_sthr1.40 <- rels_exp_upBRCA1_chr17_sthr1.40 %>%
  mutate(
    from = paste0(from, "_expression"),
    to   = paste0(to, "_expression")
  )

# --- Filter for Biologically Relevant Edges (UC Status) ---
# Difference from UD: Different genes will be classified as UP_REG/DOWN_REG
# when compared to control vs when compared to downBRCA1
rels_exp_upbrca1_UC_sthr1.40 <- rels_exp_upBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_exp_UC_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_UC_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_exp_UC_equiv_opp <- rels_exp_upbrca1_UC_sthr1.40 %>%
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

rels_exp_UC_others <- rels_exp_upbrca1_UC_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_exp_upbrca1_UC_sthr1.40 <- bind_rows(
  rels_exp_UC_others,
  rels_exp_UC_equiv_opp
)

message("Expression UC - Relation type distribution:")
print(table(rels_exp_upbrca1_UC_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_exp_UC_chr17_sthr1.40$id,
  c(rels_exp_upbrca1_UC_sthr1.40$from, rels_exp_upbrca1_UC_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_exp_upbrca1_UC_sthr1.40$from,
  nodes_exp_UC_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_exp_upbrca1_UC_sthr1.40$to,
  nodes_exp_UC_chr17_sthr1.40$id
)
message("Expression UC - Orphaned nodes: ", length(orphaned_nodes))
message("Expression UC - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 2: EXPRESSION NETWORK - downBRCA1 vs Control (DC)
# ==============================================================================

# --- Node Construction ---
nodes_exp_DC_chr17_sthr1.40 <- unique_genes_down_df %>%
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
  # Use downBRCA1 vs Control DEG results
  left_join(
    down_ctrl_brca1_deg %>%
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

message("Expression DC - Chromosome distribution:")
print(table(nodes_exp_DC_chr17_sthr1.40$chromosome))
message("Expression DC - Status distribution:")
print(table(nodes_exp_DC_chr17_sthr1.40$status))

nodes_exp_DC_chr17_sthr1.40 <- nodes_exp_DC_chr17_sthr1.40 %>%
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

# --- Filter for Biologically Relevant Edges (DC Status) ---
rels_exp_downbrca1_DC_sthr1.40 <- rels_exp_downBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_exp_DC_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_DC_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_exp_down_DC_equiv_opp <- rels_exp_downbrca1_DC_sthr1.40 %>%
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

rels_exp_down_DC_others <- rels_exp_downbrca1_DC_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_exp_downbrca1_DC_sthr1.40 <- bind_rows(
  rels_exp_down_DC_others,
  rels_exp_down_DC_equiv_opp
)

message("Expression DC - Relation type distribution:")
print(table(rels_exp_downbrca1_DC_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_exp_DC_chr17_sthr1.40$id,
  c(rels_exp_downbrca1_DC_sthr1.40$from, rels_exp_downbrca1_DC_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_exp_downbrca1_DC_sthr1.40$from,
  nodes_exp_DC_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_exp_downbrca1_DC_sthr1.40$to,
  nodes_exp_DC_chr17_sthr1.40$id
)
message("Expression DC - Orphaned nodes: ", length(orphaned_nodes))
message("Expression DC - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 3: METHYLATION NETWORK - upBRCA1 vs Control (UC)
# ==============================================================================
# Rationale: Identify CpG sites showing differential methylation specifically
# when BRCA1 is upregulated compared to control. These may represent epigenetic
# regulatory mechanisms driving the upBRCA1 phenotype.

# --- Node Construction ---
nodes_met_UC_chr17_sthr1.40 <- unique_probes_up_df %>%
  left_join(probeInfo, by = "probeID") %>%
  # Use upBRCA1 vs Control DMP results
  left_join(
    upbrca1_vs_ctrl_dmp %>%
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

message("Methylation UC - Chromosome distribution:")
print(table(nodes_met_UC_chr17_sthr1.40$chromosome))
message("Methylation UC - Status distribution:")
print(table(nodes_met_UC_chr17_sthr1.40$status))

nodes_met_UC_chr17_sthr1.40 <- nodes_met_UC_chr17_sthr1.40 %>%
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

# --- Filter for Biologically Relevant Edges (UC Status) ---
# Note: Much higher edge count than expression due to larger number of
# differentially methylated sites in UC contrast
rels_met_upbrca1_UC_sthr1.40 <- rels_met_upBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_met_UC_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_UC_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_up_UC_equiv_opp <- rels_met_upbrca1_UC_sthr1.40 %>%
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

rels_met_up_UC_others <- rels_met_upbrca1_UC_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_upbrca1_UC_sthr1.40 <- bind_rows(
  rels_met_up_UC_others,
  rels_met_up_UC_equiv_opp
)

message("Methylation UC - Relation type distribution:")
print(table(rels_met_upbrca1_UC_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_met_UC_chr17_sthr1.40$id,
  c(rels_met_upbrca1_UC_sthr1.40$from, rels_met_upbrca1_UC_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_met_upbrca1_UC_sthr1.40$from,
  nodes_met_UC_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_met_upbrca1_UC_sthr1.40$to,
  nodes_met_UC_chr17_sthr1.40$id
)
message("Methylation UC - Orphaned nodes: ", length(orphaned_nodes))
message("Methylation UC - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 4: METHYLATION NETWORK - downBRCA1 vs Control (DC)
# ==============================================================================

# --- Node Construction ---
nodes_met_DC_chr17_sthr1.40 <- unique_probes_down_df %>%
  left_join(probeInfo, by = "probeID") %>%
  # Use downBRCA1 vs Control DMP results
  left_join(
    downbrca1_vs_ctrl_dmp %>%
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

message("Methylation DC - Chromosome distribution:")
print(table(nodes_met_DC_chr17_sthr1.40$chromosome))
message("Methylation DC - Status distribution:")
print(table(nodes_met_DC_chr17_sthr1.40$status))

nodes_met_DC_chr17_sthr1.40 <- nodes_met_DC_chr17_sthr1.40 %>%
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

# --- Filter for Biologically Relevant Edges (DC Status) ---
rels_met_downbrca1_DC_sthr1.40 <- rels_met_downBRCA1_chr17_sthr1.40 %>%
  left_join(
    nodes_met_DC_chr17_sthr1.40 %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_DC_chr17_sthr1.40 %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_DC_down_equiv_opp <- rels_met_downbrca1_DC_sthr1.40 %>%
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

rels_met_DC_down_others <- rels_met_downbrca1_DC_sthr1.40 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_downbrca1_DC_sthr1.40 <- bind_rows(
  rels_met_DC_down_others,
  rels_met_DC_down_equiv_opp
)

message("Methylation DC - Relation type distribution:")
print(table(rels_met_downbrca1_DC_sthr1.40$relation_type))

# --- Network Integrity Validation ---
orphaned_nodes <- setdiff(
  nodes_met_DC_chr17_sthr1.40$id,
  c(rels_met_downbrca1_DC_sthr1.40$from, rels_met_downbrca1_DC_sthr1.40$to)
)
missing_from_nodes <- setdiff(
  rels_met_downbrca1_DC_sthr1.40$from,
  nodes_met_DC_chr17_sthr1.40$id
)
missing_to_nodes <- setdiff(
  rels_met_downbrca1_DC_sthr1.40$to,
  nodes_met_DC_chr17_sthr1.40$id
)
message("Methylation DC - Orphaned nodes: ", length(orphaned_nodes))
message("Methylation DC - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
# Output objects available in workspace:
# Expression networks:
# - nodes_exp_UC_chr17_sthr1.40, rels_exp_upbrca1_UC_sthr1.40 (upBRCA1 vs Ctrl)
# - nodes_exp_DC_chr17_sthr1.40, rels_exp_downbrca1_DC_sthr1.40 (downBRCA1 vs Ctrl)
#
# Methylation networks:
# - nodes_met_UC_chr17_sthr1.40, rels_met_upbrca1_UC_sthr1.40 (upBRCA1 vs Ctrl)
# - nodes_met_DC_chr17_sthr1.40, rels_met_downbrca1_DC_sthr1.40 (downBRCA1 vs Ctrl)
#
# Analysis Strategy:
# Compare UC/DC contrasts with UD contrast to identify:
# 1. BRCA1-specific effects (present in UC/DC but not UD)
# 2. Baseline tumor effects (common across all contrasts)
# 3. Dose-dependent effects (gradient across conditions)
# ==============================================================================
