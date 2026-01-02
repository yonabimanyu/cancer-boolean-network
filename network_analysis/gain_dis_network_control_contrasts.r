# ==============================================================================
# BIOLOGICAL NETWORK ANALYSIS: GAIN/DIS vs CONTROL CONTRASTS
# ==============================================================================
# Title: Network Construction for GAIN and DIS Conditions vs Control Comparison
# Description: Constructs gene regulatory networks comparing GAIN vs Control (GC)
#              and DIS vs Control (DC). Focuses on chromosome 17 genes showing
#              significant changes relative to control samples. This analysis
#              complements the GD (GAIN vs DIS) contrast by establishing baseline
#              differences from normal tissue, helping distinguish:
#              1. Copy number-specific effects (present in GC but not DC)
#              2. General tumor effects (common to both GC and DC)
#              3. Dosage-dependent effects (gradient across conditions)
#
# Author: [Your Name]
# Date: 2025-01-01
# Version: 3.0 (Updated 2025-11-19)
#
# Biological Context:
#   These control contrasts help separate the effects of chromosome 17 copy
#   number gain from general tumor biology. GC contrast identifies changes
#   specifically associated with amplification, while DC contrast reveals
#   changes present even in diploid tumors.
#
# Input Files:
#   - data/expression/GAIN/all_relations_gainexp_2.3.txt
#   - data/expression/DIS/all_relations_disexp_2.3.txt
#   - data/methylation/GAIN/all_relations_gainmet_2.4.txt
#   - data/methylation/DIS/all_relations_dismet_2.4.txt
#   - data/reference/gene_info_final_3.5.rds
#   - data/reference/gene_info_final_4.rds
#   - data/reference/probeInfo.rds
#   - data/deg/gain_vs_ctrl_deg.rds (GAIN vs Control DEG)
#   - data/deg/dis_vs_ctrl_deg.rds (DIS vs Control DEG)
#   - data/dmp/gain_vs_ctrl_dmp.rds (GAIN vs Control DMP)
#   - data/dmp/dis_vs_ctrl_dmp.rds (DIS vs Control DMP)
#
# Prerequisites: Run the GD contrast script first to generate unique gene/probe
#                lists (unique_genes_gain_df, unique_genes_dis_df,
#                unique_probes_met_df, unique_probes_met_dis_df)
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
gain_vs_ctrl_deg <- readRDS("data/deg/gain_vs_ctrl_deg.rds")
dis_vs_ctrl_deg <- readRDS("data/deg/dis_vs_ctrl_deg.rds")
gain_vs_ctrl_dmp <- readRDS("data/dmp/gain_vs_ctrl_dmp.rds")
dis_vs_ctrl_dmp <- readRDS("data/dmp/dis_vs_ctrl_dmp.rds")

# === LOAD UNIQUE GENE/PROBE LISTS ===
# These should be generated from the GD contrast script
# If not in workspace, uncomment to load:
# unique_genes_gain_df <- readRDS("data/processed/unique_genes_gain.rds")
# unique_genes_dis_df <- readRDS("data/processed/unique_genes_dis.rds")
# unique_probes_met_df <- readRDS("data/processed/unique_probes_gain.rds")
# unique_probes_met_dis_df <- readRDS("data/processed/unique_probes_dis.rds")

# Load relation tables (same as GD contrast, will be filtered differently)
exp_gain_chr17_sthr2.3 <- read.table(
  "data/expression/GAIN/all_relations_gainexp_2.3.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
exp_dis_chr17_sthr2.3 <- read.table(
  "data/expression/DIS/all_relations_disexp_2.3.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
met_gain_chr17_sthr2.4 <- read.table(
  "data/methylation/GAIN/all_relations_gainmet_2.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)
met_dis_chr17_sthr2.4 <- read.table(
  "data/methylation/DIS/all_relations_dismet_2.4.txt",
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = ""
)


# ==============================================================================
# PART 1: EXPRESSION NETWORK - GAIN vs Control (GC)
# ==============================================================================
# Rationale: Compare GAIN samples against control to identify genes that are
# specifically dysregulated when chromosome 17 is amplified. This establishes
# which changes are copy number-dependent vs general tumor characteristics.

# --- Node Construction ---
# Same gene list as GD contrast, but annotated with GC differential expression
# Expected count: 868 nodes (as of 2025-11-19)
nodes_exp_gain_chr17_sthr2.3_GC <- unique_genes_gain_df %>%
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
  # Key difference: Use GAIN vs Control DEG results instead of GD
  left_join(
    gain_vs_ctrl_deg %>%
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

message("Expression GAIN vs Control - Chromosome distribution:")
print(table(nodes_exp_gain_chr17_sthr2.3_GC$chromosome))
message("Expression GAIN vs Control - Status distribution:")
print(table(nodes_exp_gain_chr17_sthr2.3_GC$status))

nodes_exp_gain_chr17_sthr2.3_GC <- nodes_exp_gain_chr17_sthr2.3_GC %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
# Use same relations as GD, will be filtered with GC status
# Expected count before filtering: 18,258 edges
rels_exp_gain_chr17_sthr2.3 <- exp_gain_chr17_sthr2.3 %>%
  transmute(
    from          = source_gene,
    to            = target_gene,
    relation_type = relation_type,
    layer         = "expression"
  ) %>%
  distinct(from, to, .keep_all = TRUE)

rels_exp_gain_chr17_sthr2.3 <- rels_exp_gain_chr17_sthr2.3 %>%
  mutate(
    from = paste0(from, "_expression"),
    to   = paste0(to, "_expression")
  )

# --- Filter for Biologically Relevant Edges (GC Status) ---
# Expected count: 1,539 edges (as of 2025-11-19)
# Difference from GD: Different genes will be classified as OVER_T/DOWN_T
# when compared to control vs when compared to DIS
rels_exp_gain_GC_sthr2.3 <- rels_exp_gain_chr17_sthr2.3 %>%
  left_join(
    nodes_exp_gain_chr17_sthr2.3_GC %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_gain_chr17_sthr2.3_GC %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
# Note: Variable name inconsistency in original (used rels_exp_gain_GD_equiv_opp)
# Correcting to rels_exp_gain_GC_equiv_opp for clarity
rels_exp_gain_GC_equiv_opp <- rels_exp_gain_GC_sthr2.3 %>%
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

rels_exp_gain_GC_others <- rels_exp_gain_GC_sthr2.3 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_exp_gain_GC_sthr2.3 <- bind_rows(
  rels_exp_gain_GC_others,
  rels_exp_gain_GC_equiv_opp
)

message("Expression GAIN vs Control - Relation type distribution:")
print(table(rels_exp_gain_GC_sthr2.3$relation_type))

# --- Network Integrity Validation ---
missing_from_nodes <- setdiff(
  rels_exp_gain_GC_sthr2.3$from,
  nodes_exp_gain_chr17_sthr2.3_GC$id
)
missing_to_nodes <- setdiff(
  rels_exp_gain_GC_sthr2.3$to,
  nodes_exp_gain_chr17_sthr2.3_GC$id
)
message("Expression GAIN vs Control - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 2: EXPRESSION NETWORK - DIS vs Control (DC)
# ==============================================================================
# Rationale: Identify genes dysregulated in diploid tumor samples compared to
# control. This reveals baseline tumor effects independent of chr17 amplification.

# --- Node Construction ---
# Expected count: 1,002 nodes (as of 2025-11-19)
nodes_exp_dis_chr17_sthr2.3_DC <- unique_genes_dis_df %>%
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
  # Use DIS vs Control DEG results
  left_join(
    dis_vs_ctrl_deg %>%
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

message("Expression DIS vs Control - Chromosome distribution:")
print(table(nodes_exp_dis_chr17_sthr2.3_DC$chromosome))
message("Expression DIS vs Control - Status distribution:")
print(table(nodes_exp_dis_chr17_sthr2.3_DC$status))

nodes_exp_dis_chr17_sthr2.3_DC <- nodes_exp_dis_chr17_sthr2.3_DC %>%
  mutate(id = paste0(id, "_expression"))

# --- Edge Construction ---
# Expected count before filtering: 101,640 edges
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

# --- Filter for Biologically Relevant Edges (DC Status) ---
# Expected count: 2,026 edges (as of 2025-11-19)
rels_exp_dis_DC_sthr2.3 <- rels_exp_dis_chr17_sthr2.3 %>%
  left_join(
    nodes_exp_dis_chr17_sthr2.3_DC %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_exp_dis_chr17_sthr2.3_DC %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_exp_dis_DC_equiv_opp <- rels_exp_dis_DC_sthr2.3 %>%
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

rels_exp_dis_DC_others <- rels_exp_dis_DC_sthr2.3 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_exp_dis_DC_sthr2.3 <- bind_rows(
  rels_exp_dis_DC_others,
  rels_exp_dis_DC_equiv_opp
)

message("Expression DIS vs Control - Relation type distribution:")
print(table(rels_exp_dis_DC_sthr2.3$relation_type))

# --- Network Integrity Validation ---
missing_from_nodes <- setdiff(
  rels_exp_dis_DC_sthr2.3$from,
  nodes_exp_dis_chr17_sthr2.3_DC$id
)
missing_to_nodes <- setdiff(
  rels_exp_dis_DC_sthr2.3$to,
  nodes_exp_dis_chr17_sthr2.3_DC$id
)
message("Expression DIS vs Control - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 3: METHYLATION NETWORK - GAIN vs Control (GC)
# ==============================================================================
# Rationale: Identify CpG sites showing differential methylation specifically
# when chromosome 17 is amplified compared to control. These may represent
# epigenetic regulatory mechanisms driving the GAIN phenotype or responding
# to increased gene dosage.

# --- Node Construction ---
# Expected count: 1,496 nodes (as of 2025-11-19)
nodes_met_gain_chr17_sthr2.4_GC <- unique_probes_met_df %>%
  left_join(probeInfo, by = "probeID") %>%
  # Use GAIN vs Control DMP results
  left_join(
    gain_vs_ctrl_dmp %>%
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

message("Methylation GAIN vs Control - Chromosome distribution:")
print(table(nodes_met_gain_chr17_sthr2.4_GC$chromosome))
message("Methylation GAIN vs Control - Status distribution:")
print(table(nodes_met_gain_chr17_sthr2.4_GC$status))

nodes_met_gain_chr17_sthr2.4_GC <- nodes_met_gain_chr17_sthr2.4_GC %>%
  mutate(id = paste0(id, "_methylation"))

# --- Edge Construction ---
# Expected count before filtering: 23,170 edges
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

# --- Filter for Biologically Relevant Edges (GC Status) ---
# Expected count: 15,756 edges (as of 2025-11-19)
# Much higher edge count than expression, indicating widespread methylation
# changes associated with chr17 amplification
rels_met_gain_GC_sthr2.4 <- rels_met_gain_chr17_sthr2.4 %>%
  left_join(
    nodes_met_gain_chr17_sthr2.4_GC %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_gain_chr17_sthr2.4_GC %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_gain_GC_equiv_opp <- rels_met_gain_GC_sthr2.4 %>%
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

rels_met_gain_GC_others <- rels_met_gain_GC_sthr2.4 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_gain_GC_sthr2.4 <- bind_rows(
  rels_met_gain_GC_others,
  rels_met_gain_GC_equiv_opp
)

message("Methylation GAIN vs Control - Relation type distribution:")
print(table(rels_met_gain_GC_sthr2.4$relation_type))

# --- Network Integrity Validation ---
missing_from_nodes <- setdiff(
  rels_met_gain_GC_sthr2.4$from,
  nodes_met_gain_chr17_sthr2.4_GC$id
)
missing_to_nodes <- setdiff(
  rels_met_gain_GC_sthr2.4$to,
  nodes_met_gain_chr17_sthr2.4_GC$id
)
message("Methylation GAIN vs Control - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))


# ==============================================================================
# PART 4: METHYLATION NETWORK - DIS vs Control (DC)
# ==============================================================================

# --- Node Construction ---
# Expected count: 1,659 nodes (as of 2025-11-19)
nodes_met_dis_chr17_sthr2.4_DC <- unique_probes_met_dis_df %>%
  left_join(probeInfo, by = "probeID") %>%
  # Use DIS vs Control DMP results
  left_join(
    dis_vs_ctrl_dmp %>%
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

message("Methylation DIS vs Control - Chromosome distribution:")
print(table(nodes_met_dis_chr17_sthr2.4_DC$chromosome))
message("Methylation DIS vs Control - Status distribution:")
print(table(nodes_met_dis_chr17_sthr2.4_DC$status))

nodes_met_dis_chr17_sthr2.4_DC <- nodes_met_dis_chr17_sthr2.4_DC %>%
  mutate(id = paste0(id, "_methylation"))

# --- Edge Construction ---
# Expected count before filtering: 65,632 edges
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

# --- Filter for Biologically Relevant Edges (DC Status) ---
# Expected count: 44,413 edges (as of 2025-11-19)
# Notably higher than GAIN GC (15,756), suggesting DIS samples show more
# widespread methylation changes vs control
rels_met_dis_DC_sthr2.4 <- rels_met_dis_chr17_sthr2.4 %>%
  left_join(
    nodes_met_dis_chr17_sthr2.4_DC %>% select(id, status),
    by = c("from" = "id")
  ) %>%
  rename(from_status = status) %>%
  left_join(
    nodes_met_dis_chr17_sthr2.4_DC %>% select(id, status),
    by = c("to" = "id")
  ) %>%
  rename(to_status = status) %>%
  filter(from_status != "UNCHANGED" & to_status != "UNCHANGED") %>%
  select(-from_status, -to_status)

# --- Deduplicate Symmetric Relations ---
rels_met_dis_DC_equiv_opp <- rels_met_dis_DC_sthr2.4 %>%
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

rels_met_dis_DC_others <- rels_met_dis_DC_sthr2.4 %>%
  filter(!relation_type %in% c("Equivalent", "Opposite"))

rels_met_dis_DC_sthr2.4 <- bind_rows(
  rels_met_dis_DC_others,
  rels_met_dis_DC_equiv_opp
)

message("Methylation DIS vs Control - Relation type distribution:")
print(table(rels_met_dis_DC_sthr2.4$relation_type))

# --- Network Integrity Validation ---
missing_from_nodes <- setdiff(
  rels_met_dis_DC_sthr2.4$from,
  nodes_met_dis_chr17_sthr2.4_DC$id
)
missing_to_nodes <- setdiff(
  rels_met_dis_DC_sthr2.4$to,
  nodes_met_dis_chr17_sthr2.4_DC$id
)
message("Methylation DIS vs Control - Missing FROM/TO nodes: ",
        length(missing_from_nodes), "/", length(missing_to_nodes))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
# Output objects available in workspace:
# Expression networks:
# - nodes_exp_gain_chr17_sthr2.3_GC, rels_exp_gain_GC_sthr2.3 (GAIN vs Ctrl)
# - nodes_exp_dis_chr17_sthr2.3_DC, rels_exp_dis_DC_sthr2.3 (DIS vs Ctrl)
#
# Methylation networks:
# - nodes_met_gain_chr17_sthr2.4_GC, rels_met_gain_GC_sthr2.4 (GAIN vs Ctrl)
# - nodes_met_dis_chr17_sthr2.4_DC, rels_met_dis_DC_sthr2.4 (DIS vs Ctrl)
#
# Analysis Strategy:
# Compare GC/DC contrasts with GD contrast to identify:
# 1. Copy number-specific effects (present in GC, absent in DC)
#    - Genes with GAIN-specific expression changes
#    - CpG sites with amplification-associated methylation
# 2. Baseline tumor effects (common to both GC and DC)
#    - General tumor biology independent of chr17 status
# 3. Dosage-dependent effects (gradient: Control < DIS < GAIN)
#    - Genes showing progressive changes with copy number
#
# Key Findings:
# - Expression networks: GC (1,539 edges) > DC (2,026 edges)
# - Methylation networks: DC (44,413 edges) >> GC (15,756 edges)
#   Suggests diploid tumors have more methylation dysregulation vs control
#   compared to chr17 amplification effects
# ==============================================================================