# ==============================================================================
# MULTI-OMICS NETWORK INTEGRATION: CHROMOSOME 17 COPY NUMBER VARIATION
# ==============================================================================
# Title: Interlayer Network Construction for GAIN/DIS Groups
# Description: Integrates expression and methylation layers into multilayer networks
#              for chromosome 17 copy number variation analysis. Creates regulatory
#              connections between methylation and expression nodes to capture
#              epigenetic regulation of gene expression in the context of chromosomal
#              amplification (GAIN) vs disomic (DIS) conditions.
#
# Author: Yon Abimanyu
# Date: 2026-01-01
# Version: 1.0
#
# Biological Rationale:
#   Chromosome 17 amplification can affect both gene dosage (expression) and
#   local epigenetic patterns (methylation). By integrating these layers, we
#   can identify genes where copy number changes trigger coordinated epigenetic
#   remodeling and transcriptional responses. The stricter thresholds (2.3/2.4)
#   ensure we focus on the most robust regulatory relationships.
#
# Network Structure:
#   - Expression layer: gene-gene regulatory relationships (threshold 2.3)
#   - Methylation layer: CpG site co-methylation patterns (threshold 2.4)
#   - Interlayer edges: methylation → expression (regulatory direction)
#   - Only genes with status ≠ UNCHANGED in both layers are connected
#
# Input Requirements:
#   Prerequisites: Run GAIN/DIS network construction scripts first to generate:
#   - nodes_exp_gain_chr17_sthr2.3, nodes_exp_dis_chr17_sthr2.3
#   - nodes_met_gain_chr17_sthr2.4, nodes_met_dis_chr17_sthr2.4
#   - nodes_exp_gain_chr17_sthr2.3_GC, nodes_exp_dis_chr17_sthr2.3_DC
#   - nodes_met_gain_chr17_sthr2.4_GC, nodes_met_dis_chr17_sthr2.4_DC
#   - All corresponding edge tables (rels_*)
#
# Output Files:
#   GD Contrast (GAIN vs DIS):
#   - results/gain_dis/nodes_all_gain_GD_sthr2.3_2.4.csv
#   - results/gain_dis/nodes_all_dis_GD_sthr2.3_2.4.csv
#   - results/gain_dis/rels_all_gain_GD_sthr2.3_2.4.csv
#   - results/gain_dis/rels_all_dis_GD_sthr2.3_2.4.csv
#   
#   Control Contrasts (GC: GAIN vs Control, DC: DIS vs Control):
#   - results/gain_dis/nodes_all_gain_GC_sthr2.3_2.4.csv
#   - results/gain_dis/nodes_all_dis_DC_sthr2.3_2.4.csv
#   - results/gain_dis/rels_all_gain_GC_sthr2.3_2.4.csv
#   - results/gain_dis/rels_all_dis_DC_sthr2.3_2.4.csv
#
# Dependencies: dplyr, stringr, readr
# ==============================================================================

# === LOAD REQUIRED LIBRARIES ===
library(dplyr)
library(stringr)
library(readr)

# === CREATE OUTPUT DIRECTORY ===
dir.create("results/gain_dis", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# PART 1: INTERLAYER INTEGRATION - GAIN SAMPLES
# ==============================================================================
# Rationale: Integrate expression and methylation layers for GAIN samples
# to identify genes where copy number amplification triggers coordinated
# epigenetic and transcriptional changes.

# --- Extract Unique Genes from Each Layer ---
genes_expr_gain <- nodes_exp_gain_chr17_sthr2.3 %>%
  mutate(id = str_remove(id, "_expression")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
message("GAIN: Expression layer genes: ", nrow(genes_expr_gain))

genes_meth_gain <- nodes_met_gain_chr17_sthr2.4 %>%
  mutate(id = str_remove(id, "_methylation")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
message("GAIN: Methylation layer genes: ", nrow(genes_meth_gain))

# --- Identify Genes Present in Both Layers ---
expr_meth_gain <- intersect(genes_expr_gain$id, genes_meth_gain$id)
message("GAIN: Genes in both layers: ", length(expr_meth_gain))


# ==============================================================================
# SECTION 1.1: GD CONTRAST (GAIN vs DIS)
# ==============================================================================
# Rationale: This contrast identifies copy number-specific effects by
# comparing amplified vs disomic samples, isolating dosage-sensitive genes.

# --- Extract Status for Genes in Both Layers ---
expr_status_gain_GD <- nodes_exp_gain_chr17_sthr2.3 %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_gain) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_gain_GD <- nodes_met_gain_chr17_sthr2.4 %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_gain) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine Status Information ---
expr_meth_status_gain_GD <- expr_status_gain_GD %>%
  left_join(meth_status_gain_GD, by = "gene")

# --- Filter for Biologically Active Genes ---
# Critical: Only genes showing significant changes in BOTH layers
# Rationale: Coordinated methylation-expression changes in response to
# copy number gain indicate functional regulatory relationships
valid_genes_GD_gain <- expr_meth_status_gain_GD %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("GAIN GD: Genes with changes in both layers: ", length(valid_genes_GD_gain))
# Note: Very few genes due to strict threshold 2.3/2.4 - highly confident relationships

# --- Create Interlayer Edges ---
rels_meth_to_expr_gain_GD <- tibble(
  from = paste0(valid_genes_GD_gain, "_methylation"),
  to   = paste0(valid_genes_GD_gain, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (GD) ---
nodes_all_gain_GD_sthr2.3_2.4 <- bind_rows(
  nodes_exp_gain_chr17_sthr2.3 %>% mutate(layer = "expression"),
  nodes_met_gain_chr17_sthr2.4 %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%  # Focus on biologically active genes
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("GAIN GD: Total significant nodes: ", nrow(nodes_all_gain_GD_sthr2.3_2.4))

# Quality check
message("GAIN GD - Expression status distribution:")
print(table(nodes_exp_gain_chr17_sthr2.3$status))
message("GAIN GD - Methylation status distribution:")
print(table(nodes_met_gain_chr17_sthr2.4$status))
message("GAIN GD - Combined status distribution:")
print(table(nodes_all_gain_GD_sthr2.3_2.4$status))

# Export nodes
write_csv(
  nodes_all_gain_GD_sthr2.3_2.4,
  "results/gain_dis/nodes_all_gain_GD_sthr2.3_2.4.csv"
)

# --- Construct Complete Edge Table (GD) ---
rels_all_gain_GD_sthr2.3_2.4 <- bind_rows(
  rels_exp_gain_GD_sthr2.3,
  rels_met_gain_GD_sthr2.4,
  rels_meth_to_expr_gain_GD
) %>%
  distinct()

message("GAIN GD: Total edges: ", nrow(rels_all_gain_GD_sthr2.3_2.4))

write_csv(
  rels_all_gain_GD_sthr2.3_2.4,
  "results/gain_dis/rels_all_gain_GD_sthr2.3_2.4.csv"
)

# --- Network Integrity Validation (GD) ---
orphaned_nodes_gain_GD <- setdiff(
  nodes_all_gain_GD_sthr2.3_2.4$id,
  c(rels_all_gain_GD_sthr2.3_2.4$from, rels_all_gain_GD_sthr2.3_2.4$to)
)
message("GAIN GD - Orphaned nodes: ", length(orphaned_nodes_gain_GD))

missing_from_gain_GD <- setdiff(
  rels_all_gain_GD_sthr2.3_2.4$from,
  nodes_all_gain_GD_sthr2.3_2.4$id
)
missing_to_gain_GD <- setdiff(
  rels_all_gain_GD_sthr2.3_2.4$to,
  nodes_all_gain_GD_sthr2.3_2.4$id
)
message("GAIN GD - Missing nodes in edges: FROM=", length(missing_from_gain_GD),
        ", TO=", length(missing_to_gain_GD))


# ==============================================================================
# SECTION 1.2: GC CONTRAST (GAIN vs Control)
# ==============================================================================
# Rationale: Compare GAIN samples against normal control to identify
# tumor-specific changes associated with chr17 amplification.

# --- Extract Status for GC Contrast ---
expr_status_gain_GC <- nodes_exp_gain_chr17_sthr2.3_GC %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_gain) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_gain_GC <- nodes_met_gain_chr17_sthr2.4_GC %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_gain) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine and Filter ---
expr_meth_status_gain_GC <- expr_status_gain_GC %>%
  left_join(meth_status_gain_GC, by = "gene")

valid_genes_GC_gain <- expr_meth_status_gain_GC %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("GAIN GC: Genes with changes in both layers: ", length(valid_genes_GC_gain))

# --- Create Interlayer Edges (GC) ---
rels_meth_to_expr_gain_GC <- tibble(
  from = paste0(valid_genes_GC_gain, "_methylation"),
  to   = paste0(valid_genes_GC_gain, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (GC) ---
nodes_all_gain_GC_sthr2.3_2.4 <- bind_rows(
  nodes_exp_gain_chr17_sthr2.3_GC %>% mutate(layer = "expression"),
  nodes_met_gain_chr17_sthr2.4_GC %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("GAIN GC: Total significant nodes: ", nrow(nodes_all_gain_GC_sthr2.3_2.4))

message("GAIN GC - Status distribution:")
print(table(nodes_all_gain_GC_sthr2.3_2.4$status))

write_csv(
  nodes_all_gain_GC_sthr2.3_2.4,
  "results/gain_dis/nodes_all_gain_GC_sthr2.3_2.4.csv"
)

# --- Construct Complete Edge Table (GC) ---
rels_all_gain_GC_sthr2.3_2.4 <- bind_rows(
  rels_exp_gain_GC_sthr2.3,
  rels_met_gain_GC_sthr2.4,
  rels_meth_to_expr_gain_GC
) %>%
  distinct()

message("GAIN GC: Total edges: ", nrow(rels_all_gain_GC_sthr2.3_2.4))

write_csv(
  rels_all_gain_GC_sthr2.3_2.4,
  "results/gain_dis/rels_all_gain_GC_sthr2.3_2.4.csv"
)

# --- Network Integrity Validation (GC) ---
orphaned_nodes_gain_GC <- setdiff(
  nodes_all_gain_GC_sthr2.3_2.4$id,
  c(rels_all_gain_GC_sthr2.3_2.4$from, rels_all_gain_GC_sthr2.3_2.4$to)
)
message("GAIN GC - Orphaned nodes: ", length(orphaned_nodes_gain_GC))

missing_from_gain_GC <- setdiff(
  rels_all_gain_GC_sthr2.3_2.4$from,
  nodes_all_gain_GC_sthr2.3_2.4$id
)
missing_to_gain_GC <- setdiff(
  rels_all_gain_GC_sthr2.3_2.4$to,
  nodes_all_gain_GC_sthr2.3_2.4$id
)
message("GAIN GC - Missing nodes: FROM=", length(missing_from_gain_GC),
        ", TO=", length(missing_to_gain_GC))


# ==============================================================================
# PART 2: INTERLAYER INTEGRATION - DIS SAMPLES
# ==============================================================================

# --- Extract Unique Genes from Each Layer ---
genes_expr_dis <- nodes_exp_dis_chr17_sthr2.3 %>%
  mutate(id = str_remove(id, "_expression")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
message("DIS: Expression layer genes: ", nrow(genes_expr_dis))

genes_meth_dis <- nodes_met_dis_chr17_sthr2.4 %>%
  mutate(id = str_remove(id, "_methylation")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
message("DIS: Methylation layer genes: ", nrow(genes_meth_dis))

# --- Identify Genes Present in Both Layers ---
expr_meth_dis <- intersect(genes_expr_dis$id, genes_meth_dis$id)
message("DIS: Genes in both layers: ", length(expr_meth_dis))


# ==============================================================================
# SECTION 2.1: GD CONTRAST (DIS side of GAIN vs DIS)
# ==============================================================================

# --- Extract Status ---
expr_status_dis_GD <- nodes_exp_dis_chr17_sthr2.3 %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_dis) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_dis_GD <- nodes_met_dis_chr17_sthr2.4 %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_dis) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine and Filter ---
expr_meth_status_dis_GD <- expr_status_dis_GD %>%
  left_join(meth_status_dis_GD, by = "gene")

valid_genes_GD_dis <- expr_meth_status_dis_GD %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("DIS GD: Genes with changes in both layers: ", length(valid_genes_GD_dis))

# --- Create Interlayer Edges ---
rels_meth_to_expr_dis_GD <- tibble(
  from = paste0(valid_genes_GD_dis, "_methylation"),
  to   = paste0(valid_genes_GD_dis, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (GD) ---
nodes_all_dis_GD_sthr2.3_2.4 <- bind_rows(
  nodes_exp_dis_chr17_sthr2.3 %>% mutate(layer = "expression"),
  nodes_met_dis_chr17_sthr2.4 %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("DIS GD: Total significant nodes: ", nrow(nodes_all_dis_GD_sthr2.3_2.4))

message("DIS GD - Status distribution:")
print(table(nodes_all_dis_GD_sthr2.3_2.4$status))

write_csv(
  nodes_all_dis_GD_sthr2.3_2.4,
  "results/gain_dis/nodes_all_dis_GD_sthr2.3_2.4.csv"
)

# --- Construct Complete Edge Table (GD) ---
rels_all_dis_GD_sthr2.3_2.4 <- bind_rows(
  rels_exp_dis_GD_sthr2.3,
  rels_met_dis_GD_sthr2.4,
  rels_meth_to_expr_dis_GD
) %>%
  distinct()

message("DIS GD: Total edges: ", nrow(rels_all_dis_GD_sthr2.3_2.4))

write_csv(
  rels_all_dis_GD_sthr2.3_2.4,
  "results/gain_dis/rels_all_dis_GD_sthr2.3_2.4.csv"
)

# --- Validation (GD) ---
orphaned_nodes_dis_GD <- setdiff(
  nodes_all_dis_GD_sthr2.3_2.4$id,
  c(rels_all_dis_GD_sthr2.3_2.4$from, rels_all_dis_GD_sthr2.3_2.4$to)
)
message("DIS GD - Orphaned nodes: ", length(orphaned_nodes_dis_GD))

missing_from_dis_GD <- setdiff(
  rels_all_dis_GD_sthr2.3_2.4$from,
  nodes_all_dis_GD_sthr2.3_2.4$id
)
missing_to_dis_GD <- setdiff(
  rels_all_dis_GD_sthr2.3_2.4$to,
  nodes_all_dis_GD_sthr2.3_2.4$id
)
message("DIS GD - Missing nodes: FROM=", length(missing_from_dis_GD),
        ", TO=", length(missing_to_dis_GD))


# ==============================================================================
# SECTION 2.2: DC CONTRAST (DIS vs Control)
# ==============================================================================

# --- Extract Status ---
expr_status_dis_DC <- nodes_exp_dis_chr17_sthr2.3_DC %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_dis) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_dis_DC <- nodes_met_dis_chr17_sthr2.4_DC %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_dis) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine and Filter ---
expr_meth_status_dis_DC <- expr_status_dis_DC %>%
  left_join(meth_status_dis_DC, by = "gene")

valid_genes_DC_dis <- expr_meth_status_dis_DC %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("DIS DC: Genes with changes in both layers: ", length(valid_genes_DC_dis))

# --- Create Interlayer Edges ---
rels_meth_to_expr_dis_DC <- tibble(
  from = paste0(valid_genes_DC_dis, "_methylation"),
  to   = paste0(valid_genes_DC_dis, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (DC) ---
nodes_all_dis_DC_sthr2.3_2.4 <- bind_rows(
  nodes_exp_dis_chr17_sthr2.3_DC %>% mutate(layer = "expression"),
  nodes_met_dis_chr17_sthr2.4_DC %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("DIS DC: Total significant nodes: ", nrow(nodes_all_dis_DC_sthr2.3_2.4))

message("DIS DC - Status distribution:")
print(table(nodes_all_dis_DC_sthr2.3_2.4$status))

write_csv(
  nodes_all_dis_DC_sthr2.3_2.4,
  "results/gain_dis/nodes_all_dis_DC_sthr2.3_2.4.csv"
)

# --- Construct Complete Edge Table (DC) ---
rels_all_dis_DC_sthr2.3_2.4 <- bind_rows(
  rels_exp_dis_DC_sthr2.3,
  rels_met_dis_DC_sthr2.4,
  rels_meth_to_expr_dis_DC
) %>%
  distinct()

message("DIS DC: Total edges: ", nrow(rels_all_dis_DC_sthr2.3_2.4))

write_csv(
  rels_all_dis_DC_sthr2.3_2.4,
  "results/gain_dis/rels_all_dis_DC_sthr2.3_2.4.csv"
)

# --- Validation (DC) ---
orphaned_nodes_dis_DC <- setdiff(
  nodes_all_dis_DC_sthr2.3_2.4$id,
  c(rels_all_dis_DC_sthr2.3_2.4$from, rels_all_dis_DC_sthr2.3_2.4$to)
)
message("DIS DC - Orphaned nodes: ", length(orphaned_nodes_dis_DC))

missing_from_dis_DC <- setdiff(
  rels_all_dis_DC_sthr2.3_2.4$from,
  nodes_all_dis_DC_sthr2.3_2.4$id
)
missing_to_dis_DC <- setdiff(
  rels_all_dis_DC_sthr2.3_2.4$to,
  nodes_all_dis_DC_sthr2.3_2.4$id
)
message("DIS DC - Missing nodes: FROM=", length(missing_from_dis_DC),
        ", TO=", length(missing_to_dis_DC))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
# Output files generated in results/gain_dis/:
# 
# GAIN networks:
#   - nodes_all_gain_GD_sthr2.3_2.4.csv (632 nodes)
#   - rels_all_gain_GD_sthr2.3_2.4.csv (4,880 edges)
#   - nodes_all_gain_GC_sthr2.3_2.4.csv (1,441 nodes)
#   - rels_all_gain_GC_sthr2.3_2.4.csv (17,422 edges)
#
# DIS networks:
#   - nodes_all_dis_GD_sthr2.3_2.4.csv (675 nodes)
#   - rels_all_dis_GD_sthr2.3_2.4.csv (10,743 edges)
#   - nodes_all_dis_DC_sthr2.3_2.4.csv (1,474 nodes)
#   - rels_all_dis_DC_sthr2.3_2.4.csv (46,534 edges)
#
# Comparison with BRCA1 Networks:
#   - GAIN/DIS networks are generally sparser due to stricter thresholds
#     (2.3/2.4 vs 1.4) - higher confidence but fewer relationships
#   - Different biological focus: copy number dosage effects vs BRCA1
#     expression status
#   - Similar pattern: control contrasts have more edges than paired contrasts
# ==============================================================================
