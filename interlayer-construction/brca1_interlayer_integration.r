# ==============================================================================
# MULTI-OMICS NETWORK INTEGRATION: BRCA1 ANALYSIS
# ==============================================================================
# Title: Interlayer Network Construction for BRCA1 Expression Groups
# Description: Integrates expression and methylation layers into multilayer networks
#              by creating regulatory connections between methylation and expression
#              nodes for the same genes. Focuses on genes showing significant changes
#              in BOTH omics layers to capture active regulatory relationships.
#              Generates complete network files (nodes + edges) for visualization
#              and analysis in Cytoscape, igraph, or other network tools.
#
# Author: [Your Name]
# Date: 2025-01-01
# Version: 3.0 (Updated 2025-12-24)
#
# Biological Rationale:
#   DNA methylation is a key epigenetic mechanism that regulates gene expression.
#   By linking methylation changes to expression changes for the same gene, we
#   can identify potential regulatory relationships where epigenetic modifications
#   drive transcriptional changes. Filtering for genes significant in BOTH layers
#   ensures we capture biologically active regulatory circuits rather than random
#   co-occurrences.
#
# Network Structure:
#   - Expression layer: gene-gene regulatory relationships
#   - Methylation layer: CpG site co-methylation patterns
#   - Interlayer edges: methylation → expression (regulatory direction)
#   - Only genes with status ≠ UNCHANGED in both layers are connected
#
# Input Requirements:
#   Prerequisites: Run BRCA1 network construction scripts first to generate:
#   - nodes_exp_upBRCA1_chr17_sthr1.40, nodes_exp_downBRCA1_chr17_sthr1.40
#   - nodes_met_upBRCA1_chr17_sthr1.40, nodes_met_downBRCA1_chr17_sthr1.40
#   - nodes_exp_UC_chr17_sthr1.40, nodes_exp_DC_chr17_sthr1.40
#   - nodes_met_UC_chr17_sthr1.40, nodes_met_DC_chr17_sthr1.40
#   - All corresponding edge tables (rels_*)
#
# Output Files:
#   UD Contrast (upBRCA1 vs downBRCA1):
#   - results/brca1/nodes_all_upBRCA1_UD_sthr1.40.csv
#   - results/brca1/nodes_all_downBRCA1_UD_sthr1.40.csv
#   - results/brca1/rels_all_upBRCA1_UD_sthr1.40.csv
#   - results/brca1/rels_all_downBRCA1_UD_sthr1.40.csv
#   
#   Control Contrasts (UC: upBRCA1 vs Control, DC: downBRCA1 vs Control):
#   - results/brca1/nodes_all_upBRCA1_UC_sthr1.40.csv
#   - results/brca1/nodes_all_downBRCA1_DC_sthr1.40.csv
#   - results/brca1/rels_all_upBRCA1_UC_sthr1.40.csv
#   - results/brca1/rels_all_downBRCA1_DC_sthr1.40.csv
#
# Dependencies: dplyr, stringr, readr
# ==============================================================================

# === LOAD REQUIRED LIBRARIES ===
library(dplyr)
library(stringr)
library(readr)

# === CREATE OUTPUT DIRECTORY ===
dir.create("results/brca1", recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# PART 1: INTERLAYER INTEGRATION - upBRCA1 SAMPLES
# ==============================================================================
# Rationale: Integrate expression and methylation layers for upBRCA1 samples
# to identify genes where epigenetic changes (methylation) may regulate
# transcriptional changes (expression).

# --- Extract Unique Genes from Each Layer ---
# Remove layer suffix to get base gene symbols for matching
genes_expr_upbrca1 <- nodes_exp_upBRCA1_chr17_sthr1.40 %>%
  mutate(id = str_remove(id, "_expression")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
# Expected: 702 genes (as of 2025-12-24)

genes_meth_upbrca1 <- nodes_met_upBRCA1_chr17_sthr1.40 %>%
  mutate(id = str_remove(id, "_methylation")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
# Expected: 1,594 genes (as of 2025-12-24)

# --- Identify Genes Present in Both Layers ---
# Only genes measured in both omics can have interlayer connections
expr_meth_upbrca1 <- intersect(genes_expr_upbrca1$id, genes_meth_upbrca1$id)
message("upBRCA1: Genes in both layers: ", length(expr_meth_upbrca1))
# Expected: 535 genes (as of 2025-12-24)


# ==============================================================================
# SECTION 1.1: UD CONTRAST (upBRCA1 vs downBRCA1)
# ==============================================================================
# Rationale: This contrast identifies BRCA1 expression-specific effects by
# comparing high vs low BRCA1 groups, controlling for tumor status.

# --- Extract Status for Genes in Both Layers ---
# Status indicates biological significance: UP_REG/DOWN_REG/UNCHANGED
expr_status_upbrca1_UD <- nodes_exp_upBRCA1_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_upbrca1) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_upbrca1_UD <- nodes_met_upBRCA1_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_upbrca1) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine Status Information ---
expr_meth_status_upbrca1_UD <- expr_status_upbrca1_UD %>%
  left_join(meth_status_upbrca1_UD, by = "gene")

# --- Filter for Biologically Active Genes ---
# Critical filtering: Only connect genes showing significant changes in BOTH layers
# Rationale: Genes with coordinated methylation-expression changes are more likely
# to represent functional regulatory relationships rather than random correlations
valid_genes_UD <- expr_meth_status_upbrca1_UD %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("upBRCA1 UD: Genes with changes in both layers: ", length(valid_genes_UD))
# Expected: 13 genes (as of 2025-12-24)

# --- Create Interlayer Edges ---
# Direction: methylation → expression (reflects regulatory causality)
# Methylation changes typically precede or drive expression changes
rels_meth_to_expr_upbrca1_UD <- tibble(
  from = paste0(valid_genes_UD, "_methylation"),
  to   = paste0(valid_genes_UD, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (UD) ---
# Combine expression and methylation nodes, keeping only significant changes
# Rationale: Focus network on genes showing biological response to BRCA1 status
nodes_all_upBRCA1_UD_sthr1.40 <- bind_rows(
  nodes_exp_upBRCA1_chr17_sthr1.40 %>% mutate(layer = "expression"),
  nodes_met_upBRCA1_chr17_sthr1.40 %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%  # Critical filter: biologically active only
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("upBRCA1 UD: Total significant nodes: ", nrow(nodes_all_upBRCA1_UD_sthr1.40))
# Expected: 388 nodes (as of 2025-12-24)

# Quality check: verify status distribution
message("upBRCA1 UD - Expression status distribution:")
print(table(nodes_exp_upBRCA1_chr17_sthr1.40$status))
message("upBRCA1 UD - Methylation status distribution:")
print(table(nodes_met_upBRCA1_chr17_sthr1.40$status))
message("upBRCA1 UD - Combined status distribution:")
print(table(nodes_all_upBRCA1_UD_sthr1.40$status))

# Export nodes
write_csv(
  nodes_all_upBRCA1_UD_sthr1.40,
  "results/brca1/nodes_all_upBRCA1_UD_sthr1.40.csv"
)

# --- Construct Complete Edge Table (UD) ---
# Combine all edge types: expression-expression, methylation-methylation, interlayer
rels_all_upBRCA1_UD_sthr1.40 <- bind_rows(
  rels_exp_upbrca1_UD_sthr1.40,
  rels_met_upbrca1_UD_sthr1.40,
  rels_meth_to_expr_upbrca1_UD
) %>%
  distinct()

message("upBRCA1 UD: Total edges: ", nrow(rels_all_upBRCA1_UD_sthr1.40))
# Expected: 13,609 edges (as of 2025-12-24)

# Export edges
write_csv(
  rels_all_upBRCA1_UD_sthr1.40,
  "results/brca1/rels_all_upBRCA1_UD_sthr1.40.csv"
)

# --- Network Integrity Validation (UD) ---
orphaned_nodes_UD <- setdiff(
  nodes_all_upBRCA1_UD_sthr1.40$id,
  c(rels_all_upBRCA1_UD_sthr1.40$from, rels_all_upBRCA1_UD_sthr1.40$to)
)
message("upBRCA1 UD - Orphaned nodes: ", length(orphaned_nodes_UD))

missing_from_UD <- setdiff(
  rels_all_upBRCA1_UD_sthr1.40$from,
  nodes_all_upBRCA1_UD_sthr1.40$id
)
missing_to_UD <- setdiff(
  rels_all_upBRCA1_UD_sthr1.40$to,
  nodes_all_upBRCA1_UD_sthr1.40$id
)
message("upBRCA1 UD - Missing nodes in edges: FROM=", length(missing_from_UD),
        ", TO=", length(missing_to_UD))


# ==============================================================================
# SECTION 1.2: UC CONTRAST (upBRCA1 vs Control)
# ==============================================================================
# Rationale: Compare upBRCA1 tumors against normal control to identify
# tumor-specific changes when BRCA1 is upregulated.

# --- Extract Status for UC Contrast ---
expr_status_upbrca1_UC <- nodes_exp_UC_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_upbrca1) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_upbrca1_UC <- nodes_met_UC_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_upbrca1) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine and Filter ---
expr_meth_status_upbrca1_UC <- expr_status_upbrca1_UC %>%
  left_join(meth_status_upbrca1_UC, by = "gene")

valid_genes_UC <- expr_meth_status_upbrca1_UC %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("upBRCA1 UC: Genes with changes in both layers: ", length(valid_genes_UC))
# Expected: 121 genes (as of 2025-12-24)

# --- Create Interlayer Edges (UC) ---
rels_meth_to_expr_upbrca1_UC <- tibble(
  from = paste0(valid_genes_UC, "_methylation"),
  to   = paste0(valid_genes_UC, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (UC) ---
nodes_all_upBRCA1_UC_sthr1.40 <- bind_rows(
  nodes_exp_UC_chr17_sthr1.40 %>% mutate(layer = "expression"),
  nodes_met_UC_chr17_sthr1.40 %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("upBRCA1 UC: Total significant nodes: ", nrow(nodes_all_upBRCA1_UC_sthr1.40))
# Expected: 1,281 nodes (as of 2025-12-24)

message("upBRCA1 UC - Expression status distribution:")
print(table(nodes_exp_UC_chr17_sthr1.40$status))
message("upBRCA1 UC - Methylation status distribution:")
print(table(nodes_met_UC_chr17_sthr1.40$status))
message("upBRCA1 UC - Combined status distribution:")
print(table(nodes_all_upBRCA1_UC_sthr1.40$status))

write_csv(
  nodes_all_upBRCA1_UC_sthr1.40,
  "results/brca1/nodes_all_upBRCA1_UC_sthr1.40.csv"
)

# --- Construct Complete Edge Table (UC) ---
rels_all_upBRCA1_UC_sthr1.40 <- bind_rows(
  rels_exp_upbrca1_UC_sthr1.40,
  rels_met_upbrca1_UC_sthr1.40,
  rels_meth_to_expr_upbrca1_UC
) %>%
  distinct()

message("upBRCA1 UC: Total edges: ", nrow(rels_all_upBRCA1_UC_sthr1.40))
# Expected: 153,118 edges (as of 2025-12-24)

write_csv(
  rels_all_upBRCA1_UC_sthr1.40,
  "results/brca1/rels_all_upBRCA1_UC_sthr1.40.csv"
)

# --- Network Integrity Validation (UC) ---
orphaned_nodes_UC <- setdiff(
  nodes_all_upBRCA1_UC_sthr1.40$id,
  c(rels_all_upBRCA1_UC_sthr1.40$from, rels_all_upBRCA1_UC_sthr1.40$to)
)
message("upBRCA1 UC - Orphaned nodes: ", length(orphaned_nodes_UC))

missing_from_UC <- setdiff(
  rels_all_upBRCA1_UC_sthr1.40$from,
  nodes_all_upBRCA1_UC_sthr1.40$id
)
missing_to_UC <- setdiff(
  rels_all_upBRCA1_UC_sthr1.40$to,
  nodes_all_upBRCA1_UC_sthr1.40$id
)
message("upBRCA1 UC - Missing nodes in edges: FROM=", length(missing_from_UC),
        ", TO=", length(missing_to_UC))


# ==============================================================================
# PART 2: INTERLAYER INTEGRATION - downBRCA1 SAMPLES
# ==============================================================================

# --- Extract Unique Genes from Each Layer ---
genes_expr_downbrca1 <- nodes_exp_downBRCA1_chr17_sthr1.40 %>%
  mutate(id = str_remove(id, "_expression")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
# Expected: 937 genes

genes_meth_downbrca1 <- nodes_met_downBRCA1_chr17_sthr1.40 %>%
  mutate(id = str_remove(id, "_methylation")) %>%
  select(id) %>%
  filter(!is.na(id)) %>%
  distinct() %>%
  arrange(id)
# Expected: 1,471 genes

# --- Identify Genes Present in Both Layers ---
expr_meth_downbrca1 <- intersect(genes_expr_downbrca1$id, genes_meth_downbrca1$id)
message("downBRCA1: Genes in both layers: ", length(expr_meth_downbrca1))
# Expected: 664 genes


# ==============================================================================
# SECTION 2.1: UD CONTRAST (downBRCA1 side of upBRCA1 vs downBRCA1)
# ==============================================================================

# --- Extract Status ---
expr_status_downbrca1_UD <- nodes_exp_downBRCA1_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_downbrca1) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_downbrca1_UD <- nodes_met_downBRCA1_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_downbrca1) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine and Filter ---
expr_meth_status_downbrca1_UD <- expr_status_downbrca1_UD %>%
  left_join(meth_status_downbrca1_UD, by = "gene")

valid_genes_UD_down <- expr_meth_status_downbrca1_UD %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("downBRCA1 UD: Genes with changes in both layers: ", length(valid_genes_UD_down))
# Expected: 13 genes

# --- Create Interlayer Edges ---
rels_meth_to_expr_downbrca1_UD <- tibble(
  from = paste0(valid_genes_UD_down, "_methylation"),
  to   = paste0(valid_genes_UD_down, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (UD) ---
nodes_all_downBRCA1_UD_sthr1.40 <- bind_rows(
  nodes_exp_downBRCA1_chr17_sthr1.40 %>% mutate(layer = "expression"),
  nodes_met_downBRCA1_chr17_sthr1.40 %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("downBRCA1 UD: Total significant nodes: ", nrow(nodes_all_downBRCA1_UD_sthr1.40))
# Expected: 377 nodes

message("downBRCA1 UD - Status distribution:")
print(table(nodes_all_downBRCA1_UD_sthr1.40$status))

write_csv(
  nodes_all_downBRCA1_UD_sthr1.40,
  "results/brca1/nodes_all_downBRCA1_UD_sthr1.40.csv"
)

# --- Construct Complete Edge Table (UD) ---
rels_all_downBRCA1_UD_sthr1.40 <- bind_rows(
  rels_exp_downBRCA1_UD_sthr1.40,
  rels_met_downbrca1_UD_sthr1.40,
  rels_meth_to_expr_downbrca1_UD
) %>%
  distinct()

message("downBRCA1 UD: Total edges: ", nrow(rels_all_downBRCA1_UD_sthr1.40))
# Expected: 15,321 edges

write_csv(
  rels_all_downBRCA1_UD_sthr1.40,
  "results/brca1/rels_all_downBRCA1_UD_sthr1.40.csv"
)

# --- Validation (UD) ---
orphaned_nodes_down_UD <- setdiff(
  nodes_all_downBRCA1_UD_sthr1.40$id,
  c(rels_all_downBRCA1_UD_sthr1.40$from, rels_all_downBRCA1_UD_sthr1.40$to)
)
message("downBRCA1 UD - Orphaned nodes: ", length(orphaned_nodes_down_UD))

missing_from_down_UD <- setdiff(
  rels_all_downBRCA1_UD_sthr1.40$from,
  nodes_all_downBRCA1_UD_sthr1.40$id
)
missing_to_down_UD <- setdiff(
  rels_all_downBRCA1_UD_sthr1.40$to,
  nodes_all_downBRCA1_UD_sthr1.40$id
)
message("downBRCA1 UD - Missing nodes: FROM=", length(missing_from_down_UD),
        ", TO=", length(missing_to_down_UD))


# ==============================================================================
# SECTION 2.2: DC CONTRAST (downBRCA1 vs Control)
# ==============================================================================

# --- Extract Status ---
expr_status_downbrca1_DC <- nodes_exp_DC_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_expression")) %>%
  filter(gene %in% expr_meth_downbrca1) %>%
  transmute(
    gene,
    status_expr = status
  ) %>%
  distinct() %>%
  arrange(gene)

meth_status_downbrca1_DC <- nodes_met_DC_chr17_sthr1.40 %>%
  mutate(gene = str_remove(id, "_methylation")) %>%
  filter(gene %in% expr_meth_downbrca1) %>%
  transmute(
    gene,
    status_meth = status
  ) %>%
  distinct() %>%
  arrange(gene)

# --- Combine and Filter ---
expr_meth_status_downbrca1_DC <- expr_status_downbrca1_DC %>%
  left_join(meth_status_downbrca1_DC, by = "gene")

valid_genes_DC <- expr_meth_status_downbrca1_DC %>%
  filter(
    status_expr != "UNCHANGED",
    status_meth != "UNCHANGED"
  ) %>%
  pull(gene)

message("downBRCA1 DC: Genes with changes in both layers: ", length(valid_genes_DC))
# Expected: 106 genes

# --- Create Interlayer Edges ---
rels_meth_to_expr_downbrca1_DC <- tibble(
  from = paste0(valid_genes_DC, "_methylation"),
  to   = paste0(valid_genes_DC, "_expression"),
  relation_type = "interlayer",
  layer = "expr_meth"
)

# --- Construct Complete Node Table (DC) ---
nodes_all_downBRCA1_DC_sthr1.40 <- bind_rows(
  nodes_exp_DC_chr17_sthr1.40 %>% mutate(layer = "expression"),
  nodes_met_DC_chr17_sthr1.40 %>% mutate(layer = "methylation")
) %>%
  filter(status != "UNCHANGED") %>%
  select(any_of(c(
    "id", "symbol", "entrez_id", "ensembl_id", "probe_id",
    "type", "layer", "chromosome", "cpg_location",
    "status", "source"
  ))) %>%
  distinct(id, .keep_all = TRUE)

message("downBRCA1 DC: Total significant nodes: ", nrow(nodes_all_downBRCA1_DC_sthr1.40))
# Expected: 1,225 nodes

message("downBRCA1 DC - Status distribution:")
print(table(nodes_all_downBRCA1_DC_sthr1.40$status))

write_csv(
  nodes_all_downBRCA1_DC_sthr1.40,
  "results/brca1/nodes_all_downBRCA1_DC_sthr1.40.csv"
)

# --- Construct Complete Edge Table (DC) ---
rels_all_downBRCA1_DC_sthr1.40 <- bind_rows(
  rels_exp_downbrca1_DC_sthr1.40,
  rels_met_downbrca1_DC_sthr1.40,
  rels_meth_to_expr_downbrca1_DC
) %>%
  distinct()

message("downBRCA1 DC: Total edges: ", nrow(rels_all_downBRCA1_DC_sthr1.40))
# Expected: 154,384 edges

write_csv(
  rels_all_downBRCA1_DC_sthr1.40,
  "results/brca1/rels_all_downBRCA1_DC_sthr1.40.csv"
)

# --- Validation (DC) ---
orphaned_nodes_down_DC <- setdiff(
  nodes_all_downBRCA1_DC_sthr1.40$id,
  c(rels_all_downBRCA1_DC_sthr1.40$from, rels_all_downBRCA1_DC_sthr1.40$to)
)
message("downBRCA1 DC - Orphaned nodes: ", length(orphaned_nodes_down_DC))

missing_from_down_DC <- setdiff(
  rels_all_downBRCA1_DC_sthr1.40$from,
  nodes_all_downBRCA1_DC_sthr1.40$id
)
missing_to_down_DC <- setdiff(
  rels_all_downBRCA1_DC_sthr1.40$to,
  nodes_all_downBRCA1_DC_sthr1.40$id
)
message("downBRCA1 DC - Missing nodes: FROM=", length(missing_from_down_DC),
        ", TO=", length(missing_to_down_DC))

# ==============================================================================
# END OF SCRIPT
# ==============================================================================
# Output files generated in results/brca1/:
# 
# upBRCA1 networks:
#   - nodes_all_upBRCA1_UD_sthr1.40.csv (388 nodes)
#   - rels_all_upBRCA1_UD_sthr1.40.csv (13,609 edges)
#   - nodes_all_upBRCA1_UC_sthr1.40.csv (1,281 nodes)
#   - rels_all_upBRCA1_UC_sthr1.40.csv (153,118 edges)
#
# downBRCA1 networks:
#   - nodes_all_downBRCA1_UD_sthr1.40.csv (377 nodes)
#   - rels_all_downBRCA1_UD_sthr1.40.csv (15,321 edges)
#   - nodes_all_downBRCA1_DC_sthr1.40.csv (1,225 nodes)
#   - rels_all_downBRCA1_DC_sthr1.40.csv (154,384 edges)
#
# Network Characteristics:
#   - UD contrasts are sparser (fewer nodes/edges) due to stricter filtering
#     between similar tumor groups
#   - Control contrasts have more edges reflecting broader dysregulation
#     vs normal tissue
#   - Interlayer edges range from 13 (UD) to 121 (UC) showing varying
#     degrees of coordinated methylation-expression changes
#
# Next Steps:
#   1. Import networks into Cytoscape for visualization
#   2. Perform network topology analysis (centrality, clustering)
#   3. Identify key regulatory hubs with high interlayer connectivity
#   4. Compare UD vs UC/DC to distinguish BRCA1-specific vs general tumor effects
# ==============================================================================