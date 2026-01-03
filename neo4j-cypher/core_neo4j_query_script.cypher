// ===== FLEXIBLE STATUS BRCA1 NETWORK EXPANSION - STRICT 2-LEVEL + SINGLE SOURCE EXPRESSION =====
// Strategy: Strict 2-level methylation depth + each expression node has only 1 source
// MODIFICATION: Flexible status filters (HYPERMETH/HYPOMETH for methylation, UP_REG/DOWN_REG for expression)

// Step 1: Start from BRCA1_methylation (NO STATUS RESTRICTION)
MATCH (start:BioNode {id: 'BRCA1_methylation', chromosome: '17', layer: 'methylation'})

// Step 2: Find direct methylation neighbors (level 1) - FLEXIBLE STATUS
OPTIONAL MATCH (start)-[meth_rel1:METHYLATION_REL]-(meth_level1:BioNode)
WHERE meth_level1.layer = 'methylation' 
  AND meth_level1.chromosome = '17'
  AND meth_level1.status IN ['HYPERMETH', 'HYPOMETH']  // ✓ FLEXIBLE: Both statuses allowed
  AND meth_level1 <> start  // No self-loops

// Step 3: Find methylation neighbors of level 1 (level 2) - FLEXIBLE STATUS
OPTIONAL MATCH (meth_level1)-[meth_rel2:METHYLATION_REL]-(meth_level2:BioNode)
WHERE meth_level2.layer = 'methylation' 
  AND meth_level2.chromosome = '17'
  AND meth_level2.status IN ['HYPERMETH', 'HYPOMETH']  // ✓ FLEXIBLE: Both statuses allowed
  AND meth_level2 <> start
  AND meth_level2 <> meth_level1  // No self-loops or duplicates

// Step 4: Collect methylation nodes with strict depth control
WITH start, 
     // Level 0: Start node only
     [start] AS level0_nodes,
     // Level 1: Direct neighbors of start
     COALESCE(collect(DISTINCT meth_level1), []) AS level1_nodes,
     // Level 2: Neighbors of level 1 (but not start or level 1)
     COALESCE(collect(DISTINCT meth_level2), []) AS level2_nodes,
     // Relationships: Level 0->1 and Level 1->2 only
     COALESCE(collect(DISTINCT meth_rel1), []) AS level01_rels,
     COALESCE(collect(DISTINCT meth_rel2), []) AS level12_rels

// Step 5: Create strict depth-controlled methylation network
WITH level0_nodes + level1_nodes + level2_nodes AS all_meth_nodes,
     level01_rels + level12_rels AS valid_meth_rels,
     level0_nodes, level1_nodes, level2_nodes

// Step 6: Find expression targets with depth priority (shortest path first) - FLEXIBLE STATUS
UNWIND all_meth_nodes AS meth_node
OPTIONAL MATCH (meth_node)-[direct_rel:INTERLAYER_EXPR_METH]->(expr_node:BioNode)
WHERE expr_node.layer = 'expression' 
  AND expr_node.chromosome = '17'
  AND expr_node.status IN ['UP_REG', 'DOWN_REG']  // ✓ FLEXIBLE: Both regulation types allowed

// Step 6b: Calculate exact depth for each methylation->expression connection
WITH all_meth_nodes, valid_meth_rels, level0_nodes, level1_nodes, level2_nodes,
     collect({
       meth_node: meth_node,
       expr_node: expr_node,
       direct_rel: direct_rel,
       depth: CASE 
         WHEN meth_node IN level0_nodes THEN 0  // Direct from BRCA1_methylation (depth 0)
         WHEN meth_node IN level1_nodes THEN 1  // Via 1 methylation hop (depth 1) 
         WHEN meth_node IN level2_nodes THEN 2  // Via 2 methylation hops (depth 2)
         ELSE 999  // Should not happen, but safety check
       END
     }) AS path_candidates

// Step 6c: Filter out null connections first
WITH all_meth_nodes, valid_meth_rels, level0_nodes, level1_nodes, level2_nodes,
     [item IN path_candidates WHERE item.expr_node IS NOT NULL] AS valid_connections

// Step 6d: For each expression node, find minimum depth using simpler approach
UNWIND valid_connections AS conn
WITH all_meth_nodes, valid_meth_rels, 
     conn.expr_node AS expr_node,
     collect({meth_node: conn.meth_node, direct_rel: conn.direct_rel, depth: conn.depth}) AS all_sources_for_expr

// Step 6e: Keep only the minimum depth sources for each expression
WITH all_meth_nodes, valid_meth_rels,
     expr_node,
     [item IN all_sources_for_expr WHERE item.depth = [src IN all_sources_for_expr | src.depth][0..1][0]] AS min_depth_sources

// Step 6f: From minimum depth sources, keep only the FIRST one (single source rule)
WITH all_meth_nodes, valid_meth_rels,
     collect({
       expr_node: expr_node,
       meth_node: min_depth_sources[0].meth_node,
       direct_rel: min_depth_sources[0].direct_rel
     }) AS final_single_sources

// Step 7: Extract the final valid connections (guaranteed single source per expression)
WITH all_meth_nodes, valid_meth_rels,
     [item IN final_single_sources | item.meth_node] AS valid_meth_for_expr,
     [item IN final_single_sources | item.expr_node] AS final_expr_nodes,
     [item IN final_single_sources | item.direct_rel] AS final_expr_rels

// Step 8: Filter methylation nodes - only those in paths to expression or start
WITH [node IN all_meth_nodes WHERE 
       node.id = 'BRCA1_methylation' OR  // Always keep start
       node IN valid_meth_for_expr       // Keep only if leads to expression
     ] AS filtered_meth_nodes,
     final_expr_nodes,
     valid_meth_rels,
     final_expr_rels

// Step 9: Filter methylation relationships - only between valid nodes
WITH filtered_meth_nodes + final_expr_nodes AS all_final_nodes,
     [rel IN valid_meth_rels WHERE 
       startNode(rel) IN filtered_meth_nodes 
       AND endNode(rel) IN filtered_meth_nodes
     ] AS filtered_meth_rels,
     final_expr_rels

// Step 10: Combine all relationships with enhanced deduplication
WITH all_final_nodes, 
     filtered_meth_rels + final_expr_rels AS all_final_rels

// Step 10b: Triple-layer deduplication using multiple keys - FLEXIBLE STATUS VALIDATION
UNWIND all_final_rels AS rel
WITH all_final_nodes, rel, startNode(rel) AS source, endNode(rel) AS target
WHERE source IN all_final_nodes 
  AND target IN all_final_nodes
  AND source <> target  // No self-loops
  AND source.chromosome = '17' 
  AND target.chromosome = '17'
  AND (source.layer <> 'expression')  // Expression nodes cannot be sources
  AND (source.layer <> 'methylation' OR source.status IN ['HYPERMETH', 'HYPOMETH'])  // ✓ FLEXIBLE
  AND (target.layer <> 'methylation' OR target.status IN ['HYPERMETH', 'HYPOMETH'])  // ✓ FLEXIBLE

// Step 10c: Create multiple deduplication keys for extra safety
WITH rel, source, target,
     // Primary key: source_id + relationship_type + target_id
     source.id + "##" + type(rel) + "##" + target.id AS primary_key,
     // Secondary key: using elementId for safety
     elementId(source) + "##" + type(rel) + "##" + elementId(target) AS secondary_key,
     // Tertiary key: using both id and elementId combined
     source.id + "||" + elementId(source) + "##" + type(rel) + "##" + target.id + "||" + elementId(target) AS tertiary_key

// Step 10d: Collect with all deduplication keys
WITH collect({
  rel: rel,
  source: source,
  target: target,
  primary_key: primary_key,
  secondary_key: secondary_key,
  tertiary_key: tertiary_key
}) AS all_rel_data

// Step 10e: Multi-level deduplication - primary key first, then secondary, then tertiary
WITH [item IN all_rel_data WHERE 
  item.primary_key = [other IN all_rel_data WHERE other.primary_key = item.primary_key][0].primary_key
  AND item.secondary_key = [other IN all_rel_data WHERE other.primary_key = item.primary_key AND other.secondary_key = item.secondary_key][0].secondary_key
  AND item.tertiary_key = [other IN all_rel_data WHERE other.primary_key = item.primary_key AND other.secondary_key = item.secondary_key AND other.tertiary_key = item.tertiary_key][0].tertiary_key
] AS unique_rel_data

// Step 11: Extract final unique relationships
UNWIND unique_rel_data AS unique_item
WITH unique_item.rel AS rel,
     unique_item.source AS source,
     unique_item.target AS target

// Step 12: Pattern classification (simplified without metapathway)
WITH rel, source, target,
     CASE 
       WHEN source.layer = target.layer AND source.layer = 'methylation' THEN 'METHYLATION_INTRA'
       WHEN source.layer = 'methylation' AND target.layer = 'expression' THEN 'METH_TO_EXPR_DIRECT'
       ELSE 'OTHER'
     END AS bridging_pattern,
     
     CASE 
       WHEN split(source.id, '_')[0] = split(target.id, '_')[0] THEN 'SAME_GENE'
       ELSE 'DIFFERENT_GENE'
     END AS gene_relationship

// Step 13: Advanced pattern analysis (simplified)
WITH rel, source, target, bridging_pattern, gene_relationship,
     CASE 
       WHEN bridging_pattern = 'METH_TO_EXPR_DIRECT' AND gene_relationship = 'SAME_GENE' THEN 'SAME_GENE_DIRECT'
       WHEN bridging_pattern = 'METH_TO_EXPR_DIRECT' AND gene_relationship = 'DIFFERENT_GENE' THEN 'CROSS_GENE_DIRECT'
       WHEN bridging_pattern = 'METHYLATION_INTRA' THEN 'METHYLATION_NETWORK'
       ELSE 'OTHER_PATTERN'
     END AS detailed_pattern

// Step 14: Final output with strict validation
RETURN 
  // Node identifiers for Cytoscape
  source.id AS source_node_id,
  target.id AS target_node_id,
  
  // Edge information
  type(rel) AS edge_type,
  rel.type AS edge_subtype,
  bridging_pattern AS edge_category,
  detailed_pattern AS pattern_analysis,
  gene_relationship AS gene_relation,
  
  // Source node attributes
  split(source.id, '_')[0] AS source_gene,
  source.layer AS source_layer,
  source.status AS source_status,
  source.chromosome AS source_chromosome,
  
  // Target node attributes
  split(target.id, '_')[0] AS target_gene,
  target.layer AS target_layer,
  target.status AS target_status,
  target.chromosome AS target_chromosome,
  
  // Connection type analysis
  CASE 
    WHEN source.layer = 'methylation' AND target.layer = 'expression' THEN 'METHYLATION_TO_EXPRESSION'
    WHEN source.layer = 'methylation' AND target.layer = 'methylation' THEN 'METHYLATION_NETWORK'
    ELSE 'OTHER_CONNECTION'
  END AS connection_type,
  
  // Debug: Show the depth level for validation
  CASE 
    WHEN source.id = 'BRCA1_methylation' THEN 'LEVEL_0_START'
    WHEN EXISTS((:BioNode {id: 'BRCA1_methylation'})-[:METHYLATION_REL]-(source)) THEN 'LEVEL_1_DIRECT'
    ELSE 'LEVEL_2_INDIRECT'
  END AS methylation_depth_level

ORDER BY edge_category, pattern_analysis, source_gene, target_gene;