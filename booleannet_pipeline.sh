#!/bin/bash

################################################################################
# BooleanNet Analysis Pipeline for Multi-Omics Cancer Research
################################################################################
#
# TITLE: Boolean Network Construction from Multi-Omics Data
#
# DESCRIPTION:
#   This pipeline constructs Boolean implication networks from discretized
#   gene expression (RNA-seq) and DNA methylation data to identify regulatory
#   relationships in breast cancer with chromosome 17q amplification.
#
#   The pipeline implements two sequential algorithms:
#   1. StepMiner: Discretizes continuous omics data into categorical states
#      (Low/Intermediate/High) based on optimal threshold identification
#   2. BooleanNet: Identifies significant Boolean implication relationships
#      between gene pairs using contingency table analysis
#
# AUTHOR: Yon Abimanyu
# DATE: 2026-01-01
# VERSION: 1.0
#
# INPUT FILES (required in data/ directory):
#   - *.pcl: Gene expression/methylation matrix (tab-delimited)
#   - *.bv: BitVector file from StepMiner discretization (3-level: 0,1,2)
#   - *.ph: This file is a placeholder; what matters is that the .ph command is provided as a template
#   - extract_met.pl: Modified extraction script with 'all' command
#
# OUTPUT FILES (generated in results/ directory):
#   - file.idx: Gene index mapping for data retrieval
#   - *.rl: Boolean network file containing implication relationships
#   - all_relations_*.txt: Complete list of Boolean implications across all gene pairs
#
# DEPENDENCIES:
#   - Java Runtime Environment (JRE) 8 or higher
#   - Perl 5.42.0
#   - stepminer-1.1.jar (BooleanNet implementation)
#   - extract_met.pl (custom modified script)
#   - Minimum 16GB RAM recommended for large datasets
#
# PARAMETERS:
#   - ERROR_RATE: 0.1 (10% - max proportion of samples in sparse quadrant)
#   - STAT_THRESHOLD: 2.4 (F-statistic threshold, optimized via permutation testing)
#   - SINGLE_THRESHOLD: 0.01 (1% - min proportion of high/low expression per gene)
#
# REFERENCE:
#   Sahoo et al. (2008). Boolean implication networks derived from large scale,
#   whole genome microarray datasets. Genome Biology, 9(10), R157.
#
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable

# ============================================================================
# CONFIGURATION
# ============================================================================

# Directory structure
DATA_DIR="data"
RESULTS_DIR="results"

# Input files (relative paths)
INPUT_PCL="met_dis_chr17.pcl"           # Gene expression/methylation matrix
INPUT_BV="met_dis_chr17.bv"             # BitVector from StepMiner
INPUT_PHENOTYPE="test_output.ph"         # A dummy file, which may also be regarded as a sample phenotype file

# Output files
OUTPUT_INDEX="file.idx"                  # Gene index file
OUTPUT_NETWORK="met_dis_chr17.rl"        # Boolean network file
OUTPUT_RELATIONS="all_relations_dismet_2.4.txt"  # All Boolean implications

# Dependencies
STEPMINER_JAR="stepminer-1.1.jar"
EXTRACT_SCRIPT="extract_met.pl"

# BooleanNet parameters (tuned via permutation testing - see Methods section 2.2)
ERROR_RATE="0.1"           # Max 10% samples in sparse quadrant
STAT_THRESHOLD="2.4"       # F-statistic threshold (optimized for FDR control)
SINGLE_THRESHOLD="0.01"    # Min 1% samples with high/low expression per gene

# Java memory allocation (adjust based on system capacity)
JAVA_MIN_HEAP="8192m"      # Initial heap size
JAVA_MAX_HEAP="16384m"     # Maximum heap size

# ============================================================================
# VALIDATION
# ============================================================================

echo "=== BooleanNet Analysis Pipeline ==="
echo "Validating environment and dependencies..."

# Check required files
if [ ! -f "$DATA_DIR/$INPUT_PCL" ]; then
    echo "ERROR: Input PCL file not found: $DATA_DIR/$INPUT_PCL"
    exit 1
fi

if [ ! -f "$DATA_DIR/$INPUT_BV" ]; then
    echo "ERROR: BitVector file not found: $DATA_DIR/$INPUT_BV"
    exit 1
fi

if [ ! -f "$STEPMINER_JAR" ]; then
    echo "ERROR: StepMiner JAR not found: $STEPMINER_JAR"
    exit 1
fi

if [ ! -f "$EXTRACT_SCRIPT" ]; then
    echo "ERROR: Extract script not found: $EXTRACT_SCRIPT"
    exit 1
fi

# Check Java
if ! command -v java &> /dev/null; then
    echo "ERROR: Java not found. Please install Java 8 or higher."
    exit 1
fi

# Check Perl
if ! command -v perl &> /dev/null; then
    echo "ERROR: Perl not found. Please install Perl 5.x."
    exit 1
fi

# Create output directory
mkdir -p "$RESULTS_DIR"

echo "Environment validation complete."
echo ""

# ============================================================================
# STEP 1: GENE INDEX GENERATION
# ============================================================================
# SCIENTIFIC RATIONALE:
# Gene indexing creates a mapping structure for efficient retrieval of gene
# information from the PCL matrix. This index enables rapid lookup of gene
# symbols, probe IDs, and their positions in the data matrix, which is critical
# for subsequent Boolean network construction where thousands of gene pair
# relationships must be evaluated.
# ============================================================================

echo "=== STEP 1: Generating Gene Index ==="
echo "Input: $DATA_DIR/$INPUT_PCL"
echo "Output: $RESULTS_DIR/$OUTPUT_INDEX"
echo ""

perl "$EXTRACT_SCRIPT" index "$DATA_DIR/$INPUT_PCL" > "$RESULTS_DIR/$OUTPUT_INDEX"

if [ $? -eq 0 ]; then
    echo "Gene index generated successfully."
    echo "Index entries: $(wc -l < "$RESULTS_DIR/$OUTPUT_INDEX")"
else
    echo "ERROR: Gene index generation failed."
    exit 1
fi
echo ""

# ============================================================================
# STEP 2: BOOLEAN NETWORK CONSTRUCTION - CONTINGENCY TABLE ANALYSIS
# ============================================================================
# SCIENTIFIC RATIONALE:
# The bitMatrix step constructs 2×2 contingency tables for all gene pairs
# using discretized Boolean states (LOW=0, HIGH=1) from StepMiner output.
# For each gene pair (A,B), the algorithm evaluates six types of Boolean
# implications based on sparse quadrant patterns:
#
# ASYMMETRIC RELATIONS (one sparse quadrant):
#   Type 1: A low  → B high  (n₀₀ sparse)
#   Type 2: A low  → B low   (n₀₁ sparse)
#   Type 3: A high → B high  (n₁₀ sparse)
#   Type 4: A high → B low   (n₁₁ sparse)
#
# SYMMETRIC RELATIONS (two sparse quadrants):
#   Type 5: Equivalent (n₀₁, n₁₀ sparse - secondary diagonal)
#   Type 6: Opposite   (n₀₀, n₁₁ sparse - main diagonal)
#
# Statistical significance is determined by:
# 1. Error rate: Max proportion of samples in sparse quadrant (default: 10%)
# 2. F-statistic: Between-group vs within-group variance ratio (threshold: 2.4)
# 3. Single threshold: Min proportion of samples with high/low expression (1%)
#
# PARAMETER JUSTIFICATION:
# - stat_thr=2.4: Optimized via permutation testing to control empirical FDR
#   while maintaining sufficient network connectivity for downstream analysis
# - error_rate=0.1: Standard tolerance allowing 10% violations in sparse quadrant
# - single_thr=0.01: Ensures adequate sample representation for reliable inference
# ============================================================================

echo "=== STEP 2: Constructing Boolean Network (Contingency Analysis) ==="
echo "Input BitVector: $DATA_DIR/$INPUT_BV"
echo "Output Network: $RESULTS_DIR/$OUTPUT_NETWORK"
echo ""
echo "Parameters:"
echo "  - Error Rate: $ERROR_RATE (max proportion in sparse quadrant)"
echo "  - Statistical Threshold: $STAT_THRESHOLD (F-statistic cutoff)"
echo "  - Single Threshold: $SINGLE_THRESHOLD (min expression proportion)"
echo ""

java -Xms"$JAVA_MIN_HEAP" -Xmx"$JAVA_MAX_HEAP" \
     -XX:+UseG1GC \
     -XX:MaxGCPauseMillis=200 \
     -cp "$STEPMINER_JAR" \
     tools.CustomAnalysis boolean bitMatrix \
     "$RESULTS_DIR/$OUTPUT_NETWORK" \
     "$DATA_DIR/$INPUT_BV" \
     "$DATA_DIR/$INPUT_PHENOTYPE" \
     All \
     "$ERROR_RATE" \
     "$STAT_THRESHOLD" \
     "$SINGLE_THRESHOLD"

if [ $? -eq 0 ]; then
    echo "Contingency table analysis complete."
else
    echo "ERROR: Boolean network construction failed."
    exit 1
fi
echo ""

# ============================================================================
# STEP 3: BOOLEAN NETWORK COMPLETION - MATRIX FILLING
# ============================================================================
# SCIENTIFIC RATIONALE:
# The bitMatrixFill step completes the Boolean network adjacency matrix by
# populating all identified relationships into a comprehensive network structure.
# This creates the directed adjacency matrix A where A_ij = 1 indicates a
# significant Boolean implication from gene i to gene j (including both
# asymmetric one-way and symmetric two-way relationships).
#
# This step is essential for:
# 1. Ensuring complete network topology representation
# 2. Enabling bidirectional relationship tracking for symmetric implications
# 3. Facilitating efficient graph database import and query operations
# ============================================================================

echo "=== STEP 3: Completing Boolean Network (Matrix Filling) ==="
echo "Processing: $RESULTS_DIR/$OUTPUT_NETWORK"
echo ""

java -Xms"$JAVA_MIN_HEAP" -Xmx"$JAVA_MAX_HEAP" \
     -cp "$STEPMINER_JAR" \
     tools.CustomAnalysis boolean bitMatrixFill \
     "$RESULTS_DIR/$OUTPUT_NETWORK"

if [ $? -eq 0 ]; then
    echo "Network matrix filled successfully."
else
    echo "ERROR: Matrix filling failed."
    exit 1
fi
echo ""

# ============================================================================
# STEP 4: BOOLEAN NETWORK STATISTICS CALCULATION
# ============================================================================
# SCIENTIFIC RATIONALE:
# The bitMatrixFillStats step computes comprehensive network statistics including:
# 1. Distribution of implication types (Types 1-6) across the network
# 2. Network density and connectivity metrics
# 3. Quality control metrics for validation
#
# These statistics are crucial for:
# - Validating network construction quality
# - Assessing predominant regulatory logic patterns
# - Comparing networks across experimental conditions
# - Informing downstream topological analysis
# ============================================================================

echo "=== STEP 4: Computing Network Statistics ==="
echo "Analyzing: $RESULTS_DIR/$OUTPUT_NETWORK"
echo ""

java -Xms"$JAVA_MIN_HEAP" -Xmx"$JAVA_MAX_HEAP" \
     -cp "$STEPMINER_JAR" \
     tools.CustomAnalysis boolean bitMatrixFillStats \
     "$RESULTS_DIR/$OUTPUT_NETWORK"

if [ $? -eq 0 ]; then
    echo "Network statistics computed successfully."
else
    echo "ERROR: Statistics calculation failed."
    exit 1
fi
echo ""

# ============================================================================
# STEP 5: EXTRACT ALL BOOLEAN RELATIONSHIPS
# ============================================================================
# SCIENTIFIC RATIONALE:
# The extract_met.pl script with 'all' command retrieves all significant Boolean
# implication relationships across the entire gene pair space in a single
# execution. This modified extraction approach (vs. the default per-gene query)
# is optimized for:
# 1. Comprehensive network export for multi-omics integration
# 2. Efficient bulk data retrieval without iterative queries
# 3. Direct compatibility with graph database import formats (Neo4j)
# 4. Facilitating reproducible analysis and data sharing
#
# Output format: Tab-delimited file containing:
# - Source gene, Target gene, Implication type (1-6), Statistical metrics
# ============================================================================

echo "=== STEP 5: Extracting All Boolean Relationships ==="
echo "Network: $RESULTS_DIR/$OUTPUT_NETWORK"
echo "Index: $RESULTS_DIR/$OUTPUT_INDEX"
echo "Output: $RESULTS_DIR/$OUTPUT_RELATIONS"
echo ""

perl "$EXTRACT_SCRIPT" all \
     "$RESULTS_DIR/$OUTPUT_NETWORK" \
     "$RESULTS_DIR/$OUTPUT_INDEX" \
     > "$RESULTS_DIR/$OUTPUT_RELATIONS"

if [ $? -eq 0 ]; then
    echo "Boolean relationships extracted successfully."
    echo "Total relationships: $(wc -l < "$RESULTS_DIR/$OUTPUT_RELATIONS")"
else
    echo "ERROR: Relationship extraction failed."
    exit 1
fi
echo ""

# ============================================================================
# PIPELINE COMPLETION
# ============================================================================

echo "=== Pipeline Complete ==="
echo ""
echo "Output files:"
echo "  - Gene index: $RESULTS_DIR/$OUTPUT_INDEX"
echo "  - Network file: $RESULTS_DIR/$OUTPUT_NETWORK"
echo "  - All relationships: $RESULTS_DIR/$OUTPUT_RELATIONS"
echo ""
echo "Next steps:"
echo "  1. Import $OUTPUT_RELATIONS to Neo4j graph database"
echo "  2. Construct multi-layer networks with methylation-expression edges"
echo "  3. Perform topological analysis in Cytoscape"
echo "  4. Extract ego-networks for target genes (e.g., BRCA1)"
echo ""
echo "For details on network integration and analysis, refer to Methods section 2.3"
echo ""

################################################################################
# NOTES FOR USAGE:
################################################################################
#
# 1. PARAMETER TUNING:
#    The stat_thr parameter (currently 2.4) should be optimized for your specific
#    dataset using permutation testing (50 iterations recommended) to control
#    empirical False Discovery Rate (FDR < 0.05). See Methods section 2.3.
#
# 2. MEMORY ALLOCATION:
#    Java heap sizes (-Xms/-Xmx) should be adjusted based on:
#    - Dataset size (number of genes × number of samples)
#    - System available RAM
#    - Recommended: allocate 50-75% of total RAM for large datasets
#
# 3. MULTIPLE CONDITION ANALYSIS:
#    To analyze different experimental conditions (GAIN-Chr17q, DIS-Chr17q,
#    upBRCA1, downBRCA1), run this pipeline separately for each condition's
#    BitVector file, adjusting INPUT_BV and output filenames accordingly.
#
# 4. FILE FORMAT REQUIREMENTS:
#    - PCL: Tab-delimited matrix with genes as rows, samples as columns
#    - BV: BitVector format from StepMiner (0=Low, 1=Intermediate, 2=High)
#    - PH: A dummy file, not required in the working directory, provided solely as a command template
#
# 5. TROUBLESHOOTING:
#    - OutOfMemoryError: Increase -Xmx parameter
#    - Slow performance: Enable G1GC garbage collector (already configured)
#    - Empty output: Check BitVector quality and parameter stringency
#
################################################################################