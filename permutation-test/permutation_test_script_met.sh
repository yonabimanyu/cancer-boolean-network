#!/bin/bash

# ENHANCED BooleanNet Permutation Test Framework - Probe ID Version
# Version 8.0 - Full processing without timeouts or sampling
# Purpose: Find optimal statistical threshold through permutation testing
# Modified for DNA methylation probe IDs instead of Ensembl gene IDs

set -euo pipefail

# ============================================================================
# CONFIGURATION & PARAMETERS
# ============================================================================

# Input files (must exist in directory)
INPUT_BV="met_upbrca1.bv"
INPUT_IDX="file.idx"
EXTRACT_SCRIPT="extract_permutation.pl"
STEPMINER_JAR="stepminer-1.1.jar"

# Fixed parameters (as per guide)
FIXED_ERROR_RATE=0.1
FIXED_SINGLE_THR=0.01
OUTPUT_PH="test_output.ph"  # Template parameter for bitMatrix command

# Permutation test parameters
STAT_THRESHOLDS=(1.5 1.6 1.7 1.9 2.2)
N_PERMUTATIONS=50  # Number of permutations for FDR calculation

# Output directories and files
RESULTS_DIR="permutation_results"
ORIGINAL_DIR="${RESULTS_DIR}/original"
PERMUTED_DIR="${RESULTS_DIR}/permuted"
LOGS_DIR="${RESULTS_DIR}/logs"
FINAL_REPORT="${RESULTS_DIR}/fdr_analysis_report.tsv"
LOCK_DIR="${RESULTS_DIR}/.locks"

# Performance settings
TOTAL_CORES=$(nproc)
CORES=$(( (TOTAL_CORES * 8 + 5) / 10 ))
if [ $CORES -lt 1 ]; then CORES=1; fi
if [ $CORES -gt $TOTAL_CORES ]; then CORES=$TOTAL_CORES; fi

# Java memory settings
JAVA_OPTS="-Xms8192m -Xmx16384m -XX:+UseG1GC -XX:MaxGCPauseMillis=200"

# Colors and formatting
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
BOLD='\033[1m'
NC='\033[0m'

# ============================================================================
# DIRECTORY SETUP
# ============================================================================

# Create directories immediately
mkdir -p "$RESULTS_DIR" "$ORIGINAL_DIR" "$PERMUTED_DIR" "$LOGS_DIR" "$LOCK_DIR"

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1" | tee -a "$LOGS_DIR/main.log"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1" | tee -a "$LOGS_DIR/main.log"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1" | tee -a "$LOGS_DIR/main.log"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOGS_DIR/main.log"
}

# Simple progress tracking
show_progress() {
    local current=$1
    local total=$2
    local operation=$3
    local start_time=$4
    
    local current_time=$(date +%s)
    local elapsed=$((current_time - start_time))
    local progress_percent=$(( (current * 100) / total ))
    
    echo -e "${CYAN}📊 $operation: $current/$total (${progress_percent}%) | Elapsed: ${elapsed}s${NC}"
}

# File locking functions to prevent race conditions
acquire_lock() {
    local lockfile="$LOCK_DIR/$1.lock"
    
    while ! mkdir "$lockfile" 2>/dev/null; do
        sleep 1
    done
    
    return 0
}

release_lock() {
    local lockfile="$LOCK_DIR/$1.lock"
    rmdir "$lockfile" 2>/dev/null || true
}

# Safe file append with locking
safe_append() {
    local file=$1
    local content=$2
    local lockname=$(basename "$file")
    
    acquire_lock "$lockname"
    echo -e "$content" >> "$file"
    release_lock "$lockname"
    return 0
}

# ============================================================================
# CHECKPOINT AND RESUME FUNCTIONS
# ============================================================================

save_checkpoint() {
    local checkpoint_file="$RESULTS_DIR/checkpoint.txt"
    local current_stage=$1
    local current_item=$2
    
    echo "STAGE=$current_stage" > "$checkpoint_file"
    echo "ITEM=$current_item" >> "$checkpoint_file"
    echo "TIMESTAMP=$(date)" >> "$checkpoint_file"
    
    log_info "💾 Checkpoint saved: $current_stage - $current_item"
}

check_existing_results() {
    local result_type=$1  # "original" or "permuted"
    local threshold=$2
    local perm_id=${3:-""}  # Optional for permuted
    
    if [ "$result_type" = "original" ]; then
        local rl_file="$ORIGINAL_DIR/thresh_$threshold/original_thresh_$threshold.rl"
        if [ -f "$rl_file" ]; then
            log_info "⏭️  Skipping original threshold $threshold - already completed"
            return 0  # Already exists
        fi
    else
        local rl_file="$PERMUTED_DIR/thresh_$threshold/perm_${perm_id}_thresh_$threshold.rl"
        if [ -f "$rl_file" ]; then
            log_info "⏭️  Skipping permutation $perm_id threshold $threshold - already completed"
            return 0  # Already exists
        fi
    fi
    
    return 1  # Needs to be processed
}

# ============================================================================
# VALIDATION FUNCTIONS
# ============================================================================

validate_environment() {
    log_info "🔍 Validating environment and input files..."
    
    # Check required files
    local required_files=("$INPUT_BV" "$INPUT_IDX" "$EXTRACT_SCRIPT" "$STEPMINER_JAR")
    for file in "${required_files[@]}"; do
        if [ ! -f "$file" ]; then
            log_error "Required file not found: $file"
            exit 1
        fi
    done
    
    # Check Java availability
    if ! command -v java >/dev/null 2>&1; then
        log_error "Java not found in PATH"
        exit 1
    fi
    
    # Check Perl availability
    if ! command -v perl >/dev/null 2>&1; then
        log_error "Perl not found in PATH"
        exit 1
    fi
    
    # Check Python3 availability
    if ! command -v python3 >/dev/null 2>&1; then
        log_error "Python3 not found in PATH"
        exit 1
    fi
    
    # Validate BV file format - Updated for probe_id
    local header_line=$(head -n 1 "$INPUT_BV")
    if [[ ! "$header_line" =~ probe_id.*gene_symbol.*BitVector ]]; then
        log_error "BV file header format is incorrect. Expected: probe_id<tab>gene_symbol<tab>BitVector"
        exit 1
    fi
    
    # Count probes in BV file
    local probe_count=$(tail -n +2 "$INPUT_BV" | wc -l)
    log_info "📊 Found $probe_count probes in BitVector file"
    
    # Extract sample probes from index file for validation - Updated for probe IDs
    local sample_probes=$(tail -n +4 "$INPUT_IDX" | awk '{print $1}' | grep -E '^cg[0-9]+$' | head -5)
    if [ -z "$sample_probes" ]; then
        log_warning "No valid probe IDs (cg format) found in index file - checking for alternative formats..."
        # Check for any non-empty first column as probe IDs
        sample_probes=$(tail -n +4 "$INPUT_IDX" | awk '{print $1}' | grep -v "^$" | head -5)
        if [ -z "$sample_probes" ]; then
            log_error "No valid probe IDs found in index file"
            exit 1
        fi
    fi
    
    log_success "✅ Environment validation completed"
}

setup_directories() {
    log_info "📁 Setting up directory structure..."
    
    # Create subdirectories for each threshold
    for threshold in "${STAT_THRESHOLDS[@]}"; do
        mkdir -p "$ORIGINAL_DIR/thresh_$threshold"
        mkdir -p "$PERMUTED_DIR/thresh_$threshold"
    done
    
    # Create permutation data directory
    mkdir -p "$PERMUTED_DIR/data"
    
    log_success "✅ Directory structure completed"
}

# ============================================================================
# BITVECTOR PERMUTATION FUNCTIONS
# ============================================================================

create_permuted_bitvector() {
    local input_file=$1
    local output_file=$2
    local perm_id=$3
    local random_seed=$((RANDOM + perm_id * 1000))
    
    log_info "🔄 Creating permuted BitVector #$perm_id (seed: $random_seed)"
    
    # Python script for row-wise BitVector permutation - Updated for probe_id
    python3 -c "
import sys
import pandas as pd
import numpy as np

# Set random seed for reproducibility
np.random.seed($random_seed)

try:
    # Read the BV file
    df = pd.read_csv('$input_file', sep='\t', dtype=str)
    
    # Validate required columns - Updated for probe_id
    required_cols = ['probe_id', 'gene_symbol', 'BitVector']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f'Error: Missing required columns: {missing_cols}', file=sys.stderr)
        sys.exit(1)
    
    # Create permuted dataframe
    df_permuted = df.copy()
    
    # Permute each BitVector independently (row-wise permutation)
    permuted_count = 0
    for idx, row in df.iterrows():
        bitvector = row['BitVector']
        if pd.notna(bitvector) and len(str(bitvector)) > 0:
            # Convert to list of characters
            bv_chars = list(str(bitvector))
            # Shuffle the characters
            np.random.shuffle(bv_chars)
            # Join back to string
            df_permuted.at[idx, 'BitVector'] = ''.join(bv_chars)
            permuted_count += 1
    
    # Save permuted data
    df_permuted.to_csv('$output_file', sep='\t', index=False)
    print(f'Successfully permuted {permuted_count} BitVectors', file=sys.stderr)
    
except Exception as e:
    print(f'Error creating permuted BitVector: {e}', file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)
" 2>>"$LOGS_DIR/permutation_$perm_id.log"

    if [ $? -ne 0 ]; then
        log_error "Failed to create permuted BitVector #$perm_id"
        return 1
    fi
    
    return 0
}

generate_all_permutations() {
    log_info "🎲 Generating $N_PERMUTATIONS permuted datasets..."
    local start_time=$(date +%s)
    
    # Generate permutations sequentially
    for perm_id in $(seq 1 $N_PERMUTATIONS); do
        local output_file="$PERMUTED_DIR/data/permuted_${perm_id}.bv"
        
        if [ ! -f "$output_file" ]; then
            if create_permuted_bitvector "$INPUT_BV" "$output_file" "$perm_id"; then
                show_progress $perm_id $N_PERMUTATIONS "Permutation Generation" $start_time
            else
                log_error "❌ Permutation $perm_id failed"
                return 1
            fi
        else
            log_info "⏭️  Permutation $perm_id already exists"
        fi
    done
    
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    log_success "✅ All permutations generated in ${duration}s"
}

# ============================================================================
# BOOLEANNET PIPELINE FUNCTIONS
# ============================================================================

run_booleannet_pipeline() {
    local input_bv=$1
    local output_prefix=$2
    local stat_threshold=$3
    local log_file=$4
    
    local rl_file="${output_prefix}.rl"
    
    log_info "🔬 Running BooleanNet pipeline: threshold=$stat_threshold"
    
    # Step 1: Generate Boolean network (bitMatrix) - NO TIMEOUT
    log_info "   Step 1/3: bitMatrix - Creating Boolean network..."
    if ! java $JAVA_OPTS -cp "$STEPMINER_JAR" tools.CustomAnalysis \
        boolean bitMatrix "$rl_file" "$input_bv" "$OUTPUT_PH" All \
        $FIXED_ERROR_RATE $stat_threshold $FIXED_SINGLE_THR \
        >>"$log_file" 2>&1; then
        log_error "bitMatrix failed for threshold $stat_threshold"
        return 1
    fi
    
    # Check if RL file was created
    if [ ! -f "$rl_file" ]; then
        log_error "RL file not created: $rl_file"
        return 1
    fi
    
    # Step 2: Fill matrix (bitMatrixFill) - NO TIMEOUT
    log_info "   Step 2/3: bitMatrixFill - Filling matrix..."
    if ! java $JAVA_OPTS -cp "$STEPMINER_JAR" tools.CustomAnalysis \
        boolean bitMatrixFill "$rl_file" \
        >>"$log_file" 2>&1; then
        log_error "bitMatrixFill failed for threshold $stat_threshold"
        return 1
    fi
    
    # Step 3: Generate statistics (bitMatrixFillStats) - NO TIMEOUT
    log_info "   Step 3/3: bitMatrixFillStats - Computing statistics..."
    if ! java $JAVA_OPTS -cp "$STEPMINER_JAR" tools.CustomAnalysis \
        boolean bitMatrixFillStats "$rl_file" \
        >>"$log_file" 2>&1; then
        log_error "bitMatrixFillStats failed for threshold $stat_threshold"
        return 1
    fi
    
    log_success "✅ BooleanNet pipeline completed: $(basename $rl_file)"
    return 0
}

# ============================================================================
# SIMPLIFIED IMPLICATIONS EXTRACTION FUNCTION
# ============================================================================

extract_all_implications() {
    local rl_file=$1
    local idx_file=$2
    local log_file=$3
    
    if [ ! -f "$rl_file" ]; then
        echo "0"
        return 1
    fi
    
    log_info "🔍 Extracting ALL Boolean implications from: $(basename $rl_file)"
    
    # Create detailed implications log
    local implications_log="${log_file%.log}_implications.log"
    echo "=== Boolean Implications Extraction ===" > "$implications_log"
    echo "RL File: $(basename $rl_file)" >> "$implications_log"
    echo "Timestamp: $(date)" >> "$implications_log"
    echo "Command: perl $EXTRACT_SCRIPT all $rl_file $idx_file" >> "$implications_log"
    echo "" >> "$implications_log"
    
    # Start time for extraction
    local extraction_start=$(date +%s)
    
    # Extract ALL implications in one command - NO TIMEOUT
    log_info "   Running extraction command: perl $EXTRACT_SCRIPT all $rl_file $idx_file"
    
    # Run the extraction and capture output
    local extraction_output=$(perl "$EXTRACT_SCRIPT" all "$rl_file" "$idx_file" 2>&1)
    local extraction_status=$?
    
    # Save the raw output
    echo "$extraction_output" >> "$implications_log"
    
    if [ $extraction_status -ne 0 ]; then
        log_error "   Extraction command failed with status $extraction_status"
        echo "0"
        return 1
    fi
    
    # Count the number of implications (lines in output)
    local total_implications=$(echo "$extraction_output" | wc -l)
    
    # Remove empty lines from count
    total_implications=$(echo "$extraction_output" | grep -v "^$" | wc -l)
    
    # Calculate extraction time
    local extraction_time=$(($(date +%s) - extraction_start))
    
    # Log summary
    echo "" >> "$implications_log"
    echo "=== Extraction Summary ===" >> "$implications_log"
    echo "TOTAL RELATIONS EXTRACTED: $total_implications" >> "$implications_log"
    echo "EXTRACTION TIME: ${extraction_time} seconds" >> "$implications_log"
    
    # Simple log output
    log_success "   ${BOLD}EXTRACTED $total_implications RELATIONS${NC} (${extraction_time}s)"
    
    echo "$total_implications"
    return 0
}

# ============================================================================
# ANALYSIS ORCHESTRATION
# ============================================================================

analyze_original_data() {
    log_info "📊 Analyzing original data across all thresholds..."
    local start_time=$(date +%s)
    
    # Create results file header if not exists
    if [ ! -f "$ORIGINAL_DIR/original_results.tsv" ]; then
        echo -e "threshold\timplications_count\tprocessing_time" > "$ORIGINAL_DIR/original_results.tsv"
    fi
    
    local count=0
    local skipped=0
    for threshold in "${STAT_THRESHOLDS[@]}"; do
        count=$((count + 1))
        
        # Check if this threshold was already processed
        if check_existing_results "original" "$threshold"; then
            skipped=$((skipped + 1))
            continue
        fi
        
        show_progress $((count - skipped)) $((${#STAT_THRESHOLDS[@]} - skipped)) "Original Analysis" $start_time
        save_checkpoint "ORIGINAL" "threshold_$threshold"
        
        local output_prefix="$ORIGINAL_DIR/thresh_$threshold/original_thresh_$threshold"
        local log_file="$LOGS_DIR/original_thresh_$threshold.log"
        local analysis_start=$(date +%s)
        
        # Run the 3-step BooleanNet pipeline
        if run_booleannet_pipeline "$INPUT_BV" "$output_prefix" "$threshold" "$log_file"; then
            # Extract ALL implications using the simplified command
            local implications_count=$(extract_all_implications "${output_prefix}.rl" "$INPUT_IDX" "$log_file")
            local analysis_time=$(($(date +%s) - analysis_start))
            
            # Update results file with locking
            safe_append "$ORIGINAL_DIR/original_results.tsv" "$threshold\t$implications_count\t$analysis_time"
            
            log_success "✅ Original analysis completed for threshold $threshold: $implications_count implications"
        else
            safe_append "$ORIGINAL_DIR/original_results.tsv" "$threshold\t0\tFAILED"
            log_error "❌ Original analysis failed for threshold $threshold"
        fi
    done
    
    local total_time=$(($(date +%s) - start_time))
    log_success "✅ Original data analysis completed in ${total_time}s"
}

analyze_permuted_data() {
    log_info "🎲 Analyzing permuted data across all thresholds..."
    local start_time=$(date +%s)
    
    # Create results file header if not exists
    if [ ! -f "$PERMUTED_DIR/permuted_results.tsv" ]; then
        echo -e "threshold\tpermutation_id\timplications_count\tprocessing_time" > "$PERMUTED_DIR/permuted_results.tsv"
    fi
    
    local total_jobs=$((${#STAT_THRESHOLDS[@]} * N_PERMUTATIONS))
    local completed_jobs=0
    local skipped_jobs=0
    
    for threshold in "${STAT_THRESHOLDS[@]}"; do
        log_info "🔬 Processing threshold $threshold with $N_PERMUTATIONS permutations..."
        
        # Process permutations sequentially
        for perm_id in $(seq 1 $N_PERMUTATIONS); do
            # Check if this combination was already processed
            if check_existing_results "permuted" "$threshold" "$perm_id"; then
                skipped_jobs=$((skipped_jobs + 1))
                completed_jobs=$((completed_jobs + 1))
                continue
            fi
            
            save_checkpoint "PERMUTED" "threshold_${threshold}_perm_${perm_id}"
            
            local input_bv="$PERMUTED_DIR/data/permuted_${perm_id}.bv"
            local output_prefix="$PERMUTED_DIR/thresh_$threshold/perm_${perm_id}_thresh_$threshold"
            local log_file="$LOGS_DIR/perm_${perm_id}_thresh_$threshold.log"
            
            if [ ! -f "$input_bv" ]; then
                log_error "❌ Permuted file not found: $(basename $input_bv)"
                safe_append "$PERMUTED_DIR/permuted_results.tsv" "$threshold\t$perm_id\t0\tMISSING_FILE"
                completed_jobs=$((completed_jobs + 1))
                continue
            fi
            
            local analysis_start=$(date +%s)
            
            # Run the 3-step BooleanNet pipeline
            if run_booleannet_pipeline "$input_bv" "$output_prefix" "$threshold" "$log_file"; then
                # Extract ALL implications using the simplified command
                local implications_count=$(extract_all_implications "${output_prefix}.rl" "$INPUT_IDX" "$log_file")
                local analysis_time=$(($(date +%s) - analysis_start))
                
                safe_append "$PERMUTED_DIR/permuted_results.tsv" "$threshold\t$perm_id\t$implications_count\t$analysis_time"
                
                log_success "✅ Permutation $perm_id (threshold $threshold): $implications_count implications"
            else
                safe_append "$PERMUTED_DIR/permuted_results.tsv" "$threshold\t$perm_id\t0\tFAILED"
                log_error "❌ Permutation $perm_id (threshold $threshold): FAILED"
            fi
            
            completed_jobs=$((completed_jobs + 1))
            show_progress $((completed_jobs - skipped_jobs)) $((total_jobs - skipped_jobs)) "Permutation Analysis" $start_time
        done
    done
    
    local total_time=$(($(date +%s) - start_time))
    log_success "✅ Permuted data analysis completed in ${total_time}s"
}

# ============================================================================
# FDR CALCULATION AND REPORTING
# ============================================================================

# Function to extract implications count from log files
get_implications_from_log() {
    local log_file=$1
    if [ -f "$log_file" ]; then
        local count=$(grep "TOTAL RELATIONS EXTRACTED:" "$log_file" | head -1 | awk '{print $4}' | grep -o '[0-9]*')
        if [ ! -z "$count" ]; then
            echo "$count"
        else
            echo "0"
        fi
    else
        echo "0"
    fi
}

calculate_fdr_statistics() {
    log_info "📈 Calculating False Discovery Rate (FDR) statistics..."
    
    # Create comprehensive FDR report
    cat > "$FINAL_REPORT" << 'EOF'
threshold	original_implications	mean_permuted_implications	std_permuted_implications	empirical_fdr	empirical_pvalue	recommended
EOF
    
    # Process each threshold
    for threshold in "${STAT_THRESHOLDS[@]}"; do
        log_info "📊 Processing FDR for threshold: $threshold"
        
        # Get original count from log file first, then TSV as fallback
        local original_log_file="$LOGS_DIR/original_thresh_${threshold}_implications.log"
        local original_count=$(get_implications_from_log "$original_log_file")
        
        # If not found in log, try TSV file
        if [ "$original_count" = "0" ] && [ -f "$ORIGINAL_DIR/original_results.tsv" ]; then
            local tsv_count=$(grep "^$threshold" "$ORIGINAL_DIR/original_results.tsv" | cut -f2)
            if [ ! -z "$tsv_count" ] && [ "$tsv_count" != "0" ] && [ "$tsv_count" != "FAILED" ]; then
                original_count="$tsv_count"
            fi
        fi
        
        if [ "$original_count" = "0" ]; then
            echo -e "$threshold\t0\t0\t0\tNA\tNA\tNO" >> "$FINAL_REPORT"
            log_warning "No original data found for threshold $threshold"
            continue
        fi
        
        # Get permuted counts from log files
        local permuted_counts=""
        for perm_id in $(seq 1 $N_PERMUTATIONS); do
            local perm_log_file="$LOGS_DIR/perm_${perm_id}_thresh_${threshold}_implications.log"
            local count=$(get_implications_from_log "$perm_log_file")
            
            # If not found in log, check if RL file exists (means 0 implications)
            if [ "$count" = "0" ]; then
                local rl_file="$PERMUTED_DIR/thresh_$threshold/perm_${perm_id}_thresh_$threshold.rl"
                if [ -f "$rl_file" ]; then
                    count="0"  # RL file exists but no implications
                else
                    continue  # Skip missing analyses
                fi
            fi
            
            if [ ! -z "$permuted_counts" ]; then
                permuted_counts="${permuted_counts}\n${count}"
            else
                permuted_counts="${count}"
            fi
        done
        
        if [ -z "$permuted_counts" ]; then
            echo -e "$threshold\t$original_count\t0\t0\tNA\tNA\tNO" >> "$FINAL_REPORT"
            log_warning "No permuted data found for threshold $threshold"
            continue
        fi
        
        # Calculate statistics using Python
        python3 -c "
import sys
import numpy as np

# Original count
original_count = float('$original_count')

# Read permuted counts
permuted_data = '''$permuted_counts'''
if permuted_data.strip():
    permuted_counts = []
    for line in permuted_data.strip().split('\n'):
        try:
            val = float(line.strip())
            permuted_counts.append(val)
        except:
            continue
    
    if len(permuted_counts) > 0:
        mean_permuted = np.mean(permuted_counts)
        std_permuted = np.std(permuted_counts, ddof=1)
        
        # Calculate empirical FDR = average(permuted) / original
        empirical_fdr = mean_permuted / original_count if original_count > 0 else float('inf')
        
        # Calculate empirical p-value
        empirical_pvalue = sum(1 for x in permuted_counts if x >= original_count) / len(permuted_counts)
        
        # Recommendation: FDR < 0.05 (5%) and sufficient implications
        recommended = 'YES' if empirical_fdr < 0.05 and original_count >= 100 else 'NO'
        
        print(f'$threshold\t{int(original_count)}\t{mean_permuted:.2f}\t{std_permuted:.2f}\t{empirical_fdr:.4f}\t{empirical_pvalue:.4f}\t{recommended}')
    else:
        print(f'$threshold\t{int(original_count)}\t0\t0\tNA\tNA\tNO')
else:
    print(f'$threshold\t{int(original_count)}\t0\t0\tNA\tNA\tNO')
" >> "$FINAL_REPORT"
    done
    
    log_success "✅ FDR statistics calculated"
}

generate_summary_report() {
    log_info "📋 Generating comprehensive summary report..."
    
    local summary_file="$RESULTS_DIR/analysis_summary.txt"
    
    cat > "$summary_file" << EOF
================================================================================
BooleanNet Permutation Test Analysis Summary
DNA Methylation Probe ID Version
================================================================================
Generated: $(date)

PURPOSE:
Find optimal statistical threshold for Boolean implication analysis of
DNA methylation data through permutation testing and False Discovery Rate 
(FDR) control.

INPUT FILES:
  - BitVector File: $INPUT_BV (probe-based methylation data)
  - Index File: $INPUT_IDX
  - Extract Script: $EXTRACT_SCRIPT
  - StepMiner JAR: $STEPMINER_JAR

PARAMETERS:
  - Statistical Thresholds Tested: ${STAT_THRESHOLDS[*]}
  - Error Rate (fixed): $FIXED_ERROR_RATE
  - Single Threshold (fixed): $FIXED_SINGLE_THR
  - Number of Permutations: $N_PERMUTATIONS
  - Processing: ALL probes (no sampling)
  - Total Analyses: $((${#STAT_THRESHOLDS[@]} * (N_PERMUTATIONS + 1)))

METHODOLOGY:
1. Original Analysis: Run BooleanNet on original methylation data for each threshold
2. Permutation Test: Randomly permute methylation BitVectors and re-analyze
3. FDR Calculation: Compare original vs permuted implications
4. Threshold Selection: Choose threshold with FDR < 0.05

RESULTS SUMMARY:
================================================================================
EOF
    
    # Add FDR results to summary
    echo "" >> "$summary_file"
    cat "$FINAL_REPORT" | column -t -s $'\t' >> "$summary_file"
    
    echo "" >> "$summary_file"
    echo "RECOMMENDED THRESHOLDS (FDR < 0.05):" >> "$summary_file"
    echo "------------------------------------" >> "$summary_file"
    if grep -q "YES" "$FINAL_REPORT"; then
        grep "YES$" "$FINAL_REPORT" | while IFS=$'\t' read -r threshold original mean std fdr pvalue recommended; do
            echo "  ✓ Threshold $threshold: FDR=$fdr, Implications=$original" >> "$summary_file"
        done
    else
        echo "  ✗ No thresholds meet FDR < 0.05 criteria" >> "$summary_file"
    fi
    
    echo "" >> "$summary_file"
    echo "================================================================================" >> "$summary_file"
    echo "OUTPUT FILES:" >> "$summary_file"
    echo "  - FDR Analysis Report: $FINAL_REPORT" >> "$summary_file"
    echo "  - Original Results: $ORIGINAL_DIR/original_results.tsv" >> "$summary_file"
    echo "  - Permuted Results: $PERMUTED_DIR/permuted_results.tsv" >> "$summary_file"
    echo "  - Analysis Logs: $LOGS_DIR/" >> "$summary_file"
    echo "================================================================================" >> "$summary_file"
    
    log_success "✅ Summary report generated: $summary_file"
}

# ============================================================================
# MAIN EXECUTION FUNCTIONS
# ============================================================================

cleanup_on_exit() {
    log_warning "🛑 Process interrupted! Cleaning up..."
    
    # Kill any background Java processes
    pkill -f "stepminer-1.1.jar" 2>/dev/null || true
    
    # Clean up lock files
    rm -rf "$LOCK_DIR" 2>/dev/null || true
    
    log_info "🧹 Cleanup completed"
    exit 1
}

display_final_results() {
    # Show recommended thresholds
    echo -e "${GREEN}🎯 OPTIMAL STATISTICAL THRESHOLDS (FDR < 0.05):${NC}"
    if [ -f "$FINAL_REPORT" ] && grep -q "YES" "$FINAL_REPORT"; then
        grep "YES$" "$FINAL_REPORT" | tail -n +2 | while IFS=$'\t' read -r threshold original mean std fdr pvalue recommended; do
            echo -e "${GREEN}   ✓ Threshold: $threshold | FDR: $fdr | P-value: $pvalue | Implications: $original${NC}"
        done
    else
        echo -e "${YELLOW}   ✗ No thresholds meet FDR < 0.05 criteria${NC}"
        echo -e "${YELLOW}   Consider using more permutations or adjusting parameters${NC}"
    fi
}

main() {
    local script_start=$(date +%s)
    
    echo -e "${MAGENTA}================================================================================${NC}"
    echo -e "${MAGENTA}🧬 BooleanNet Permutation Test Framework v8.0 - Probe ID Version${NC}"
    echo -e "${MAGENTA}     DNA Methylation Data Analysis${NC}"
    echo -e "${MAGENTA}================================================================================${NC}"
    
    log_info "🚀 Starting Boolean implication analysis for methylation data with permutation testing..."
    log_info "📊 Statistical thresholds to test: ${STAT_THRESHOLDS[*]}"
    log_info "🎲 Permutations per threshold: $N_PERMUTATIONS"
    log_info "🧬 Processing: ALL methylation probes (no sampling)"
    log_info "💻 Using $CORES cores (${TOTAL_CORES} available)"
    
    # Setup and validation
    validate_environment
    setup_directories
    
    # Generate permuted datasets
    generate_all_permutations
    
    # Analyze original data
    analyze_original_data
    
    # Analyze permuted data
    analyze_permuted_data
    
    # Calculate FDR and generate reports
    calculate_fdr_statistics
    generate_summary_report
    
    local total_time=$(($(date +%s) - script_start))
    local hours=$((total_time / 3600))
    local minutes=$(( (total_time % 3600) / 60 ))
    local seconds=$((total_time % 60))
    
    echo -e "${MAGENTA}================================================================================${NC}"
    log_success "🎉 ANALYSIS COMPLETED SUCCESSFULLY!"
    log_success "⏱️  Total execution time: ${hours}h ${minutes}m ${seconds}s"
    log_success "📁 Results directory: $RESULTS_DIR"
    log_success "📋 Summary report: $RESULTS_DIR/analysis_summary.txt"
    log_success "📈 FDR report: $FINAL_REPORT"
    echo -e "${MAGENTA}================================================================================${NC}"
    
    # Display final results
    display_final_results
    
    echo -e "${MAGENTA}================================================================================${NC}"
    echo -e "${GREEN}🎊 BooleanNet Permutation Test for methylation data completed successfully!${NC}"
}

# Export functions for potential use
export -f log_info log_success log_warning log_error
export -f run_booleannet_pipeline extract_all_implications create_permuted_bitvector
export -f get_implications_from_log acquire_lock release_lock safe_append
export JAVA_OPTS STEPMINER_JAR OUTPUT_PH FIXED_ERROR_RATE FIXED_SINGLE_THR
export INPUT_IDX EXTRACT_SCRIPT PERMUTED_DIR LOGS_DIR LOCK_DIR
export GREEN BLUE YELLOW RED CYAN MAGENTA BOLD NC

# Set trap for cleanup
trap cleanup_on_exit INT TERM

# Run main function
main

# Script completed
exit 0
