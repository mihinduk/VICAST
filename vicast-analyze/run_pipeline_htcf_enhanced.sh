#!/bin/bash
# Enhanced Viral Pipeline Wrapper
# This script runs the consolidated pipeline and generates next_steps commands
# The consolidated pipeline handles all virus configuration and database setup

# Function to display usage
usage() {
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <threads> [--large-files]"
    echo ""
    echo "Arguments:"
    echo "  R1_fastq      : Path to R1 FASTQ file"
    echo "  R2_fastq      : Path to R2 FASTQ file"
    echo "  accession     : GenBank accession number (e.g., NC_001477.1)"
    echo "  threads       : Number of threads to use"
    echo "  --large-files : Optional flag for large file handling"
    echo ""
    echo "Prerequisites:"
    echo "  Activate conda before running: source /ref/sahlab/software/anaconda3/bin/activate"
    echo ""
    echo "Configuration:"
    echo "  Create pipeline_config.sh from pipeline_config.template.sh to customize paths"
    echo ""
    echo "This enhanced version automatically:"
    echo "  - Downloads virus reference and annotation"
    echo "  - Adds new viruses to configuration database"
    echo "  - Verifies proper gene annotation exists"
    echo "  - Runs complete assembly and variant calling"
    exit 1
}

# Check if correct number of arguments provided
if [ $# -lt 4 ]; then
    usage
fi

# Assign arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=$4
LARGE_FILES_FLAG=""

# Check for optional flags
if [ $# -eq 5 ] && [ "$5" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
fi

# Extract sample name from R1 filename
SAMPLE_NAME=$(basename "$R1" | sed -E "s/_R?[12](_[0-9]+)?..*$//")

# Get the directory where this script is located (auto-detection)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Default PIPELINE_BASE to script directory (works for new users out-of-the-box)
PIPELINE_BASE="${SCRIPT_DIR}"

# Default MAMBA_CMD (can be overridden by config)
MAMBA_CMD="conda run -n viral_genomics"

# Source configuration file if it exists (can override defaults)
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"
if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from: $CONFIG_FILE"
    source "$CONFIG_FILE"
    # Config can override PIPELINE_BASE and MAMBA_CMD if needed
else
    echo "ℹ️  No configuration file found, using auto-detected paths"
    echo "   Pipeline base: $PIPELINE_BASE"
fi

# Verify critical paths exist, prompt if missing
CONSOLIDATED_PIPELINE="${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh"
if [ ! -f "$CONSOLIDATED_PIPELINE" ]; then
    echo ""
    echo "❌ Error: Consolidated pipeline not found at:"
    echo "   $CONSOLIDATED_PIPELINE"
    echo ""
    read -p "Enter the full path to your VICAST/vicast-analyze directory: " USER_PIPELINE_BASE

    if [ -d "$USER_PIPELINE_BASE" ] && [ -f "$USER_PIPELINE_BASE/run_pipeline_htcf_consolidated.sh" ]; then
        PIPELINE_BASE="$USER_PIPELINE_BASE"
        CONSOLIDATED_PIPELINE="${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh"
        echo "✅ Found consolidated pipeline at: $CONSOLIDATED_PIPELINE"
        echo ""
        echo "💡 Tip: Save this path in pipeline_config.sh to avoid this prompt:"
        echo "   PIPELINE_BASE=\"$PIPELINE_BASE\""
    else
        echo "❌ Error: Consolidated pipeline still not found. Exiting."
        exit 1
    fi
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found in PATH"
    echo ""
    echo "Please activate conda first:"
    echo "  source /ref/sahlab/software/anaconda3/bin/activate"
    echo ""
    echo "Then re-run this script."
    exit 1
fi

# Set up HTCF environment
echo "========================================="
echo "Enhanced Viral Pipeline - Module 1"
echo "========================================="
echo "Sample: $SAMPLE_NAME"
echo "Accession: $ACCESSION"
echo "Working directory: $(pwd)"
echo "Time: $(date)"
echo "========================================="

# Run the consolidated pipeline
# The consolidated pipeline handles:
# - snpEff database checking and building
# - Reference genome download from NCBI if needed
# - Complete assembly and variant calling workflow
echo ""
echo "Running viral assembly pipeline..."
echo "========================================="

bash "${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh" "$R1" "$R2" "$ACCESSION" "$THREADS" $LARGE_FILES_FLAG

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "✅ Pipeline completed successfully!"
    echo "========================================="
    
    # Create next steps file
    cat > "next_steps_${SAMPLE_NAME}.txt" << EOF
# Next steps for ${SAMPLE_NAME} (${ACCESSION})

## Pipeline already completed:#   - QC and trimming (fastp)#   - Alignment to reference (BWA)#   - Variant calling (LoFreq)#   - Variant annotation (SnpEff)#   - Depth file generation (output: ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt)## Optional post-processing steps:

# 1. Parse mutations for filtering
${MAMBA_CMD} \\
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/parse_snpeff_tsv.py \\
  "cleaned_seqs/variants/${SAMPLE_NAME}.snpEFF.ann.tsv" \\
  "${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv" \\
  --quality 1000 --depth 200 --freq 0.01


# 2. Diagnostic Report (Optional)
sbatch ${PIPELINE_BASE}/viral_pipeline/analysis/submit_viral_diagnostic.sh \\
  "${R1}" \\
  "${R2}" \\
  ${ACCESSION} \\
  diagnostic_${SAMPLE_NAME} \\
  4

# 3. Generate consensus
${MAMBA_CMD} \\
  python3 ${PIPELINE_BASE}/viral_pipeline/utils/generate_filtered_consensus.py \\
  --vcf ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \\
  --reference cleaned_seqs/${ACCESSION}.fasta \\
  --accession ${ACCESSION} \\
  --quality 1000 --depth 200 --freq 0.50 \\
  --output-prefix ${SAMPLE_NAME}_results/${SAMPLE_NAME}_consensus
EOF
    
    echo ""
    echo "📝 Visualization commands saved to: next_steps_${SAMPLE_NAME}.txt"
    echo ""
else
    echo ""
    echo "❌ Pipeline failed with exit code: $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi