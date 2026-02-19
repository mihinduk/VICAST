#!/bin/bash
# =============================================================================
# VICAST-Analyze: Annotation Workflow Only (Steps 7-9)
# =============================================================================
# Runs: variant filtering, snpEff annotation, and parsing
# REQUIRES: QC workflow must be run first
# Uses --resume-from-vcf for true efficiency - no work is duplicated!
#
# Prerequisites:
#   - SNPEFF_HOME or SNPEFF_JAR environment variable set
#   - QC workflow completed (run_vicast_analyze_qc_only.sh)
# =============================================================================

set -e  # Exit on error

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [OPTIONS]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo ""
    echo "PREREQUISITES:"
    echo "  Run QC workflow first: run_vicast_analyze_qc_only.sh <R1> <R2> <accession>"
    echo ""
    echo "This script runs annotation workflow (Steps 7-9) EFFICIENTLY:"
    echo "  Discovers existing VCF files from previous QC run (--resume-from-vcf)"
    echo "  Skips Steps 1-6 completely - no duplication of work!"
    echo "  7. Filter variants (lofreq filter with depth + quality thresholds)"
    echo "  8. Annotate variants with snpEff"
    echo "  9. Parse annotations"
    echo ""
    echo "Options:"
    echo "  --min-depth N    Minimum read depth (default: 200)"
    echo "  --min-qual N     Minimum variant quality, phred (default: 90)"
    echo "  --large-files    High-memory mode for large files (1-5GB)"
    echo "  --extremely-large-files  Extreme memory mode for >5GB files"
    echo ""
    echo "Environment variables (set before running):"
    echo "  SNPEFF_HOME     - Path to SnpEff installation"
    echo "  SNPEFF_JAR      - Path to snpEff.jar (optional if SNPEFF_HOME set)"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=4
LARGE_FILES_FLAG=""
MIN_DEPTH=200
MIN_QUAL=90

# Parse optional flags from all remaining arguments
shift 3  # Remove R1, R2, ACCESSION
while [ $# -gt 0 ]; do
    case "$1" in
        --large-files)
            LARGE_FILES_FLAG="--large-files"
            ;;
        --extremely-large-files)
            LARGE_FILES_FLAG="--extremely-large-files"
            ;;
        --min-depth)
            shift
            MIN_DEPTH="$1"
            ;;
        --min-qual)
            shift
            MIN_QUAL="$1"
            ;;
        *)
            # Treat first bare number as threads
            if [[ "$1" =~ ^[0-9]+$ ]] && [ "$THREADS" -eq 4 ]; then
                THREADS="$1"
            else
                echo "Unknown option: $1"
                exit 1
            fi
            ;;
    esac
    shift
done

# =============================================================================
# Configuration Loading
# =============================================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VICAST_HOME="$( cd "${SCRIPT_DIR}/.." && pwd )"

# Source user configuration if it exists
if [ -f "$HOME/.vicast/config.sh" ]; then
    echo "Loading user configuration from: $HOME/.vicast/config.sh"
    source "$HOME/.vicast/config.sh"
elif [ -f "${VICAST_HOME}/vicast_config.template.sh" ]; then
    source "${VICAST_HOME}/vicast_config.template.sh"
fi

# Source local pipeline config if it exists
if [ -f "${SCRIPT_DIR}/pipeline_config.sh" ]; then
    source "${SCRIPT_DIR}/pipeline_config.sh"
fi

# =============================================================================
# Environment Detection
# =============================================================================

# Detect SnpEff paths
if [ -z "$SNPEFF_JAR" ]; then
    if [ -n "$SNPEFF_HOME" ]; then
        SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
        SNPEFF_DIR="${SNPEFF_HOME}"
    elif [ -n "$SNPEFF_DIR" ]; then
        SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
    else
        echo ""
        echo "ERROR: SnpEff not configured"
        echo "Please set SNPEFF_HOME or SNPEFF_JAR environment variable"
        exit 1
    fi
else
    SNPEFF_DIR="$(dirname "$SNPEFF_JAR")"
fi

# Detect Java
if [ -z "$JAVA_PATH" ]; then
    if [ -n "$JAVA_HOME" ]; then
        JAVA_PATH="${JAVA_HOME}/bin/java"
    else
        JAVA_PATH="java"
    fi
fi

# Detect conda environment
VICAST_CONDA_ENV="${VICAST_CONDA_ENV:-vicast_analyze}"

# Check if conda/micromamba is available (optional in Docker)
CONDA_AVAILABLE=false
if command -v conda &> /dev/null; then
    CONDA_AVAILABLE=true
    CONDA_CMD="conda"
elif command -v micromamba &> /dev/null; then
    CONDA_AVAILABLE=true
    CONDA_CMD="micromamba"
    VICAST_CONDA_ENV="base"
fi

PIPELINE_DIR="$(cd "$(dirname "$0")"; pwd)"

# Validate input files
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "============================================="
echo "VICAST-ANALYZE: ANNOTATION WORKFLOW (Steps 7-9)"
echo "============================================="
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Accession: $ACCESSION"
echo "  Threads: $THREADS"
echo "  Min depth: $MIN_DEPTH"
echo "  Min qual:  $MIN_QUAL (phred)"
echo "============================================="
echo ""

# Verify QC workflow outputs exist
echo "Verifying QC workflow outputs..."
SAMPLE_NAME=$(basename "$R1" | sed -E "s/_R?[12](_[0-9]+)?..*$//")

if [ ! -d "./cleaned_seqs/variants" ]; then
    echo ""
    echo "❌ ERROR: variants directory not found!"
    echo "Please run QC workflow first:"
    echo "  run_vicast_analyze_qc_only.sh $R1 $R2 $ACCESSION"
    exit 1
fi

# Check for VCF files
VCF_COUNT=$(find ./cleaned_seqs/variants -name "*.vcf" 2>/dev/null | wc -l)
if [ "$VCF_COUNT" -eq 0 ]; then
    echo ""
    echo "❌ ERROR: No VCF files found in ./cleaned_seqs/variants/"
    echo "Please run QC workflow first:"
    echo "  run_vicast_analyze_qc_only.sh $R1 $R2 $ACCESSION"
    exit 1
fi

echo "✓ Found $VCF_COUNT VCF file(s) from QC workflow"
echo ""

# Activate conda/micromamba environment if available
if [ "$CONDA_AVAILABLE" = true ]; then
    if [ "$CONDA_CMD" = "conda" ]; then
        eval "$(conda shell.bash hook)"
        if conda activate "$VICAST_CONDA_ENV" 2>/dev/null; then
            echo "Activated conda environment: $VICAST_CONDA_ENV"
        else
            echo "Warning: Could not activate conda environment '$VICAST_CONDA_ENV'"
            echo "Attempting to continue with current environment..."
        fi
    elif [ "$CONDA_CMD" = "micromamba" ]; then
        echo "Detected micromamba (Docker environment) - using current environment"
    fi
else
    echo "Warning: conda/micromamba not found - assuming tools are in PATH"
fi

# Set PYTHONUNBUFFERED for real-time output
export PYTHONUNBUFFERED=1

echo "Running annotation workflow (Steps 7-9) ONLY..."
echo ""

# NEW: Using --resume-from-vcf flag to skip Steps 1-6 completely
# This discovers existing VCF files from previous QC run and jumps straight to annotation.
# No more re-running or re-validating - true efficiency!

echo "✓ Using new --resume-from-vcf mode for maximum efficiency"
echo "✓ Discovering existing VCF files from previous QC run..."
echo "✓ Running Steps 7-9 only (filter, annotate, parse)..."
echo ""

PIPELINE_CMD="stdbuf -oL -eL python -u ${PIPELINE_DIR}/viral_pipeline.py \
    --r1 \"$R1\" \
    --r2 \"$R2\" \
    --accession \"$ACCESSION\" \
    --threads $THREADS \
    --snpeff-jar \"$SNPEFF_JAR\" \
    --java-path \"$JAVA_PATH\" \
    --min-depth $MIN_DEPTH \
    --min-qual $MIN_QUAL \
    --resume-from-vcf"

if [ -n "$LARGE_FILES_FLAG" ]; then
    PIPELINE_CMD="$PIPELINE_CMD $LARGE_FILES_FLAG"
fi

eval $PIPELINE_CMD

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "============================================="
    echo "✓ ANNOTATION WORKFLOW COMPLETED (Steps 7-9)"
    echo "============================================="
    echo ""
    echo "OUTPUT FILES:"
    echo "  Annotation files:"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*.snpEFF.ann.tsv" 2>/dev/null | head -5 | sed 's/^/    /'
    echo "  Parsed files (min depth 200):"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*_200.tsv" 2>/dev/null | head -5 | sed 's/^/    /'
    echo "  Summary files:"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*summary*" 2>/dev/null | head -5 | sed 's/^/    /'

    # Show all annotation files if sample name match fails
    if [ -z "$(find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*.snpEFF.ann.tsv" 2>/dev/null)" ]; then
        echo "  All recent annotation files:"
        find ./cleaned_seqs/variants -name "*.snpEFF.ann.tsv" -newer ${0} 2>/dev/null | head -5 | sed 's/^/    /'
    fi
    echo ""
    echo "Pipeline completed successfully!"
    echo "============================================="
else
    echo ""
    echo "❌ Annotation workflow failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi
