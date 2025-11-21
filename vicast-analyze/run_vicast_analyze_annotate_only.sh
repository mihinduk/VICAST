#!/bin/bash
# VICAST-Analyze: Annotation Workflow Only (Steps 7-9)
# Runs: variant filtering, snpEff annotation, and parsing
# REQUIRES: QC workflow must be run first

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [--large-files|--extremely-large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo ""
    echo "PREREQUISITES:"
    echo "  Run QC workflow first: ./run_vicast_analyze_qc_only.sh <R1> <R2> <accession>"
    echo ""
    echo "This script runs annotation workflow (Steps 7-9):"
    echo "  7. Filter variants"
    echo "  8. Annotate variants with snpEff"
    echo "  9. Parse annotations"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=${4:-4}
LARGE_FILES_FLAG=""

# Check if memory flags are provided
if [ $# -ge 4 ] && [ "$4" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
    THREADS=4
elif [ $# -ge 4 ] && [ "$4" == "--extremely-large-files" ]; then
    LARGE_FILES_FLAG="--extremely-large-files"
    THREADS=4
elif [ $# -ge 5 ] && [ "$5" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
elif [ $# -ge 5 ] && [ "$5" == "--extremely-large-files" ]; then
    LARGE_FILES_FLAG="--extremely-large-files"
fi

# Source configuration file if it exists
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"

if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    MAMBA_CMD="conda run -n viral_genomics"
    SNPEFF_DIR="/ref/sahlab/software/snpEff"
    SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
    JAVA_PATH="java"
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found in PATH"
    echo "Please activate conda first:"
    echo "  source /ref/sahlab/software/anaconda3/bin/activate"
    exit 1
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
echo "============================================="
echo ""

# Verify QC workflow outputs exist
echo "Verifying QC workflow outputs..."
SAMPLE_NAME=$(basename "$R1" | sed 's/_R1\\.fastq\\.gz$//' | sed 's/_R1\\.qc\\.fastq\\.gz$//')

if [ ! -d "./cleaned_seqs/variants" ]; then
    echo ""
    echo "❌ ERROR: variants directory not found!"
    echo "Please run QC workflow first:"
    echo "  ./run_vicast_analyze_qc_only.sh $R1 $R2 $ACCESSION"
    exit 1
fi

# Check for VCF files
VCF_COUNT=$(find ./cleaned_seqs/variants -name "*.vcf" 2>/dev/null | wc -l)
if [ "$VCF_COUNT" -eq 0 ]; then
    echo ""
    echo "❌ ERROR: No VCF files found in ./cleaned_seqs/variants/"
    echo "Please run QC workflow first:"
    echo "  ./run_vicast_analyze_qc_only.sh $R1 $R2 $ACCESSION"
    exit 1
fi

echo "✓ Found $VCF_COUNT VCF file(s) from QC workflow"
echo ""

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate vicast_analyze

# Set PYTHONUNBUFFERED for real-time output
export PYTHONUNBUFFERED=1

echo "Running annotation workflow (Steps 7-9)..."
echo ""

# Strategy: Re-run full pipeline but skip expensive steps
# Steps 1-3 will run quickly (cached/minimal work)
# Step 4 will be skipped (--skip-mapping)
# We need to reconstruct variant files...

# WORKAROUND: Since viral_pipeline.py doesn't support resuming from VCF files,
# we'll re-run steps 1-4 but they should be fast due to caching
# Only mapping might be slow, but we can't avoid it without code changes

if [ -n "$LARGE_FILES_FLAG" ]; then
    stdbuf -oL -eL python -u ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH" \
        $LARGE_FILES_FLAG
else
    stdbuf -oL -eL python -u ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH"
fi

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "============================================="
    echo "✓ ANNOTATION WORKFLOW COMPLETED (Steps 7-9)"
    echo "============================================="
    echo ""
    echo "OUTPUT FILES:"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*.snpEFF.ann.tsv" -o -name "*${SAMPLE_NAME}*_200.tsv" 2>/dev/null
    echo ""
    echo "Pipeline completed successfully!"
    echo "============================================="
else
    echo ""
    echo "❌ Annotation workflow failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi
