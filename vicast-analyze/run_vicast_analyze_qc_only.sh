#!/bin/bash
# VICAST-Analyze: QC Workflow Only (Steps 1-6)
# Runs: mapping, variant calling, depth file, and diagnostic report
# STOPS before annotation to allow user review

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [--large-files|--extremely-large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo ""
    echo "This script runs QC workflow (Steps 1-6):"
    echo "  1. Prepare reference genome"
    echo "  2. Calculate read statistics"
    echo "  3. Clean reads (fastp QC)"
    echo "  4. Map reads and call variants"
    echo "  5. Generate depth file"
    echo "  6. Run diagnostic report"
    echo ""
    echo "After completion, review QC outputs and run annotation if satisfied:"
    echo "  ./run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>"
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
    echo "Loading configuration from: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "⚠️  No configuration file found. Using hardcoded paths."
    MAMBA_CMD="conda run -n viral_genomics"
    SNPEFF_DIR="/ref/sahlab/software/snpEff"
    SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
    JAVA_PATH="java"
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found in PATH"
    echo ""
    echo "Please activate conda first:"
    echo "  source /ref/sahlab/software/anaconda3/bin/activate"
    exit 1
fi

# Set paths
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

echo "========================================"
echo "VICAST-ANALYZE: QC WORKFLOW (Steps 1-6)"
echo "========================================"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Accession: $ACCESSION"
echo "  Threads: $THREADS"
echo "  Memory mode: $(if [ "$LARGE_FILES_FLAG" == "--extremely-large-files" ]; then echo "extreme (256GB)"; elif [ "$LARGE_FILES_FLAG" == "--large-files" ]; then echo "large (64GB)"; else echo "standard (32GB)"; fi)"
echo "========================================"
echo ""

# Pre-flight check: Verify snpEff database exists
echo "Pre-flight check: snpEff database"
echo "Checking if $ACCESSION is in snpEff database..."

eval "$(conda shell.bash hook)"
conda activate vicast_analyze

if java -jar "$SNPEFF_JAR" databases 2>/dev/null | grep -q "$ACCESSION"; then
    echo "✓ snpEff database found for $ACCESSION"
    if [ -f "$SNPEFF_DIR/data/$ACCESSION/snpEffectPredictor.bin" ]; then
        echo "✓ Database is built and ready"
    else
        echo ""
        echo "❌ ERROR: Database exists in config but is not built!"
        echo "Please use VICAST-annotate to build the database:"
        echo "  cd /path/to/VICAST/vicast-annotate"
        echo "  python3 vicast_annotate.py $ACCESSION"
        exit 1
    fi
else
    echo ""
    echo "❌ ERROR: Genome $ACCESSION not found in snpEff database!"
    echo "Please use VICAST-annotate to add this genome first."
    exit 1
fi

echo ""
echo "Running QC workflow (Steps 1-6)..."
echo ""

# Set PYTHONUNBUFFERED for real-time output
export PYTHONUNBUFFERED=1

# Run pipeline with --skip-annotation flag
if [ -n "$LARGE_FILES_FLAG" ]; then
    stdbuf -oL -eL python -u ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH" \
        --skip-annotation \
        $LARGE_FILES_FLAG
else
    stdbuf -oL -eL python -u ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH" \
        --skip-annotation
fi

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    # Extract sample name
    SAMPLE_NAME=$(basename "$R1" | sed 's/_R1\\.fastq\\.gz$//' | sed 's/_R1\\.qc\\.fastq\\.gz$//')

    echo ""
    echo "========================================"
    echo "✓ QC WORKFLOW COMPLETED (Steps 1-6)"
    echo "========================================"
    echo ""
    echo "REVIEW THE FOLLOWING QC OUTPUTS:"
    echo ""
    echo "1. Depth File:"
    echo "   ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt"
    echo ""
    echo "2. Diagnostic Report:"
    echo "   diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt"
    echo "   diagnostic_${SAMPLE_NAME}/diagnostic_${SAMPLE_NAME}_presentation_ready_report.html"
    echo ""
    echo "3. VCF Files:"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*.vcf" 2>/dev/null | head -5
    echo ""
    echo "NEXT STEPS:"
    echo ""
    echo "  If QC results are satisfactory, continue with annotation:"
    echo "    ./run_vicast_analyze_annotate_only.sh $R1 $R2 $ACCESSION"
    echo ""
    echo "  Or submit as SLURM job:"
    echo "    sbatch --wrap=\"./run_vicast_analyze_annotate_only.sh $R1 $R2 $ACCESSION\""
    echo ""
    echo "========================================"
else
    echo ""
    echo "❌ QC workflow failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi
