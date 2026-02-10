#!/bin/bash
# Consolidated viral pipeline for HTCF
# Requires snpEff database to be pre-configured via VICAST-annotate

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [--large-files|--extremely-large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --large-files"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --extremely-large-files"
    echo ""
    echo "Prerequisites:"
    echo "  Activate conda before running (if not already active)"
    echo "  Or set CONDA_BASE environment variable"
    echo ""
    echo "Configuration:"
    echo "  Create pipeline_config.sh from pipeline_config.template.sh to customize paths"
    echo ""
    echo "This script will:"
    echo "  1. Verify snpEff database exists for the accession (must be pre-configured)"
    echo "  2. Run the complete viral genomics pipeline"
    echo ""
    echo "Note: Use VICAST-annotate to add new genomes to snpEff database:"
    echo "  python3 vicast_annotate.py <accession>"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=${4:-4}  # Default to 4 threads if not specified
LARGE_FILES_FLAG=""

# Check if memory flags are provided (can be in position 4 or 5)
if [ $# -ge 4 ] && [ "$4" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
    THREADS=4  # Reset to default since flag was in threads position
elif [ $# -ge 4 ] && [ "$4" == "--extremely-large-files" ]; then
    LARGE_FILES_FLAG="--extremely-large-files"
    THREADS=4  # Reset to default since flag was in threads position
elif [ $# -ge 5 ] && [ "$5" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
elif [ $# -ge 5 ] && [ "$5" == "--extremely-large-files" ]; then
    LARGE_FILES_FLAG="--extremely-large-files"
fi

# Source configuration file if it exists, with fallbacks
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"

if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "⚠️  No configuration file found. Using environment variables or defaults."
    # Fallback defaults - use environment variables if set
    MAMBA_CMD="${MAMBA_CMD:-conda run -n viral_genomics}"
    SNPEFF_DIR="${SNPEFF_DIR:-}"
    if [ -z "$SNPEFF_DIR" ] && [ -d "/ref/sahlab/software/snpEff" ]; then
    SNPEFF_DIR="${SNPEFF_DIR:-/ref/sahlab/software/snpEff}"  # Use env var or default
    fi
    SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
    JAVA_PATH="${JAVA_PATH:-java}"
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found in PATH"
    echo ""
    echo "Please activate conda first:"
    if [ -n "$CONDA_ACTIVATE" ]; then
        echo "  $CONDA_ACTIVATE"
    else
        echo "  source <path-to-conda>/bin/activate"
        echo "  Or set CONDA_BASE environment variable"
    fi
    echo ""
    echo "Then re-run this script."
    exit 1
fi

# Set up environment
echo "Setting up HTCF environment..."

# Set paths
PIPELINE_DIR="$(cd "$(dirname "$0")"; pwd)"

# Validate input files exist
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "Running viral pipeline with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Accession: $ACCESSION"
echo "  Threads: $THREADS"
echo "  Memory mode: $(if [ "$LARGE_FILES_FLAG" == "--extremely-large-files" ]; then echo "extreme (256GB)"; elif [ "$LARGE_FILES_FLAG" == "--large-files" ]; then echo "large (64GB)"; else echo "standard (32GB)"; fi)"
echo "  Working directory: $(pwd)"

# =============================================================================
# PRE-FLIGHT CHECK: Verify snpEff database exists
# =============================================================================
# This pipeline requires snpEff database to be pre-configured via VICAST-annotate
# Database building is intentionally NOT done here to maintain separation of concerns

echo "========================================="
echo "Pre-flight check: snpEff database"
echo "========================================="
echo "Checking if $ACCESSION is in snpEff database..."

# Activate conda environment for snpEff check
eval "$(conda shell.bash hook)"
conda activate vicast_analyze

if java -jar "$SNPEFF_JAR" databases 2>/dev/null | grep -q "$ACCESSION"; then
    echo "✓ snpEff database found for $ACCESSION in config"

    # Verify the database is actually built (has snpEffectPredictor.bin)
    if [ -f "$SNPEFF_DIR/data/$ACCESSION/snpEffectPredictor.bin" ]; then
        echo "✓ Database is built and ready"
        echo ""
    else
        echo ""
        echo "❌ ERROR: Database exists in config but is not built!"
        echo ""
        echo "Please use VICAST-annotate to properly build the database:"
        echo "  cd /path/to/VICAST/vicast-annotate"
        echo "  python3 vicast_annotate.py $ACCESSION"
        echo ""
        exit 1
    fi
else
    echo ""
    echo "❌ ERROR: Genome $ACCESSION not found in snpEff database!"
    echo ""
    echo "This pipeline requires the snpEff database to be set up first."
    echo "Please use VICAST-annotate to add this genome:"
    echo ""
    echo "  cd /path/to/VICAST/vicast-annotate"
    echo "  python3 vicast_annotate.py $ACCESSION"
    echo ""
    echo "This will:"
    echo "  - Download the reference genome and GenBank annotation"
    echo "  - Add the genome to snpEff configuration"
    echo "  - Build the snpEff database with viral-friendly settings"
    echo ""
    exit 1
fi

echo "Running the viral pipeline..."

# Conda environment already activated from pre-flight check

# Set PYTHONUNBUFFERED for real-time output visibility
# This is critical for monitoring progress in SLURM jobs
export PYTHONUNBUFFERED=1

# Run the pipeline with unbuffered I/O
# stdbuf forces line-buffered output, python -u disables Python buffering
echo "Starting pipeline with real-time progress indicators..."
if [ -n "$LARGE_FILES_FLAG" ]; then
    echo "Using $LARGE_FILES_FLAG for increased memory allocation"
    stdbuf -oL -eL python -u ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH" \
        --large-files
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
    echo "Pipeline completed successfully!"
    echo ""
    
    # Extract sample name from R1 filename to only show files for current sample
    SAMPLE_NAME=$(basename "$R1" | sed -E "s/_R?[12](_[0-9]+)?..*$//")
    
    echo "Output files:"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*.tsv" -o -name "*${SAMPLE_NAME}*.vcf" 2>/dev/null | sort
    echo ""
    echo "Key results:"
    if ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv >/dev/null 2>&1; then
        echo "  Annotation TSV: $(ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv | head -1)"
        echo "  Total annotated variants: $(wc -l < $(ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv | head -1))"
    fi
else
    echo "Pipeline failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi