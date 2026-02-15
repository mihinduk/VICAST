#!/bin/bash
# =============================================================================
# VICAST-Analyze: Full Workflow (All Steps 1-9)
# =============================================================================
# Runs complete pipeline from QC through annotation without stopping
#
# Prerequisites:
#   - SNPEFF_HOME or SNPEFF_JAR environment variable set
#   - Conda environment with required tools (bwa, samtools, lofreq, fastp)
#   - Target genome added to SnpEff database via VICAST-annotate
# =============================================================================

set -e  # Exit on error

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [--large-files|--extremely-large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo ""
    echo "This script runs the complete workflow (Steps 1-9):"
    echo "  1. Prepare reference genome"
    echo "  2. Calculate read statistics"
    echo "  3. Clean reads (fastp QC)"
    echo "  4. Map reads and call variants"
    echo "  5. Generate depth file"
    echo "  6. Run diagnostic report"
    echo "  7. Filter variants"
    echo "  8. Annotate variants with snpEff"
    echo "  9. Parse annotations"
    echo ""
    echo "For two-part workflow with manual review:"
    echo "  Part 1: run_vicast_analyze_qc_only.sh <R1> <R2> <accession>"
    echo "  Part 2: run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>"
    echo ""
    echo "Prerequisites:"
    echo "  - snpEff database must exist (use VICAST-annotate to add genomes)"
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
    echo "Loading default configuration from: ${VICAST_HOME}/vicast_config.template.sh"
    source "${VICAST_HOME}/vicast_config.template.sh"
fi

# Source local pipeline config if it exists (for backwards compatibility)
if [ -f "${SCRIPT_DIR}/pipeline_config.sh" ]; then
    echo "Loading pipeline configuration from: ${SCRIPT_DIR}/pipeline_config.sh"
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
        echo "==============================="
        echo ""
        echo "Please set one of the following environment variables:"
        echo "  export SNPEFF_HOME=/path/to/snpEff"
        echo "  export SNPEFF_JAR=/path/to/snpEff/snpEff.jar"
        echo ""
        exit 1
    fi
else
    SNPEFF_DIR="$(dirname "$SNPEFF_JAR")"
fi

# Set SNPEFF_DATA if not set
if [ -z "$SNPEFF_DATA" ]; then
    SNPEFF_DATA="${SNPEFF_DIR}/data"
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
    # In Docker with micromamba, use base environment
    VICAST_CONDA_ENV="base"
fi

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
        # In Docker, tools are already in base environment - no activation needed
        echo "Detected micromamba (Docker environment) - using current environment"
    fi
else
    echo "Warning: conda/micromamba not found - assuming tools are in PATH"
fi

# Check for database files directly (works offline, faster)
if [ -f "$SNPEFF_DATA/$ACCESSION/snpEffectPredictor.bin" ]; then
    echo "✓ Database found: $SNPEFF_DATA/$ACCESSION/"
    echo "✓ Database is built and ready"
    echo ""
else
    echo ""
    echo "❌ ERROR: Genome $ACCESSION database not found!"
    echo ""
    echo "Checked location: $SNPEFF_DATA/$ACCESSION/"
    echo ""
    echo "Please install the database first:"
    echo "  Option 1 (pre-built database):"
    echo "    install_prebuilt_database.sh --install $ACCESSION"
    echo ""
    echo "  Option 2 (custom annotation):"
    echo "    python3 /opt/vicast/vicast-annotate/step1_parse_viral_genome.py $ACCESSION"
    echo "    # Edit the TSV file, then:"
    echo "    python3 /opt/vicast/vicast-annotate/step2_add_to_snpeff.py $ACCESSION ${ACCESSION}.tsv"
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