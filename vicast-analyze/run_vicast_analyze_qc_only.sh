#!/bin/bash
# =============================================================================
# VICAST-Analyze: QC Workflow Only (Steps 1-6)
# =============================================================================
# Runs: mapping, variant calling, depth file, and diagnostic report
# STOPS before annotation to allow user review
#
# Prerequisites:
#   - SNPEFF_HOME or SNPEFF_JAR environment variable set
#   - Conda environment with required tools (bwa, samtools, lofreq, fastp)
#   - Target genome added to SnpEff database via VICAST-annotate
# =============================================================================

set -e  # Exit on error

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [OPTIONS]"
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
    echo "Options:"
    echo "  --large-files    High-memory mode for large files (1-5GB)"
    echo "  --extremely-large-files  Extreme memory mode for >5GB files"
    echo ""
    echo "Environment variables (set before running):"
    echo "  SNPEFF_HOME     - Path to SnpEff installation"
    echo "  SNPEFF_JAR      - Path to snpEff.jar (optional if SNPEFF_HOME set)"
    echo "  JAVA_HOME       - Path to Java installation (optional)"
    echo "  VICAST_CONDA_ENV - Conda environment name (default: vicast_analyze)"
    echo ""
    echo "After completion, review QC outputs and run annotation if satisfied:"
    echo "  ./run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=4
LARGE_FILES_FLAG=""

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
    else
        echo ""
        echo "ERROR: SnpEff not configured"
        echo "==============================="
        echo ""
        echo "Please set one of the following environment variables:"
        echo "  export SNPEFF_HOME=/path/to/snpEff"
        echo "  export SNPEFF_JAR=/path/to/snpEff/snpEff.jar"
        echo ""
        echo "Or copy and configure the VICAST configuration template:"
        echo "  mkdir -p ~/.vicast"
        echo "  cp ${VICAST_HOME}/vicast_config.template.sh ~/.vicast/config.sh"
        echo "  # Edit ~/.vicast/config.sh with your paths"
        echo "  source ~/.vicast/config.sh"
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

# =============================================================================
# Validation
# =============================================================================

# Validate input files
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

# Validate SnpEff JAR exists
if [ ! -f "$SNPEFF_JAR" ]; then
    echo ""
    echo "ERROR: SnpEff JAR not found: $SNPEFF_JAR"
    echo "========================================"
    echo ""
    echo "Please verify your SNPEFF_HOME or SNPEFF_JAR environment variable."
    echo ""
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
echo ""
echo "Configuration:"
echo "  SNPEFF_JAR: $SNPEFF_JAR"
echo "  SNPEFF_DATA: $SNPEFF_DATA"
echo "  JAVA_PATH: $JAVA_PATH"
echo "  Conda env: $VICAST_CONDA_ENV"
echo "========================================"
echo ""

# =============================================================================
# Pre-flight Checks
# =============================================================================

echo "Pre-flight check: snpEff database"
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

# Check for database files directly first (faster and works offline)
if [ -f "${SNPEFF_DATA}/${ACCESSION}/snpEffectPredictor.bin" ]; then
    echo "✓ Database found: ${SNPEFF_DATA}/${ACCESSION}/"
    echo "✓ Database is built and ready"
else
    echo ""
    echo "ERROR: Genome $ACCESSION database not found!"
    echo "========================================================"
    echo ""
    echo "Checked location: ${SNPEFF_DATA}/${ACCESSION}/"
    echo ""
    echo "Please install the database first:"
    echo "  Option 1 (pre-built database):"
    echo "    install_prebuilt_database.sh --install $ACCESSION"
    echo ""
    echo "  Option 2 (custom annotation):"
    echo "    cd ${VICAST_HOME}/vicast-annotate"
    echo "    python3 step1_parse_viral_genome.py $ACCESSION"
    echo "    # Edit the TSV file, then:"
    echo "    python3 step2_add_to_snpeff.py $ACCESSION ${ACCESSION}.tsv"
    echo ""
    exit 1
fi

# =============================================================================
# Run Pipeline
# =============================================================================

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

# =============================================================================
# Results Summary
# =============================================================================

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    # Extract sample name (match Python's logic: strip _R1/_R2/_1/_2 and _001 suffix)
    SAMPLE_NAME=$(basename "$R1" | sed -E 's/(_R?[12])?(_001)?\.fastq\.gz$//')

    # Use HOST_PWD for display paths (set via -e HOST_PWD=$(pwd) in docker run)
    DISPLAY_DIR="${HOST_PWD:-$(pwd)}"

    echo ""
    echo "========================================"
    echo "OK: QC WORKFLOW COMPLETED (Steps 1-6)"
    echo "========================================"
    echo ""
    echo "COMPLETED STEPS:"
    echo "  1. Prepare reference - Extract genome from SnpEff database"
    echo "  2. Read statistics - Count total reads, calculate coverage"
    echo "  3. Quality control - fastp trimming/filtering"
    echo "  4. Map & call variants - BWA mapping → lofreq variant calling"
    echo "  5. Coverage depth - Generate depth profile across genome"
    echo "  6. Diagnostic report - HTML report with QC metrics, contamination check"
    echo ""
    echo "REVIEW THE FOLLOWING QC OUTPUTS:"
    echo ""
    echo "1. Depth File:"
    echo "   ${DISPLAY_DIR}/${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt"
    echo ""
    echo "2. Diagnostic Report:"
    echo "   ${DISPLAY_DIR}/diagnostic_${SAMPLE_NAME}/${SAMPLE_NAME}_diagnostic_report.txt"
    echo "   ${DISPLAY_DIR}/diagnostic_${SAMPLE_NAME}/diagnostic_${SAMPLE_NAME}_presentation_ready_report.html"
    echo ""
    echo "3. VCF Files:"
    # Extract base sample name without _R1/_R2/_1/_2 suffix for VCF matching
    BASE_SAMPLE=$(echo "$SAMPLE_NAME" | sed -E 's/_[12]$//')
    find "$(pwd)/cleaned_seqs/variants" -name "${BASE_SAMPLE}*.vcf" 2>/dev/null | sed "s|$(pwd)|${DISPLAY_DIR}|" | head -10
    if [ -f "$(pwd)/cleaned_seqs/variants/${BASE_SAMPLE}_vars.vcf" ]; then
        echo "   Unfiltered: ${DISPLAY_DIR}/cleaned_seqs/variants/${BASE_SAMPLE}_vars.vcf"
    fi
    if [ -f "$(pwd)/cleaned_seqs/variants/${BASE_SAMPLE}_vars.filt.vcf" ]; then
        echo "   Filtered:   ${DISPLAY_DIR}/cleaned_seqs/variants/${BASE_SAMPLE}_vars.filt.vcf"
        # Convert filtered VCF to TSV with expanded INFO columns
        FILT_VCF="$(pwd)/cleaned_seqs/variants/${BASE_SAMPLE}_vars.filt.vcf"
        FILT_TSV="$(pwd)/cleaned_seqs/variants/${BASE_SAMPLE}_vars.filt.tsv"
        echo ""
        echo "4. Variant TSV (INFO fields expanded):"
        python "${PIPELINE_DIR}/vcf_to_tsv.py" "$FILT_VCF" -o "$FILT_TSV" --split-dp4 2>&1 || \
            echo "   WARNING: VCF-to-TSV conversion failed"
        if [ -f "$FILT_TSV" ]; then
            echo "   ${DISPLAY_DIR}/cleaned_seqs/variants/${BASE_SAMPLE}_vars.filt.tsv"
        fi
    fi
    echo ""
    echo "NEXT STEPS:"
    echo ""
    echo "  Review QC outputs, then continue with annotation (Steps 7-9):"
    echo ""
    echo "  Step 7 applies lofreq filter before annotation with these defaults:"
    echo "    --min-depth 200   Minimum read depth (lofreq -v)"
    echo "    --min-qual  1000  Minimum variant quality, phred (lofreq -Q/-K)"
    echo ""
    echo "  Basic annotation (recommended defaults):"
    echo "    run_vicast_analyze_annotate_only.sh $R1 $R2 $ACCESSION"
    echo ""
    echo "  With custom filter thresholds:"
    echo "    run_vicast_analyze_annotate_only.sh $R1 $R2 $ACCESSION --min-depth 100 --min-qual 50"
    echo ""
    echo "  Or submit as SLURM job:"
    echo "    sbatch --wrap=\"run_vicast_analyze_annotate_only.sh $R1 $R2 $ACCESSION\""
    echo ""
    echo "========================================"
else
    echo ""
    echo "ERROR: QC workflow failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi
