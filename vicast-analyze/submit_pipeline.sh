#!/bin/bash
# VICAST Pipeline Submission Script
# Provides option to run interactively or as SLURM job

usage() {
    cat << 'USAGE'
Usage: ./submit_pipeline.sh [OPTIONS] <R1_fastq> <R2_fastq> <accession> <threads>

Run vicast-analyze pipeline either interactively or as a SLURM job.

Arguments:
  R1_fastq      : Path to R1 FASTQ file
  R2_fastq      : Path to R2 FASTQ file
  accession     : GenBank accession number (e.g., NC_001477.1)
  threads       : Number of threads to use

Options:
  --slurm       : Submit as SLURM job (RECOMMENDED for large files)
  --interactive : Run interactively (WARNING: may run out of memory on login node)
  --mem SIZE    : Memory to request for SLURM job (default: 64G)
  --time TIME   : Time limit for SLURM job (default: 4:00:00)
  --large-files : Enable large file handling mode

Examples:
  # Submit as SLURM job (recommended):
  ./submit_pipeline.sh --slurm sample_R1.fastq.gz sample_R2.fastq.gz NC_009942.1 8

  # Run interactively (small files only):
  ./submit_pipeline.sh --interactive sample_R1.fastq.gz sample_R2.fastq.gz NC_009942.1 4

  # Custom resources:
  ./submit_pipeline.sh --slurm --mem 128G --time 8:00:00 sample_R1.fastq.gz sample_R2.fastq.gz NC_009942.1 16

USAGE
    exit 1
}

# Default values
MODE=""
MEM="64G"
TIME="4:00:00"
LARGE_FILES_FLAG=""

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --slurm)
            MODE="slurm"
            shift
            ;;
        --interactive)
            MODE="interactive"
            shift
            ;;
        --mem)
            MEM="$2"
            shift 2
            ;;
        --time)
            TIME="$2"
            shift 2
            ;;
        --large-files)
            LARGE_FILES_FLAG="--large-files"
            shift
            ;;
        --help|-h)
            usage
            ;;
        *)
            break
            ;;
    esac
done

# Check remaining arguments
if [ $# -lt 4 ]; then
    echo "ERROR: Missing required arguments"
    echo ""
    usage
fi

R1=$1
R2=$2
ACCESSION=$3
THREADS=$4

# Validate input files
if [ ! -f "$R1" ]; then
    echo "ERROR: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "ERROR: R2 file not found: $R2"
    exit 1
fi

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PIPELINE_SCRIPT="${SCRIPT_DIR}/run_pipeline_htcf_enhanced.sh"

if [ ! -f "$PIPELINE_SCRIPT" ]; then
    echo "ERROR: Pipeline script not found: $PIPELINE_SCRIPT"
    exit 1
fi

# Extract sample name
SAMPLE_NAME=$(basename "$R1" | sed -E "s/_R?[12](_[0-9]+)?..*$//")

# Ask user for mode if not specified
if [ -z "$MODE" ]; then
    echo "========================================="
    echo "VICAST Pipeline Submission"
    echo "========================================="
    echo ""
    echo "How would you like to run the pipeline?"
    echo ""
    echo "  1 - SLURM job (RECOMMENDED)"
    echo "     Runs on compute node with dedicated resources"
    echo "     Won't get killed for memory issues"
    echo "     Can monitor with: tail -f vicast_${SAMPLE_NAME}_*.out"
    echo ""
    echo "  2 - Interactive"
    echo "     Runs immediately on current node"
    echo "     WARNING: May fail with large files - OOM"
    echo "     Only use for small test files"
    echo ""
    read -p "Enter choice [1-2]: " choice

    case $choice in
        1)
            MODE="slurm"
            ;;
        2)
            MODE="interactive"
            echo ""
            echo "WARNING: Running interactively may cause Out of Memory errors!"
            read -p "Continue? [y/N]: " confirm
            if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
                echo "Aborted."
                exit 0
            fi
            ;;
        *)
            echo "Invalid choice. Exiting."
            exit 1
            ;;
    esac
fi

# Execute based on mode
if [ "$MODE" == "slurm" ]; then
    echo "========================================="
    echo "Submitting pipeline as SLURM job"
    echo "========================================="
    echo "Sample: $SAMPLE_NAME"
    echo "Memory: $MEM"
    echo "Time: $TIME"
    echo "Threads: $THREADS"
    echo "Working directory: $(pwd)"
    echo ""

    # Create job script
    JOB_SCRIPT="vicast_${SAMPLE_NAME}_job.sh"

    cat > "$JOB_SCRIPT" << 'EOFSLURM'
#!/bin/bash
#SBATCH --job-name=vicast_SAMPLE
#SBATCH --output=vicast_SAMPLE_%j.out
#SBATCH --error=vicast_SAMPLE_%j.err
#SBATCH --time=TIME_PLACEHOLDER
#SBATCH --mem=MEM_PLACEHOLDER
#SBATCH --cpus-per-task=THREADS_PLACEHOLDER

cd WORKDIR_PLACEHOLDER

echo "========================================="
echo "VICAST Pipeline - SLURM Job"
echo "========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo "========================================="
echo ""

# Activate conda/micromamba environment
echo "Activating environment..."
if command -v micromamba &> /dev/null; then
    echo "Detected micromamba (Docker environment) - using current environment"
elif command -v conda &> /dev/null; then
    eval "\$(conda shell.bash hook)"
    conda activate vicast_analyze 2>/dev/null || echo "Using current environment"
elif [ -n "\$CONDA_BASE" ] && [ -f "\$CONDA_BASE/bin/activate" ]; then
    source "\$CONDA_BASE/bin/activate"
    conda activate vicast_analyze 2>/dev/null || echo "Using current environment"
fi
echo ""

bash PIPELINE_SCRIPT_PLACEHOLDER \
    R1_PLACEHOLDER \
    R2_PLACEHOLDER \
    ACCESSION_PLACEHOLDER \
    THREADS_PLACEHOLDER \
    LARGE_FILES_PLACEHOLDER

EXIT_CODE=$?

echo ""
echo "========================================="
echo "Job completed: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================="

exit $EXIT_CODE
EOFSLURM

    # Replace placeholders
    sed -i "s/SAMPLE/${SAMPLE_NAME}/g" "$JOB_SCRIPT"
    sed -i "s|TIME_PLACEHOLDER|${TIME}|" "$JOB_SCRIPT"
    sed -i "s|MEM_PLACEHOLDER|${MEM}|" "$JOB_SCRIPT"
    sed -i "s|THREADS_PLACEHOLDER|${THREADS}|" "$JOB_SCRIPT"
    sed -i "s|WORKDIR_PLACEHOLDER|$(pwd)|" "$JOB_SCRIPT"
    sed -i "s|PIPELINE_SCRIPT_PLACEHOLDER|${PIPELINE_SCRIPT}|" "$JOB_SCRIPT"
    sed -i "s|R1_PLACEHOLDER|${R1}|" "$JOB_SCRIPT"
    sed -i "s|R2_PLACEHOLDER|${R2}|" "$JOB_SCRIPT"
    sed -i "s|ACCESSION_PLACEHOLDER|${ACCESSION}|" "$JOB_SCRIPT"
    sed -i "s|LARGE_FILES_PLACEHOLDER|${LARGE_FILES_FLAG}|" "$JOB_SCRIPT"

    # Submit job
    JOB_ID=$(sbatch "$JOB_SCRIPT" | awk '{print $NF}')

    echo "Job submitted: $JOB_ID"
    echo ""
    echo "Monitor progress:"
    echo "  tail -f vicast_${SAMPLE_NAME}_${JOB_ID}.out"
    echo ""
    echo "Check job status:"
    echo "  squeue -j $JOB_ID"
    echo ""
    echo "Cancel job:"
    echo "  scancel $JOB_ID"
    echo ""

elif [ "$MODE" == "interactive" ]; then
    echo "========================================="
    echo "Running pipeline interactively"
    echo "========================================="
    echo "Sample: $SAMPLE_NAME"
    echo "Threads: $THREADS"
    echo "WARNING: Running on: $(hostname)"
    echo ""

    bash "$PIPELINE_SCRIPT" "$R1" "$R2" "$ACCESSION" "$THREADS" $LARGE_FILES_FLAG

    EXIT_CODE=$?
    echo ""
    echo "========================================="
    echo "Pipeline completed with exit code: $EXIT_CODE"
    echo "========================================="
    exit $EXIT_CODE
else
    echo "ERROR: Invalid mode: $MODE"
    exit 1
fi
