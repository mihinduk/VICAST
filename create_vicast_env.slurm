#!/bin/bash
#SBATCH --job-name=create_vicast_env
#SBATCH --output=create_vicast_env_%j.out
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=4

echo "=========================================="
echo "Creating VICAST conda environment"
echo "=========================================="
echo "Start time: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Memory: 32GB"
echo ""

# Activate conda
echo "Activating conda..."
source /ref/sahlab/software/anaconda3/bin/activate

# Extract pipeline directory from original sbatch command
# SLURM preserves the original command in job details
ORIGINAL_SCRIPT_PATH=$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K[^ ]+')
PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT_PATH")"

echo "Original script: $ORIGINAL_SCRIPT_PATH"
echo "Pipeline directory: $PIPELINE_DIR"
echo ""

cd "$PIPELINE_DIR"
echo "Working directory: $(pwd)"
echo ""

# Use minimal environment for faster creation (only VICAST-annotate essentials)
ENV_FILE="environment_minimal.yml"

# Check if environment file exists
if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: $ENV_FILE not found in $(pwd)"
    echo "Trying environment.yml instead..."
    ENV_FILE="environment.yml"
    if [ ! -f "$ENV_FILE" ]; then
        echo "ERROR: No environment file found"
        exit 1
    fi
fi

echo "Found $ENV_FILE"
echo ""

# Show what will be installed
echo "Environment specification:"
echo "------------------------"
head -20 "$ENV_FILE"
echo "..."
echo ""

# Check if mamba is available (much faster than conda)
if command -v mamba &> /dev/null; then
    echo "Found mamba - using faster solver"
    ENV_CREATOR="mamba"
else
    echo "Using conda (consider installing mamba for faster solving)"
    ENV_CREATOR="conda"
fi
echo ""

# Create the environment
echo "Creating environment with $ENV_CREATOR (this may take 10-30 minutes)..."
$ENV_CREATOR env create -f "$ENV_FILE"

if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "SUCCESS! Environment created"
    echo "=========================================="
    echo ""

    # Activate and verify
    echo "Activating vicast environment..."
    conda activate vicast

    echo "Verifying installation..."
    echo ""

    echo "Python version:"
    python --version
    echo ""

    echo "Testing package imports..."
    python -c "import Bio; print('✓ Biopython:', Bio.__version__)"
    python -c "import pandas; print('✓ Pandas:', pandas.__version__)"
    python -c "import numpy; print('✓ NumPy:', numpy.__version__)"
    echo ""

    echo "Checking BLAST installation..."
    which blastx && blastx -version | head -3
    echo ""

    echo "Environment location:"
    conda info --envs | grep vicast
    echo ""

    echo "To use this environment:"
    echo "  source /ref/sahlab/software/anaconda3/bin/activate"
    echo "  conda activate vicast"

else
    echo ""
    echo "=========================================="
    echo "ERROR: Environment creation failed"
    echo "=========================================="
    echo "Check the output above for errors"
    exit 1
fi

echo ""
echo "End time: $(date)"
