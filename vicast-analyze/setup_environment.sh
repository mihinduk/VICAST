#!/bin/bash
# Setup script for viral_genomics_analyze environment
# This script helps users create the conda environment in an appropriate location

set -e

echo "========================================="
echo "VICAST-ANALYZE Environment Setup"
echo "========================================="
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found"
    echo ""
    echo "Please install conda/mamba first or activate it:"
    echo "  source /path/to/anaconda3/bin/activate"
    echo ""
    echo "For HTCF users:"
    echo "  source $CONDA_BASE/bin/activate (or your conda installation)"
    exit 1
fi

echo "✅ Conda found: $(which conda)"
echo ""

# Warn about home directory
echo "⚠️  IMPORTANT: Conda environments are LARGE (>2GB)"
echo "   DO NOT install in your home directory on shared servers!"
echo ""

# Provide generic recommendations
echo "📍 Installation Location Recommendations:"
echo ""
echo "   Shared servers: Use lab/group space or scratch directory"
echo "   Local systems:  Home directory is acceptable"
echo ""

# Ask user for installation location
echo "Where should the environment be installed?"
echo ""
echo "Options:"
echo "  1) Shared lab/group space (enter path)"
echo "  2) Personal scratch/work space (enter path)"
echo "  3) Custom location (you specify)"
echo "  4) Default conda location (NOT recommended for shared servers)"
echo ""

read -p "Enter choice [1-4]: " choice

case $choice in
    1)
        read -p "Enter shared lab/group path (e.g., /shared/lab/envs): " ENV_DIR
        if [ -z "$ENV_DIR" ]; then
            echo "Error: Path cannot be empty"
            exit 1
        fi
        USE_NAMED=true
        ;;
    2)
        read -p "Enter scratch/work path (e.g., /scratch/$USER/envs): " ENV_DIR
        if [ -z "$ENV_DIR" ]; then
            echo "Error: Path cannot be empty"
            exit 1
        fi
        USE_NAMED=true
        ;;
    3)
        read -p "Enter custom path: " ENV_DIR
        if [ -z "$ENV_DIR" ]; then
            echo "Error: Path cannot be empty"
            exit 1
        fi
        read -p "Use named environment (y) or prefix path (n)? [y/n]: " use_named_input
        if [[ "$use_named_input" =~ ^[Yy]$ ]]; then
            USE_NAMED=true
        else
            USE_NAMED=false
        fi
        ;;
    4)
        echo ""
        echo "⚠️  WARNING: Installing to default location (~/.conda/envs/)"
        echo "   This may exceed home directory quota on shared servers!"
        read -p "Continue anyway? [y/N]: " confirm
        if [[ ! "$confirm" =~ ^[Yy]$ ]]; then
            echo "Aborted."
            exit 0
        fi
        ENV_DIR="$HOME/.conda/envs"
        USE_NAMED=true
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

# Create directory if it doesn't exist
if [ "$USE_NAMED" = true ]; then
    mkdir -p "$ENV_DIR"
    echo ""
    echo "Environment will be created at: $ENV_DIR/vicast_analyze"
    echo ""
    echo "Setting CONDA_ENVS_DIRS=$ENV_DIR:$HOME/.conda/envs"
    export CONDA_ENVS_DIRS="$ENV_DIR:$HOME/.conda/envs"

    # Create environment
    echo "Creating environment from vicast_analyze.yml..."
    conda env create -f vicast_analyze.yml

    ENV_ACTIVATE="conda activate vicast_analyze"
    MAMBA_CMD="conda run -n vicast_analyze"
else
    # Using --prefix
    mkdir -p "$ENV_DIR"
    FULL_ENV_PATH="$ENV_DIR/vicast_analyze"

    echo ""
    echo "Environment will be created at: $FULL_ENV_PATH"
    echo ""

    # Create environment with prefix
    echo "Creating environment from vicast_analyze.yml..."
    conda env create --prefix "$FULL_ENV_PATH" -f vicast_analyze.yml

    ENV_ACTIVATE="conda activate $FULL_ENV_PATH"
    MAMBA_CMD="conda run --prefix $FULL_ENV_PATH"
fi

echo ""
echo "========================================="
echo "✅ Environment created successfully!"
echo "========================================="
echo ""
echo "To activate this environment:"
echo "  $ENV_ACTIVATE"
echo ""
echo "To configure the pipeline, edit pipeline_config.sh:"
echo "  MAMBA_CMD=\"$MAMBA_CMD\""
echo ""
echo "Verify installation:"
echo "  $ENV_ACTIVATE"
echo "  bwa"
echo "  samtools --version"
echo "  lofreq version"
echo "  megahit --version"
echo "  blastn -version"
echo "  python -c \"from Bio import Entrez; print('BioPython OK')\""
echo ""
