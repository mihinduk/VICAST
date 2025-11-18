#!/bin/bash
# Setup script for viral_genomics_analyze environment
# This script helps users create the conda environment in an appropriate location

set -e

echo "========================================="
echo "VICAST-Analyze Environment Setup"
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
    echo "  source /ref/sahlab/software/anaconda3/bin/activate"
    exit 1
fi

echo "✅ Conda found: $(which conda)"
echo ""

# Warn about home directory
echo "⚠️  IMPORTANT: Conda environments are LARGE (>2GB)"
echo "   DO NOT install in your home directory on shared servers!"
echo ""

# Detect if on HTCF
if [ -d "/ref/sahlab/software" ]; then
    echo "📍 HTCF detected - Recommending shared lab installation"
    echo ""
    DEFAULT_ENV_DIR="/ref/sahlab/software/envs"
    IS_HTCF=true
else
    echo "📍 Non-HTCF system detected"
    echo ""
    DEFAULT_ENV_DIR="/scratch/your_lab/your_username/envs"
    IS_HTCF=false
fi

# Ask user for installation location
echo "Where should the environment be installed?"
echo ""
echo "Options:"
if [ "$IS_HTCF" = true ]; then
    echo "  1) Shared lab space (RECOMMENDED): /ref/sahlab/software/envs"
    echo "  2) Your scratch space: /scratch/sahlab/\$USER/envs"
else
    echo "  1) Shared lab space: /path/to/shared/envs"
    echo "  2) Your scratch space: /scratch/your_lab/\$USER/envs"
fi
echo "  3) Custom location (you specify)"
echo "  4) Default conda location (NOT recommended for shared servers)"
echo ""

read -p "Enter choice [1-4] (default: 1): " choice
choice=${choice:-1}

case $choice in
    1)
        if [ "$IS_HTCF" = true ]; then
            ENV_DIR="/ref/sahlab/software/envs"
        else
            read -p "Enter shared lab path: " ENV_DIR
        fi
        USE_NAMED=true
        ;;
    2)
        if [ "$IS_HTCF" = true ]; then
            ENV_DIR="/scratch/sahlab/$USER/envs"
        else
            read -p "Enter scratch path: " ENV_DIR
        fi
        USE_NAMED=true
        ;;
    3)
        read -p "Enter custom path: " ENV_DIR
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
    echo "Environment will be created at: $ENV_DIR/viral_genomics_analyze"
    echo ""
    echo "Setting CONDA_ENVS_DIRS=$ENV_DIR:$HOME/.conda/envs"
    export CONDA_ENVS_DIRS="$ENV_DIR:$HOME/.conda/envs"

    # Create environment
    echo "Creating environment from viral_genomics_analyze.yml..."
    conda env create -f viral_genomics_analyze.yml

    ENV_ACTIVATE="conda activate viral_genomics_analyze"
    MAMBA_CMD="conda run -n viral_genomics_analyze"
else
    # Using --prefix
    mkdir -p "$ENV_DIR"
    FULL_ENV_PATH="$ENV_DIR/viral_genomics_analyze"

    echo ""
    echo "Environment will be created at: $FULL_ENV_PATH"
    echo ""

    # Create environment with prefix
    echo "Creating environment from viral_genomics_analyze.yml..."
    conda env create --prefix "$FULL_ENV_PATH" -f viral_genomics_analyze.yml

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
echo "  which efetch"
echo ""
