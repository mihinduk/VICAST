#!/bin/bash
# Install Java 21 using conda on HTCF

echo "========================================="
echo "Installing Java 21 via Conda on HTCF"
echo "========================================="

# Activate conda on HTCF
source /ref/sahlab/software/anaconda3/bin/activate

# Create conda environment with Java 21
ENV_NAME="java21_snpeff"

echo "Creating conda environment '$ENV_NAME' with Java 21..."
conda create -n $ENV_NAME -c conda-forge openjdk=21 -y

echo ""
echo "========================================="
echo "Installation Complete!"
echo "========================================="
echo ""
echo "Java 21 has been installed in conda environment: $ENV_NAME"
echo ""
echo "To use this Java installation:"
echo "  1. Activate conda: source /ref/sahlab/software/anaconda3/bin/activate"
echo "  2. Activate environment: conda activate $ENV_NAME"
echo "  3. Java path will be: \$CONDA_PREFIX"
echo ""
echo "Your paths for setup_snpeff_custom_paths.sh:"
echo "  Java path: \$CONDA_PREFIX (after activating the environment)"
echo "  snpEff path: /home/mihindu/software/snpEff"
echo ""
echo "Full workflow:"
echo "  source /ref/sahlab/software/anaconda3/bin/activate"
echo "  conda activate $ENV_NAME"
echo "  cd /scratch/sahlab/kathie/viral_genome_snpeff_update"
echo "  ./setup_snpeff_custom_paths.sh"
echo "  # Enter Java path: \$CONDA_PREFIX"
echo "  # Enter snpEff path: /home/mihindu/software/snpEff"
