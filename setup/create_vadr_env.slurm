#!/bin/bash
#SBATCH --job-name=create_vicast_vadr_env
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=create_vadr_env_%j.out
#SBATCH --error=create_vadr_env_%j.err

# Create vicast-vadr conda environment with sufficient memory

echo "Starting vicast-vadr environment creation at $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# Activate conda
source /ref/sahlab/software/anaconda3/bin/activate

# Create environment
echo "Creating vicast-vadr environment..."
conda env create -f environment_vadr.yml

# Check if successful
if [ $? -eq 0 ]; then
    echo "✓ Successfully created vicast-vadr environment"
    echo "To use: conda activate vicast-vadr"
else
    echo "✗ Failed to create environment"
    exit 1
fi

echo "Completed at $(date)"
