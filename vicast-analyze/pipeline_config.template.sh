#!/bin/bash
# Viral Genomics Pipeline Configuration Template
# Copy this file to pipeline_config.sh and update the paths for your system
# Usage: cp pipeline_config.template.sh pipeline_config.sh
#        Edit pipeline_config.sh with your actual paths
#        Source it in your scripts: source pipeline_config.sh

# =============================================================================
# CONDA/MAMBA ENVIRONMENT
# =============================================================================
# Conda/Mamba run command for running tools in the viral genomics environment
# This should use conda run or mamba run with your environment name
MAMBA_CMD="/path/to/conda_or_mamba run -n your_viral_env"

# Example for HTCF (using shared anaconda3):
# MAMBA_CMD="conda run -n viral_genomics"
# Note: Must have already run: source /ref/sahlab/software/anaconda3/bin/activate

# =============================================================================
# SNPEFF CONFIGURATION
# =============================================================================
# Directory where snpEff is installed
SNPEFF_DIR="/path/to/snpEff"

# Example for HTCF:
# SNPEFF_DIR="/ref/sahlab/software/snpEff"

# snpEff JAR file (usually in the snpEff directory)
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"

# Java executable (usually just 'java' if it's in your PATH)
JAVA_PATH="java"

# =============================================================================
# PIPELINE PATHS
# =============================================================================
# Base directory for the viral genomics pipeline
PIPELINE_BASE="/path/to/viral-genomics-pipeline"

# Example for HTCF:
# PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"

# Path to the consolidated pipeline script
CONSOLIDATED_PIPELINE="${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh"

# Path to the main viral_pipeline.py script
VIRAL_PIPELINE_SCRIPT="${PIPELINE_BASE}/viral_pipeline.py"

# =============================================================================
# OPTIONAL: WORKING DIRECTORIES
# =============================================================================
# Default working directory for pipeline outputs (can be overridden)
# WORK_DIR="/scratch/your_lab/your_username/viral_analysis"

# =============================================================================
# VALIDATION
# =============================================================================
# Function to validate configuration
validate_config() {
    local errors=0
    
    echo "Validating pipeline configuration..."
    
    # Check mamba
    if ! command -v ${MAMBA_CMD%% *} &> /dev/null; then
        echo "❌ Error: Mamba executable not found at ${MAMBA_CMD%% *}"
        errors=$((errors + 1))
    fi
    
    # Check snpEff
    if [ ! -f "$SNPEFF_JAR" ]; then
        echo "❌ Error: snpEff JAR not found at $SNPEFF_JAR"
        errors=$((errors + 1))
    fi
    
    # Check pipeline base
    if [ ! -d "$PIPELINE_BASE" ]; then
        echo "❌ Error: Pipeline directory not found at $PIPELINE_BASE"
        errors=$((errors + 1))
    fi
    
    # Check consolidated pipeline
    if [ ! -f "$CONSOLIDATED_PIPELINE" ]; then
        echo "❌ Error: Consolidated pipeline not found at $CONSOLIDATED_PIPELINE"
        errors=$((errors + 1))
    fi
    
    if [ $errors -eq 0 ]; then
        echo "✅ Configuration validated successfully!"
        return 0
    else
        echo "⚠️  Found $errors configuration error(s). Please update pipeline_config.sh"
        return 1
    fi
}

# Uncomment the line below to validate configuration when sourcing this file
# validate_config
