#!/bin/bash
# Viral Genomics Pipeline Configuration
# Created from pipeline_config.template.sh
# This file contains HTCF-specific paths

# =============================================================================
# CONDA/MAMBA ENVIRONMENT
# =============================================================================
# Conda run command for running tools in the viral_genomics environment
# Must source anaconda3 first, then use conda run
MAMBA_CMD="conda run -n viral_genomics"

# =============================================================================
# SNPEFF CONFIGURATION
# =============================================================================
# Directory where snpEff is installed
SNPEFF_DIR="/ref/sahlab/software/snpEff"

# snpEff JAR file (in the snpEff directory)
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"

# Java executable
JAVA_PATH="java"

# =============================================================================
# PIPELINE PATHS
# =============================================================================
# Base directory for the viral genomics pipeline
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"

# Path to the consolidated pipeline script
CONSOLIDATED_PIPELINE="${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh"

# Path to the main viral_pipeline.py script
VIRAL_PIPELINE_SCRIPT="${PIPELINE_BASE}/viral_pipeline.py"

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
