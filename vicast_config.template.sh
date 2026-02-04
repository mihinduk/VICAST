#!/bin/bash
# =============================================================================
# VICAST Configuration Template
# =============================================================================
#
# Copy this file to ~/.vicast/config.sh or source it in your ~/.bashrc
#
# Usage:
#   cp vicast_config.template.sh ~/.vicast/config.sh
#   # Edit ~/.vicast/config.sh with your paths
#   source ~/.vicast/config.sh
#
# Or add to ~/.bashrc:
#   source /path/to/VICAST/vicast_config.template.sh
#
# =============================================================================

# -----------------------------------------------------------------------------
# SNPEFF CONFIGURATION (Required for annotation)
# -----------------------------------------------------------------------------
# Set SNPEFF_HOME to your SnpEff installation directory
# The JAR and data paths will be derived automatically if not set

# export SNPEFF_HOME="/path/to/snpEff"
# export SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
# export SNPEFF_DATA="${SNPEFF_HOME}/data"

# Example paths for common setups:
#   Linux (manual install):  SNPEFF_HOME="/opt/snpEff"
#   Linux (conda):           SNPEFF_HOME="${CONDA_PREFIX}/share/snpeff-5.2-0"
#   macOS (homebrew):        SNPEFF_HOME="/usr/local/opt/snpeff/libexec"

# -----------------------------------------------------------------------------
# JAVA CONFIGURATION (Required - Java 21+ for SnpEff 5.2+)
# -----------------------------------------------------------------------------
# Only needed if java is not in PATH or you need a specific version

# export JAVA_HOME="/path/to/jdk-21"
# export PATH="${JAVA_HOME}/bin:${PATH}"

# Example paths:
#   Linux (manual):    JAVA_HOME="/usr/lib/jvm/java-21-openjdk"
#   Linux (conda):     JAVA_HOME="${CONDA_PREFIX}"
#   macOS (homebrew):  JAVA_HOME="$(/usr/libexec/java_home -v 21)"

# -----------------------------------------------------------------------------
# NCBI ENTREZ CONFIGURATION (Required for genome downloads)
# -----------------------------------------------------------------------------
# NCBI requires an email address for Entrez queries
# Set this to your institutional email

export NCBI_EMAIL="${NCBI_EMAIL:-your.email@institution.edu}"

# Optional: NCBI API key for higher rate limits
# export NCBI_API_KEY="your_api_key_here"

# -----------------------------------------------------------------------------
# WORKING DIRECTORIES (Recommended for HPC environments)
# -----------------------------------------------------------------------------
# Scratch directory for temporary files and large outputs
# On HPC systems, this should be a fast local or parallel filesystem

# export VICAST_SCRATCH="/scratch/${USER}/vicast"
# export VICAST_SCRATCH="/tmp/vicast_${USER}"

# Ensure scratch directory exists
if [ -n "${VICAST_SCRATCH}" ]; then
    mkdir -p "${VICAST_SCRATCH}" 2>/dev/null
fi

# -----------------------------------------------------------------------------
# RESOURCE LIMITS
# -----------------------------------------------------------------------------
# Default number of threads for parallel operations
export VICAST_THREADS="${VICAST_THREADS:-4}"

# Memory settings for Java (SnpEff)
# export _JAVA_OPTIONS="-Xmx8g"

# -----------------------------------------------------------------------------
# CONDA ENVIRONMENT (Optional)
# -----------------------------------------------------------------------------
# If using conda, specify the environment name or path

# export VICAST_CONDA_ENV="vicast"
# export VICAST_CONDA_ENV="/path/to/envs/vicast"

# Auto-activate conda environment (optional)
# if command -v conda &> /dev/null && [ -n "${VICAST_CONDA_ENV}" ]; then
#     conda activate "${VICAST_CONDA_ENV}"
# fi

# -----------------------------------------------------------------------------
# BLAST DATABASE (Optional - for Pathway 3 and contamination detection)
# -----------------------------------------------------------------------------
# Path to BLAST databases for BLASTx annotation

# export BLASTDB="/path/to/blast/databases"

# -----------------------------------------------------------------------------
# VICAST INSTALLATION PATH (Auto-detected)
# -----------------------------------------------------------------------------
# This is automatically set when sourcing this file

if [ -n "${BASH_SOURCE[0]}" ]; then
    export VICAST_HOME="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
elif [ -n "${0}" ]; then
    export VICAST_HOME="$(cd "$(dirname "${0}")" && pwd)"
fi

# Add VICAST to PATH (optional)
# export PATH="${VICAST_HOME}/vicast-annotate:${VICAST_HOME}/vicast-analyze:${PATH}"

# -----------------------------------------------------------------------------
# VALIDATION
# -----------------------------------------------------------------------------
# Uncomment to validate configuration on source

# vicast_validate_config() {
#     local errors=0
#
#     if [ -z "${SNPEFF_JAR}" ] || [ ! -f "${SNPEFF_JAR}" ]; then
#         echo "Warning: SNPEFF_JAR not configured or not found" >&2
#         errors=$((errors + 1))
#     fi
#
#     if ! command -v java &> /dev/null; then
#         echo "Warning: java not found in PATH" >&2
#         errors=$((errors + 1))
#     fi
#
#     if [ "${NCBI_EMAIL}" = "your.email@institution.edu" ]; then
#         echo "Warning: NCBI_EMAIL not configured" >&2
#         errors=$((errors + 1))
#     fi
#
#     return ${errors}
# }
#
# vicast_validate_config

# -----------------------------------------------------------------------------
# SITE-SPECIFIC CONFIGURATION
# -----------------------------------------------------------------------------
# Add your institution-specific settings below this line
#


