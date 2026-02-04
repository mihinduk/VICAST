#!/usr/bin/env bash
#
# VICAST SnpEff Setup Script
#
# This script downloads and configures SnpEff for use with VICAST.
# It can also build custom genome databases from user-provided files.
#
# Usage:
#   ./tools/setup_snpeff.sh [--install-only]
#   ./tools/setup_snpeff.sh --build-genome GENOME_NAME FASTA GFF3
#
# Examples:
#   # Install SnpEff only
#   ./tools/setup_snpeff.sh --install-only
#
#   # Install and build a custom genome database
#   ./tools/setup_snpeff.sh --build-genome my_virus genome.fasta genes.gff3
#

set -e

# Defaults
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VICAST_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SNPEFF_DIR="${SNPEFF_HOME:-$SCRIPT_DIR/snpEff}"
SNPEFF_VERSION="5.2"
SNPEFF_URL="https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip/download"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }
log_step() { echo -e "${BLUE}[STEP]${NC} $1"; }

usage() {
    echo "VICAST SnpEff Setup Script"
    echo ""
    echo "Usage:"
    echo "  $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  --install-only              Only install SnpEff, don't build any genomes"
    echo "  --build-genome NAME FA GFF  Build a custom genome database"
    echo "  --snpeff-dir DIR            Install SnpEff to DIR (default: $SNPEFF_DIR)"
    echo "  --help, -h                  Show this help message"
    echo ""
    echo "Environment Variables:"
    echo "  SNPEFF_HOME                 Override default SnpEff installation directory"
    echo ""
    echo "Examples:"
    echo "  # Install SnpEff"
    echo "  $0 --install-only"
    echo ""
    echo "  # Build a custom genome database"
    echo "  $0 --build-genome dengue_1 genome.fasta genes.gff3"
    echo ""
}

check_java() {
    if ! command -v java &> /dev/null; then
        log_error "Java not found. Please install Java 8 or higher."
        echo "  On macOS: brew install openjdk"
        echo "  On Ubuntu: sudo apt install default-jdk"
        exit 1
    fi

    JAVA_VERSION=$(java -version 2>&1 | head -1)
    log_info "Java: $JAVA_VERSION"
}

install_snpeff() {
    log_step "Installing SnpEff..."

    if [ -f "$SNPEFF_DIR/snpEff.jar" ]; then
        log_info "SnpEff already installed at $SNPEFF_DIR"
        return 0
    fi

    mkdir -p "$(dirname "$SNPEFF_DIR")"
    cd "$(dirname "$SNPEFF_DIR")"

    log_info "Downloading SnpEff..."
    curl -L -o snpEff_latest_core.zip "$SNPEFF_URL"

    log_info "Extracting..."
    unzip -q snpEff_latest_core.zip
    rm snpEff_latest_core.zip

    # Handle case where extraction creates nested directory
    if [ -d "snpEff/snpEff" ]; then
        mv snpEff/snpEff/* snpEff/
        rmdir snpEff/snpEff
    fi

    log_info "SnpEff installed to $SNPEFF_DIR"

    # Test installation
    if java -jar "$SNPEFF_DIR/snpEff.jar" -version 2>&1 | head -1; then
        log_info "SnpEff installation verified"
    else
        log_error "SnpEff installation failed"
        exit 1
    fi
}

build_genome() {
    local GENOME_NAME="$1"
    local FASTA_PATH="$2"
    local GFF_PATH="$3"

    log_step "Building genome database: $GENOME_NAME"

    # Validate inputs
    if [ ! -f "$FASTA_PATH" ]; then
        log_error "FASTA file not found: $FASTA_PATH"
        exit 1
    fi

    if [ ! -f "$GFF_PATH" ]; then
        log_error "GFF3 file not found: $GFF_PATH"
        exit 1
    fi

    # Validate VCF REF bases if vicast validation is available
    if command -v python3 &> /dev/null; then
        PYTHONPATH="$VICAST_ROOT/src" python3 -c "
from vicast.validation import validate_gff_for_snpeff
is_valid, errors, warnings = validate_gff_for_snpeff('$GFF_PATH', '$FASTA_PATH')
if not is_valid:
    print('GFF3 validation FAILED:')
    for e in errors:
        print(f'  ERROR: {e}')
    exit(1)
for w in warnings:
    print(f'  WARN: {w}')
print('GFF3 validation passed')
" 2>/dev/null || log_warn "Skipping VICAST validation (not installed)"
    fi

    # Add to snpEff.config if not present
    if ! grep -q "^${GENOME_NAME}\.genome" "$SNPEFF_DIR/snpEff.config" 2>/dev/null; then
        echo "" >> "$SNPEFF_DIR/snpEff.config"
        echo "# VICAST custom genome" >> "$SNPEFF_DIR/snpEff.config"
        echo "${GENOME_NAME}.genome : ${GENOME_NAME}" >> "$SNPEFF_DIR/snpEff.config"
        log_info "Added $GENOME_NAME to snpEff.config"
    else
        log_info "$GENOME_NAME already in snpEff.config"
    fi

    # Create data directory
    DATA_DIR="$SNPEFF_DIR/data/$GENOME_NAME"
    mkdir -p "$DATA_DIR"

    # Copy files
    cp "$FASTA_PATH" "$DATA_DIR/sequences.fa"
    cp "$GFF_PATH" "$DATA_DIR/genes.gff"

    log_info "Copied reference files to $DATA_DIR"

    # Build database
    log_info "Building SnpEff database (this may take a moment)..."
    java -jar "$SNPEFF_DIR/snpEff.jar" build -gff3 -v "$GENOME_NAME" 2>&1 | \
        grep -E "(Done|Error|Warning|Protein|CDS)" || true

    log_info "Database built successfully"

    # Verify
    log_info "Verifying database..."
    java -jar "$SNPEFF_DIR/snpEff.jar" databases 2>/dev/null | grep -q "$GENOME_NAME" && \
        log_info "Database $GENOME_NAME is available" || \
        log_warn "Database verification inconclusive"
}

# Parse arguments
INSTALL_ONLY=false
BUILD_GENOME=""
GENOME_FASTA=""
GENOME_GFF=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --install-only)
            INSTALL_ONLY=true
            shift
            ;;
        --build-genome)
            BUILD_GENOME="$2"
            GENOME_FASTA="$3"
            GENOME_GFF="$4"
            shift 4
            ;;
        --snpeff-dir)
            SNPEFF_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Main
echo "=============================================="
echo "VICAST SnpEff Setup"
echo "=============================================="
echo ""

check_java
install_snpeff

if [ "$INSTALL_ONLY" = true ]; then
    echo ""
    log_info "Installation complete"
    echo ""
    echo "To use SnpEff, set the following environment variable:"
    echo "  export SNPEFF_HOME=\"$SNPEFF_DIR\""
    echo ""
    echo "Or add to your shell profile (~/.bashrc or ~/.zshrc):"
    echo "  echo 'export SNPEFF_HOME=\"$SNPEFF_DIR\"' >> ~/.bashrc"
    exit 0
fi

if [ -n "$BUILD_GENOME" ]; then
    build_genome "$BUILD_GENOME" "$GENOME_FASTA" "$GENOME_GFF"
fi

echo ""
log_info "Setup complete"
echo ""
echo "To annotate variants:"
echo "  java -jar $SNPEFF_DIR/snpEff.jar GENOME_NAME input.vcf > output.vcf"
echo ""
