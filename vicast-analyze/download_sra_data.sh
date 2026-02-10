#!/bin/bash

# =============================================================================
# VICAST SRA Data Download Utility
# Downloads sequencing data from NCBI SRA for testing VICAST pipelines
# =============================================================================

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Usage message
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] SRR_ACCESSION [OUTPUT_DIR]

Download sequencing data from NCBI SRA and convert to FASTQ format.

Arguments:
    SRR_ACCESSION   SRA run accession (e.g., SRR5992153)
    OUTPUT_DIR      Output directory (default: current directory)

Options:
    -h, --help      Show this help message
    -k, --keep-sra  Keep the .sra file after conversion (default: delete)
    -t, --threads   Number of threads for fasterq-dump (default: 4)

Examples:
    # Download DENV-2 strain 16681 dataset
    $(basename "$0") SRR5992153

    # Download to specific directory with 8 threads
    $(basename "$0") -t 8 SRR5992153 /data/test_data

    # Download and keep the .sra file
    $(basename "$0") --keep-sra SRR5992153

Recommended DENV Test Datasets:
    SRR5992153 - DENV-2 strain 16681, Vero cells, paired-end, high coverage

EOF
    exit 0
}

# Default values
KEEP_SRA=false
THREADS=4

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            ;;
        -k|--keep-sra)
            KEEP_SRA=true
            shift
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -*)
            echo -e "${RED}Error: Unknown option: $1${NC}" >&2
            usage
            ;;
        *)
            break
            ;;
    esac
done

# Check required argument
if [[ $# -lt 1 ]]; then
    echo -e "${RED}Error: SRR accession required${NC}" >&2
    usage
fi

SRR_ACC="$1"
OUTPUT_DIR="${2:-.}"

# Validate SRR accession format
if [[ ! $SRR_ACC =~ ^[SED]RR[0-9]+$ ]]; then
    echo -e "${RED}Error: Invalid SRR accession format: $SRR_ACC${NC}" >&2
    echo "Expected format: SRR###### or ERR###### or DRR######" >&2
    exit 1
fi

# Check if sra-tools is installed
if ! command -v prefetch &> /dev/null || ! command -v fasterq-dump &> /dev/null; then
    echo -e "${RED}Error: sra-tools not installed${NC}" >&2
    echo "Install with: conda install -c bioconda sra-tools" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo -e "${GREEN}=== VICAST SRA Download Utility ===${NC}"
echo "Accession: $SRR_ACC"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

# Step 1: Prefetch the SRA file
echo -e "${YELLOW}[1/3] Downloading SRA file...${NC}"
if prefetch "$SRR_ACC" --max-size 100G; then
    echo -e "${GREEN}✓ Download complete${NC}"
else
    echo -e "${RED}✗ Download failed${NC}" >&2
    exit 1
fi

# Step 2: Convert to FASTQ
echo -e "${YELLOW}[2/3] Converting to FASTQ...${NC}"
if fasterq-dump "$SRR_ACC" --threads "$THREADS" --split-files --progress; then
    echo -e "${GREEN}✓ Conversion complete${NC}"
else
    echo -e "${RED}✗ Conversion failed${NC}" >&2
    exit 1
fi

# Step 3: Compress FASTQ files
echo -e "${YELLOW}[3/3] Compressing FASTQ files...${NC}"
if ls "${SRR_ACC}"*.fastq 1> /dev/null 2>&1; then
    gzip -v "${SRR_ACC}"*.fastq
    echo -e "${GREEN}✓ Compression complete${NC}"
else
    echo -e "${RED}✗ No FASTQ files found${NC}" >&2
    exit 1
fi

# Clean up .sra file unless --keep-sra specified
if [[ "$KEEP_SRA" == false ]]; then
    echo -e "${YELLOW}Removing .sra file...${NC}"
    rm -rf ~/ncbi/public/sra/"${SRR_ACC}.sra"
    echo -e "${GREEN}✓ Cleanup complete${NC}"
fi

# Display results
echo ""
echo -e "${GREEN}=== Download Complete ===${NC}"
echo "Output files:"
ls -lh "${SRR_ACC}"*.fastq.gz

# Count reads
echo ""
echo "Read counts:"
for file in "${SRR_ACC}"*.fastq.gz; do
    if [[ -f "$file" ]]; then
        count=$(zcat "$file" | wc -l | awk '{print $1/4}')
        printf "  %-30s %'d reads\n" "$file" "$count"
    fi
done

echo ""
echo -e "${GREEN}Ready to run VICAST analysis!${NC}"
echo ""
echo "Example commands:"
echo "  # Quality check only"
echo "  run_vicast_analyze_qc_only.sh ${SRR_ACC}_1.fastq.gz ${SRR_ACC}_2.fastq.gz NC_001474.2 $THREADS"
echo ""
echo "  # Full analysis pipeline"
echo "  run_vicast_analyze_full.sh ${SRR_ACC}_1.fastq.gz ${SRR_ACC}_2.fastq.gz NC_001474.2 $THREADS"
