#!/bin/bash
################################################################################
# Example Usage: check_read_cooccurrence.py
#
# This script demonstrates common usage patterns for the read co-occurrence
# analysis tool.
################################################################################

# Exit on error
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo "=================================="
echo "Read Co-occurrence Analysis Examples"
echo "=================================="
echo ""

# Check if files exist
if [ ! -f "$1" ]; then
    echo -e "${RED}ERROR: Please provide a BAM file as first argument${NC}"
    echo "Usage: $0 <bam_file> <vcf_file>"
    echo ""
    echo "Example:"
    echo "  $0 sample_aligned_sorted.bam sample_annotated_filtered.tsv"
    exit 1
fi

if [ ! -f "$2" ]; then
    echo -e "${RED}ERROR: Please provide a VCF/TSV file as second argument${NC}"
    echo "Usage: $0 <bam_file> <vcf_file>"
    exit 1
fi

BAM_FILE="$1"
VCF_FILE="$2"
BASENAME=$(basename "$BAM_FILE" .bam)

echo -e "${GREEN}Input files:${NC}"
echo "  BAM: $BAM_FILE"
echo "  VCF: $VCF_FILE"
echo ""

# Check if BAM is indexed
if [ ! -f "${BAM_FILE}.bai" ]; then
    echo -e "${YELLOW}BAM index not found. Creating index...${NC}"
    samtools index "$BAM_FILE"
    echo -e "${GREEN}✓ Index created${NC}"
fi

# Example 1: Basic usage
echo "=================================="
echo "Example 1: Basic Analysis"
echo "=================================="
echo "Analyze all variant pairs with default settings"
echo ""
echo "Command:"
echo "  python check_read_cooccurrence.py \\"
echo "    --bam $BAM_FILE \\"
echo "    --vcf $VCF_FILE \\"
echo "    --output ${BASENAME}_cooccurrence_basic.tsv"
echo ""
read -p "Press Enter to run (or Ctrl+C to skip)..."
python check_read_cooccurrence.py \
    --bam "$BAM_FILE" \
    --vcf "$VCF_FILE" \
    --output "${BASENAME}_cooccurrence_basic.tsv"
echo -e "${GREEN}✓ Complete${NC}"
echo ""

# Example 2: High-quality variants only
echo "=================================="
echo "Example 2: High-Quality Variants Only"
echo "=================================="
echo "Analyze only high-confidence variants"
echo ""
echo "Command:"
echo "  python check_read_cooccurrence.py \\"
echo "    --bam $BAM_FILE \\"
echo "    --vcf $VCF_FILE \\"
echo "    --min-qual 1000 \\"
echo "    --min-depth 200 \\"
echo "    --min-freq 0.03 \\"
echo "    --output ${BASENAME}_cooccurrence_highqual.tsv"
echo ""
read -p "Press Enter to run (or Ctrl+C to skip)..."
python check_read_cooccurrence.py \
    --bam "$BAM_FILE" \
    --vcf "$VCF_FILE" \
    --min-qual 1000 \
    --min-depth 200 \
    --min-freq 0.03 \
    --output "${BASENAME}_cooccurrence_highqual.tsv"
echo -e "${GREEN}✓ Complete${NC}"
echo ""

# Example 3: Only very close variants (within read length)
echo "=================================="
echo "Example 3: Close Variants Only"
echo "=================================="
echo "Analyze only variants within read length (150bp)"
echo ""
echo "Command:"
echo "  python check_read_cooccurrence.py \\"
echo "    --bam $BAM_FILE \\"
echo "    --vcf $VCF_FILE \\"
echo "    --max-distance 150 \\"
echo "    --no-pairs \\"
echo "    --output ${BASENAME}_cooccurrence_close.tsv"
echo ""
read -p "Press Enter to run (or Ctrl+C to skip)..."
python check_read_cooccurrence.py \
    --bam "$BAM_FILE" \
    --vcf "$VCF_FILE" \
    --max-distance 150 \
    --no-pairs \
    --output "${BASENAME}_cooccurrence_close.tsv"
echo -e "${GREEN}✓ Complete${NC}"
echo ""

# Example 4: Verbose output
echo "=================================="
echo "Example 4: Verbose Output"
echo "=================================="
echo "Run with detailed logging"
echo ""
echo "Command:"
echo "  python check_read_cooccurrence.py \\"
echo "    --bam $BAM_FILE \\"
echo "    --vcf $VCF_FILE \\"
echo "    --verbose \\"
echo "    --output ${BASENAME}_cooccurrence_verbose.tsv"
echo ""
read -p "Press Enter to run (or Ctrl+C to skip)..."
python check_read_cooccurrence.py \
    --bam "$BAM_FILE" \
    --vcf "$VCF_FILE" \
    --verbose \
    --output "${BASENAME}_cooccurrence_verbose.tsv"
echo -e "${GREEN}✓ Complete${NC}"
echo ""

# Summary
echo "=================================="
echo "Summary"
echo "=================================="
echo -e "${GREEN}All examples completed successfully!${NC}"
echo ""
echo "Output files:"
ls -lh ${BASENAME}_cooccurrence_*.tsv 2>/dev/null || echo "No output files found"
echo ""
echo "To view results:"
echo "  column -t -s $'\\t' ${BASENAME}_cooccurrence_basic.tsv | less -S"
echo ""
echo "Next steps:"
echo "  1. Review co-occurrence results"
echo "  2. Compare with frequency-based haplotype predictions"
echo "  3. Look for proven linkage (linkage_proven = TRUE)"
echo "  4. Investigate any unexpected co-occurrence patterns"
echo ""
