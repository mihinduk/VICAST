#!/bin/bash
################################################################################
# Validation Script for Read Co-occurrence Analysis Module
#
# This script validates that all dependencies are installed and the module
# is ready to use.
################################################################################

set -e

# Colors
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "=========================================="
echo "Read Co-occurrence Module Validation"
echo "=========================================="
echo ""

# Check Python
echo -e "${BLUE}Checking Python...${NC}"
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version)
    echo -e "${GREEN}✓ Python found: $PYTHON_VERSION${NC}"
else
    echo -e "${RED}✗ Python 3 not found${NC}"
    exit 1
fi
echo ""

# Check pysam
echo -e "${BLUE}Checking pysam...${NC}"
if python3 -c "import pysam; print(f'pysam version: {pysam.__version__}')" 2>/dev/null; then
    PYSAM_VERSION=$(python3 -c "import pysam; print(pysam.__version__)")
    echo -e "${GREEN}✓ pysam found: version $PYSAM_VERSION${NC}"
else
    echo -e "${RED}✗ pysam not found${NC}"
    echo -e "${YELLOW}  Install with: conda install -c bioconda pysam${NC}"
    echo -e "${YELLOW}  Or: pip install pysam${NC}"
    exit 1
fi
echo ""

# Check pandas
echo -e "${BLUE}Checking pandas...${NC}"
if python3 -c "import pandas; print(f'pandas version: {pandas.__version__}')" 2>/dev/null; then
    PANDAS_VERSION=$(python3 -c "import pandas; print(pandas.__version__)")
    echo -e "${GREEN}✓ pandas found: version $PANDAS_VERSION${NC}"
else
    echo -e "${RED}✗ pandas not found${NC}"
    echo -e "${YELLOW}  Install with: conda install pandas${NC}"
    exit 1
fi
echo ""

# Check samtools (optional, but recommended)
echo -e "${BLUE}Checking samtools (optional)...${NC}"
if command -v samtools &> /dev/null; then
    SAMTOOLS_VERSION=$(samtools --version | head -1)
    echo -e "${GREEN}✓ samtools found: $SAMTOOLS_VERSION${NC}"
else
    echo -e "${YELLOW}⚠ samtools not found (recommended for BAM indexing)${NC}"
    echo -e "${YELLOW}  Install with: conda install -c bioconda samtools${NC}"
fi
echo ""

# Check if module files exist
echo -e "${BLUE}Checking module files...${NC}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

FILES=(
    "check_read_cooccurrence.py"
    "test_read_cooccurrence.py"
    "READ_COOCCURRENCE_GUIDE.md"
    "README_READ_COOCCURRENCE.md"
    "example_cooccurrence_usage.sh"
)

for FILE in "${FILES[@]}"; do
    if [ -f "$SCRIPT_DIR/$FILE" ]; then
        echo -e "${GREEN}✓ $FILE${NC}"
    else
        echo -e "${RED}✗ $FILE missing${NC}"
        exit 1
    fi
done
echo ""

# Test import
echo -e "${BLUE}Testing module import...${NC}"
cd "$SCRIPT_DIR"
if python3 -c "from check_read_cooccurrence import VariantCooccurrenceAnalyzer; print('Import successful')" 2>/dev/null; then
    echo -e "${GREEN}✓ Module imports successfully${NC}"
else
    echo -e "${RED}✗ Module import failed${NC}"
    exit 1
fi
echo ""

# Test help
echo -e "${BLUE}Testing command-line interface...${NC}"
if python3 check_read_cooccurrence.py --help > /dev/null 2>&1; then
    echo -e "${GREEN}✓ Command-line interface works${NC}"
else
    echo -e "${RED}✗ Command-line interface failed${NC}"
    exit 1
fi
echo ""

# Summary
echo "=========================================="
echo -e "${GREEN}✓ All validation checks passed!${NC}"
echo "=========================================="
echo ""
echo "The read co-occurrence analysis module is ready to use."
echo ""
echo "Quick start:"
echo "  python3 check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf"
echo ""
echo "Run tests:"
echo "  python3 test_read_cooccurrence.py"
echo ""
echo "Documentation:"
echo "  - Quick start: README_READ_COOCCURRENCE.md"
echo "  - Full guide: READ_COOCCURRENCE_GUIDE.md"
echo "  - Examples: example_cooccurrence_usage.sh"
echo ""
