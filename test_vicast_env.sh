#!/bin/bash
# Test script to verify VICAST environment is working correctly

echo "=========================================="
echo "Testing VICAST Environment"
echo "=========================================="
echo ""

# Check we're in the right environment
if [ "$CONDA_DEFAULT_ENV" != "vicast" ]; then
    echo "ERROR: Not in vicast environment"
    echo "Current environment: $CONDA_DEFAULT_ENV"
    echo ""
    echo "Please activate first:"
    echo "  source $CONDA_BASE/bin/activate (or your conda installation)"
    echo "  conda activate vicast"
    exit 1
fi

echo "✓ Environment: $CONDA_DEFAULT_ENV"
echo ""

# Test 1: Python packages
echo "Test 1: Python Package Imports"
echo "-------------------------------"
python3 << 'EOF'
import sys
try:
    import Bio
    print(f"✓ Biopython {Bio.__version__}")
except ImportError as e:
    print(f"✗ Biopython import failed: {e}")
    sys.exit(1)

try:
    import pandas
    print(f"✓ Pandas {pandas.__version__}")
except ImportError as e:
    print(f"✗ Pandas import failed: {e}")
    sys.exit(1)

try:
    import numpy
    print(f"✓ NumPy {numpy.__version__}")
except ImportError as e:
    print(f"✗ NumPy import failed: {e}")
    sys.exit(1)
EOF

if [ $? -ne 0 ]; then
    echo ""
    echo "✗ Python package test FAILED"
    exit 1
fi

echo ""

# Test 2: BLAST
echo "Test 2: BLAST Installation"
echo "---------------------------"
if command -v blastx &> /dev/null; then
    echo "✓ blastx found"
    blastx -version | head -1
else
    echo "✗ blastx not found"
    exit 1
fi

echo ""

# Test 3: VICAST scripts
echo "Test 3: VICAST Script Execution"
echo "--------------------------------"

# Get script directory
SCRIPT_DIR="$(dirname "$0")"

# Test step0
echo "Testing step0_check_snpeff.py..."
python3 "$SCRIPT_DIR/vicast-annotate/step0_check_snpeff.py" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ step0_check_snpeff.py loads"
else
    echo "✗ step0_check_snpeff.py failed"
    exit 1
fi

# Test step1
echo "Testing step1_parse_viral_genome.py..."
python3 "$SCRIPT_DIR/vicast-annotate/step1_parse_viral_genome.py" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ step1_parse_viral_genome.py loads"
else
    echo "✗ step1_parse_viral_genome.py failed"
    exit 1
fi

# Test step1_blastx
echo "Testing step1_blastx_annotate.py..."
python3 "$SCRIPT_DIR/vicast-annotate/step1_blastx_annotate.py" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ step1_blastx_annotate.py loads"
else
    echo "✗ step1_blastx_annotate.py failed"
    exit 1
fi

# Test step2
echo "Testing step2_add_to_snpeff.py..."
python3 "$SCRIPT_DIR/vicast-annotate/step2_add_to_snpeff.py" --help > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo "✓ step2_add_to_snpeff.py loads"
else
    echo "✗ step2_add_to_snpeff.py failed"
    exit 1
fi

echo ""

# Test 4: Quick functional test
echo "Test 4: Pathway Detection (Quick Functional Test)"
echo "--------------------------------------------------"
cd /tmp
python3 "$SCRIPT_DIR/vicast-annotate/step0_check_snpeff.py" NC_045512 2>&1 | head -10

echo ""
echo "=========================================="
echo "All Tests Passed! ✓"
echo "=========================================="
echo ""
echo "VICAST environment is ready to use."
echo ""
echo "Next steps:"
echo "1. Test Pathway 2 (well-annotated genomes)"
echo "2. Test Pathway 3 (BLASTx) with a BLAST database"
echo "3. Test Pathway 4 (segmented viruses)"
