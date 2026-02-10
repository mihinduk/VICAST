#!/bin/bash

# =============================================================================
# Package SnpEff Database for Distribution
# Creates tar.gz archive of existing SnpEff database
# =============================================================================

set -euo pipefail

usage() {
    cat << EOF
Usage: $(basename "$0") ACCESSION [OUTPUT_DIR]

Package an existing SnpEff database for distribution.

Arguments:
    ACCESSION   Genome accession (e.g., NC_001474.2)
    OUTPUT_DIR  Output directory (default: current directory)

Environment:
    SNPEFF_DATA SnpEff data directory (default: \$SNPEFF_HOME/data)

Example:
    # Package DENV-2 database
    $(basename "$0") NC_001474.2

    # Package to specific directory
    $(basename "$0") NC_001474.2 /path/to/output

Output:
    Creates: {ACCESSION}.tar.gz and {ACCESSION}.md5

EOF
    exit 1
}

if [[ $# -lt 1 ]]; then
    usage
fi

ACCESSION="$1"
OUTPUT_DIR="${2:-.}"
SNPEFF_DATA="${SNPEFF_DATA:-${SNPEFF_HOME}/data}"

# Check if database exists
if [[ ! -d "$SNPEFF_DATA/$ACCESSION" ]]; then
    echo "Error: Database $ACCESSION not found in $SNPEFF_DATA" >&2
    exit 1
fi

if [[ ! -f "$SNPEFF_DATA/$ACCESSION/snpEffectPredictor.bin" ]]; then
    echo "Error: $ACCESSION appears incomplete (missing snpEffectPredictor.bin)" >&2
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Package database
echo "Packaging $ACCESSION..."
OUTPUT_FILE="$OUTPUT_DIR/${ACCESSION}.tar.gz"

cd "$SNPEFF_DATA"
tar -czf "$OUTPUT_FILE" "$ACCESSION"

# Calculate MD5
echo "Calculating checksum..."
if command -v md5sum &> /dev/null; then
    md5sum "$OUTPUT_FILE" > "${OUTPUT_FILE}.md5"
    MD5=$(md5sum "$OUTPUT_FILE" | cut -d' ' -f1)
elif command -v md5 &> /dev/null; then
    md5 -q "$OUTPUT_FILE" > "${OUTPUT_FILE}.md5"
    MD5=$(md5 -q "$OUTPUT_FILE")
fi

# Get size
SIZE_BYTES=$(stat -f%z "$OUTPUT_FILE" 2>/dev/null || stat -c%s "$OUTPUT_FILE")
SIZE_MB=$(echo "scale=1; $SIZE_BYTES / 1048576" | bc)

echo ""
echo "✓ Package created successfully"
echo ""
echo "File: $OUTPUT_FILE"
echo "Size: ${SIZE_MB} MB"
echo "MD5:  $MD5"
echo ""
echo "Add to manifest.json:"
cat << EOF_JSON
    {
      "accession": "$ACCESSION",
      "name": "VIRUS_NAME",
      "description": "Virus description",
      "size_mb": $SIZE_MB,
      "features": NUM_FEATURES,
      "tags": ["tag1", "tag2"],
      "url": "https://github.com/shandley/VICAST/raw/main/prebuilt_databases/${ACCESSION}.tar.gz",
      "md5": "$MD5"
    }
EOF_JSON
