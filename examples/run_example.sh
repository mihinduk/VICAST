#!/usr/bin/env bash
#
# VICAST Example Workflow
#
# This script demonstrates the VICAST pipeline using synthetic test data.
# It shows both the annotation and analysis components.
#
# Prerequisites:
#   - VICAST installed (pip install -e .)
#   - SnpEff installed and configured (SNPEFF_HOME set)
#   - Python 3.9+
#
# Usage:
#   ./run_example.sh [--skip-snpeff]
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTPUT_DIR="${SCRIPT_DIR}/output"
EXPECTED_DIR="${SCRIPT_DIR}/expected_output"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse arguments
SKIP_SNPEFF=false
for arg in "$@"; do
    case $arg in
        --skip-snpeff)
            SKIP_SNPEFF=true
            shift
            ;;
    esac
done

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=============================================="
echo "VICAST Example Workflow"
echo "=============================================="
echo ""

# Step 1: Validate GFF file
log_info "Step 1: Validating GFF annotation file..."

python3 -c "
import sys
sys.path.insert(0, '${SCRIPT_DIR}/../src')
from vicast.validation import validate_gff_for_snpeff

gff_file = '${DATA_DIR}/test_virus.gff3'
fasta_file = '${DATA_DIR}/test_virus.fasta'

is_valid, errors, warnings = validate_gff_for_snpeff(gff_file, fasta_file)

print(f'  GFF Valid: {is_valid}')
print(f'  Errors: {len(errors)}')
print(f'  Warnings: {len(warnings)}')

if errors:
    for e in errors:
        print(f'    - {e}')
    sys.exit(1)

if warnings:
    for w in warnings:
        print(f'    (warning) {w}')
"

log_info "GFF validation passed!"
echo ""

# Step 2: Convert annotation TSV to GFF
log_info "Step 2: Converting annotation TSV to GFF format..."

python3 -c "
import sys
sys.path.insert(0, '${SCRIPT_DIR}/../vicast-annotate')
from step2_add_to_snpeff import tsv_to_gff

tsv_file = '${DATA_DIR}/test_virus_annotations.tsv'
output_gff = '${OUTPUT_DIR}/converted_annotations.gff3'

success, count, summary = tsv_to_gff(tsv_file, output_gff, force=True)

print(f'  Converted: {count} features')
print(f'  Output: {output_gff}')

if not success:
    print('  Conversion failed!')
    sys.exit(1)
"

log_info "TSV to GFF conversion complete!"
echo ""

# Step 3: Parse and filter variants
log_info "Step 3: Parsing VCF and filtering variants..."

python3 -c "
import sys
import re
import pandas as pd

vcf_file = '${DATA_DIR}/test_virus_variants.vcf'
output_tsv = '${OUTPUT_DIR}/filtered_variants.tsv'

# Parse VCF
variants = []
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 8:
            info = parts[7]

            # Extract INFO fields
            af_match = re.search(r'AF=([0-9.]+)', info)
            dp_match = re.search(r'DP=([0-9]+)', info)

            variants.append({
                'CHROM': parts[0],
                'POS': int(parts[1]),
                'REF': parts[3],
                'ALT': parts[4],
                'QUAL': float(parts[5]),
                'INFO': info,
                'AF': float(af_match.group(1)) if af_match else 0.0,
                'DP': int(dp_match.group(1)) if dp_match else 0
            })

df = pd.DataFrame(variants)
print(f'  Total variants: {len(df)}')

# Apply filters
quality_cutoff = 1000
depth_cutoff = 200
freq_cutoff = 0.10

filtered = df[
    (df['QUAL'] >= quality_cutoff) &
    (df['DP'] >= depth_cutoff) &
    (df['AF'] >= freq_cutoff)
]

print(f'  After quality filter (>={quality_cutoff}): {len(df[df[\"QUAL\"] >= quality_cutoff])}')
print(f'  After depth filter (>={depth_cutoff}): {len(df[df[\"DP\"] >= depth_cutoff])}')
print(f'  After frequency filter (>={freq_cutoff}): {len(df[df[\"AF\"] >= freq_cutoff])}')
print(f'  Variants passing all filters: {len(filtered)}')

# Save filtered variants
filtered.to_csv(output_tsv, sep='\t', index=False)
print(f'  Output: {output_tsv}')
"

log_info "Variant filtering complete!"
echo ""

# Step 4: Generate summary statistics
log_info "Step 4: Generating summary statistics..."

python3 -c "
import pandas as pd

variants_file = '${OUTPUT_DIR}/filtered_variants.tsv'
df = pd.read_csv(variants_file, sep='\t')

print('  === Variant Summary ===')
print(f'  Total filtered variants: {len(df)}')
print(f'  Mean allele frequency: {df[\"AF\"].mean():.3f}')
print(f'  Mean depth: {df[\"DP\"].mean():.1f}')
print(f'  Mean quality: {df[\"QUAL\"].mean():.1f}')
print('')
print('  Positions:')
for _, row in df.iterrows():
    print(f'    {row[\"CHROM\"]}:{row[\"POS\"]} {row[\"REF\"]}>{row[\"ALT\"]} (AF={row[\"AF\"]:.2f}, DP={row[\"DP\"]})')
"

echo ""

# Optional: SnpEff annotation
if [ "$SKIP_SNPEFF" = false ]; then
    if [ -n "${SNPEFF_HOME}" ] && [ -f "${SNPEFF_HOME}/snpEff.jar" ]; then
        log_info "Step 5: SnpEff annotation available"
        log_info "To add this genome to SnpEff, run:"
        echo "    python vicast-annotate/step2_add_to_snpeff.py ${DATA_DIR}/test_virus.fasta ${DATA_DIR}/test_virus_annotations.tsv --force"
    else
        log_warn "SnpEff not configured (SNPEFF_HOME not set or snpEff.jar not found)"
        log_info "Skipping SnpEff integration steps"
    fi
else
    log_info "Step 5: Skipped (--skip-snpeff flag)"
fi

echo ""
echo "=============================================="
log_info "Example workflow complete!"
echo "=============================================="
echo ""
echo "Output files:"
ls -la "${OUTPUT_DIR}/"
