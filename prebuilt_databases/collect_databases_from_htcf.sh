#!/bin/bash

# =============================================================================
# Collect and Package SnpEff Databases from HTCF
# Helper script to find, package, and prepare databases for distribution
# =============================================================================

set -euo pipefail

echo "=== VICAST Database Collection Tool ==="
echo ""

# Find SnpEff installation
if [[ -z "${SNPEFF_DATA:-}" ]]; then
    echo "SNPEFF_DATA not set. Searching for SnpEff installations..."

    # Common locations on HTCF
    POSSIBLE_LOCATIONS=(
        "/scratch/sahlab/kathie/snpEff/data"
        "/opt/snpEff/data"
        "$HOME/snpEff/data"
        "/ref/sahlab/software/snpEff/data"
    )

    for loc in "${POSSIBLE_LOCATIONS[@]}"; do
        if [[ -d "$loc" ]]; then
            echo "Found SnpEff data at: $loc"
            SNPEFF_DATA="$loc"
            break
        fi
    done

    if [[ -z "${SNPEFF_DATA:-}" ]]; then
        echo "Error: Could not find SnpEff data directory" >&2
        echo "Please set SNPEFF_DATA environment variable" >&2
        exit 1
    fi
fi

echo "Using SnpEff data directory: $SNPEFF_DATA"
echo ""

# Find all custom viral databases
echo "Scanning for viral genome databases..."
echo ""

# Look for directories with snpEffectPredictor.bin (built databases)
DATABASES=()
while IFS= read -r genome_dir; do
    genome=$(basename "$genome_dir")

    # Skip standard SnpEff databases (we want custom viral ones)
    # Custom viral genomes typically start with NC_, KY_, etc.
    if [[ $genome =~ ^(NC_|KY_|KU_|MN_|MK_|MT_|JX_|AF_|DQ_) ]]; then
        DATABASES+=("$genome")
    fi
done < <(find "$SNPEFF_DATA" -maxdepth 1 -type d -name "*" -exec test -f {}/snpEffectPredictor.bin \; -print)

if [[ ${#DATABASES[@]} -eq 0 ]]; then
    echo "No custom viral databases found in $SNPEFF_DATA"
    echo ""
    echo "Have you built any genomes with VICAST? Try:"
    echo "  step1_parse_viral_genome.py NC_001474.2"
    echo "  step2_add_to_snpeff.py NC_001474.2 NC_001474.2_features.tsv"
    exit 0
fi

# Display found databases
echo "Found ${#DATABASES[@]} custom viral database(s):"
echo ""
printf "%-20s %-15s %-10s\n" "ACCESSION" "SIZE" "FEATURES"
printf "%s\n" "--------------------------------------------------------"

for genome in "${DATABASES[@]}"; do
    size=$(du -sh "$SNPEFF_DATA/$genome" 2>/dev/null | cut -f1)

    # Count features (genes in genes.txt if it exists)
    features="?"
    if [[ -f "$SNPEFF_DATA/$genome/genes.txt" ]]; then
        features=$(wc -l < "$SNPEFF_DATA/$genome/genes.txt" | tr -d ' ')
    fi

    printf "%-20s %-15s %-10s\n" "$genome" "$size" "$features"
done

echo ""
echo "=== Next Steps ==="
echo ""
echo "To package these databases for distribution:"
echo ""
echo "1. Create output directory (use scratch space):"
echo "   mkdir -p /scratch/sahlab/kathie/vicast_prebuilt_dbs"
echo ""
echo "2. Package each database:"
for genome in "${DATABASES[@]}"; do
    echo "   bash prebuilt_databases/package_database.sh $genome /scratch/sahlab/kathie/vicast_prebuilt_dbs"
done
echo ""
echo "3. Download packaged files to your local machine:"
echo "   scp htcf:/scratch/sahlab/kathie/vicast_prebuilt_dbs/*.tar.gz /path/to/local/VICAST/prebuilt_databases/"
echo "   scp htcf:/scratch/sahlab/kathie/vicast_prebuilt_dbs/*.md5 /path/to/local/VICAST/prebuilt_databases/"
echo ""
echo "4. Update manifest.json with the MD5 checksums and metadata"
echo ""
echo "5. Commit and push to GitHub"
echo ""
echo "=== Quick Package All ==="
echo ""
echo "Run this to package all databases at once:"
echo ""
cat << 'EOF'
mkdir -p /scratch/sahlab/kathie/vicast_prebuilt_dbs
for genome in $(find $SNPEFF_DATA -maxdepth 1 -type d -name "NC_*" -exec test -f {}/snpEffectPredictor.bin \; -print | xargs -n1 basename); do
    echo "Packaging $genome..."
    bash prebuilt_databases/package_database.sh $genome /scratch/sahlab/kathie/vicast_prebuilt_dbs
done
EOF
