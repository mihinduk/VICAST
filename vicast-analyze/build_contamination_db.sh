#!/bin/bash

# =============================================================================
# Build VICAST Contamination Screening BLAST Database
# Downloads sequences from NCBI and builds a combined nucleotide database
# =============================================================================
#
# Reads contaminant organisms from contaminants.json config file.
# To add a new organism: edit contaminants.json and re-run this script.
#
# Use setup_blast_db.sh to install a pre-built version instead.
#
# Requirements: makeblastdb, efetch (entrez-direct), curl, python3
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

DB_NAME="microbial_contaminants"
EMAIL="${NCBI_EMAIL:-vicast@example.com}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_CONFIG="${SCRIPT_DIR}/blast_db/contaminants.json"

usage() {
    cat << 'EOF'
Usage: build_contamination_db.sh [OPTIONS] [OUTPUT_DIR]

Build the VICAST contamination screening BLAST database from NCBI sources.
Reads contaminant organisms from contaminants.json config file.

Arguments:
    OUTPUT_DIR          Directory for database output (default: current directory)

Options:
    --config FILE       Path to contaminants.json (default: auto-detect)
    --skip-viral        Skip viral RefSeq download (use existing viral_genomes.fasta in work dir)
    --viral-fasta FILE  Use pre-downloaded viral FASTA instead of downloading
    --package           Also create tar.gz for GitHub Release distribution
    --dry-run           Show what would be downloaded without actually downloading
    --help              Show this help message

Environment:
    NCBI_EMAIL          Email for NCBI Entrez queries (default: vicast@example.com)

Examples:
    # Full build
    build_contamination_db.sh /opt/vicast/blast_db

    # Build and package for distribution
    build_contamination_db.sh --package /tmp/blast_db_build

    # Dry run to verify config
    build_contamination_db.sh --dry-run

    # Use custom contaminants config
    build_contamination_db.sh --config my_contaminants.json /tmp/blast_db_build

    # Skip re-downloading viral genomes (use cached copy)
    build_contamination_db.sh --viral-fasta /tmp/viral_genomes.fasta /opt/vicast/blast_db

EOF
    exit 0
}

# Python helper to parse contaminants.json and download sequences
download_contaminants() {
    local config_file="$1"
    local output_fasta="$2"
    local dry_run="$3"

    python3 << PYEOF
import json
import subprocess
import sys
import os
import time

config_file = "$config_file"
output_fasta = "$output_fasta"
dry_run = "$dry_run" == "true"
email = "$EMAIL"

with open(config_file) as f:
    config = json.load(f)

contaminants = config.get("contaminants", [])
if not contaminants:
    print("No contaminants defined in config", file=sys.stderr)
    sys.exit(1)

print(f"Config: {len(contaminants)} contaminant organisms defined")
print()

if dry_run:
    for c in contaminants:
        acc_type = c.get("accession_type", "nucleotide")
        print(f"  [{c['category']}] {c['name']}")
        print(f"    Accession: {c['accession']} (type: {acc_type})")
        if c.get("notes"):
            print(f"    Notes: {c['notes']}")
    sys.exit(0)

total_seqs = 0
failed = []

with open(output_fasta, "w") as out:
    for c in contaminants:
        name = c["name"]
        accession = c["accession"]
        acc_type = c.get("accession_type", "nucleotide")
        category = c["category"]

        print(f"  Downloading {name} ({accession})...", end=" ", flush=True)

        if acc_type == "nucleotide":
            # Direct nucleotide accession - use efetch
            try:
                result = subprocess.run(
                    ["efetch", "-db", "nucleotide", "-id", accession,
                     "-format", "fasta"],
                    capture_output=True, text=True, timeout=120
                )
                if result.returncode == 0 and result.stdout.strip():
                    out.write(result.stdout)
                    if not result.stdout.endswith("\n"):
                        out.write("\n")
                    seq_count = result.stdout.count(">")
                    total_seqs += seq_count
                    print(f"OK ({seq_count} seq)")
                else:
                    print(f"FAILED")
                    failed.append(name)
            except Exception as e:
                print(f"ERROR: {e}")
                failed.append(name)

        elif acc_type == "assembly":
            # Assembly accession - need to resolve to nucleotide sequences
            # Use esearch + esummary to find linked nucleotide records
            try:
                # Method: search assembly DB, get linked nucleotide IDs
                result = subprocess.run(
                    ["esearch", "-db", "assembly", "-query", f"{accession}[Assembly Accession]"],
                    capture_output=True, text=True, timeout=60
                )
                if result.returncode != 0:
                    print(f"FAILED (esearch)")
                    failed.append(name)
                    continue

                # Link to nucleotide and fetch
                link_result = subprocess.run(
                    ["elink", "-target", "nucleotide"],
                    input=result.stdout, capture_output=True, text=True, timeout=60
                )
                if link_result.returncode != 0:
                    print(f"FAILED (elink)")
                    failed.append(name)
                    continue

                fetch_result = subprocess.run(
                    ["efetch", "-format", "fasta"],
                    input=link_result.stdout, capture_output=True, text=True, timeout=300
                )
                if fetch_result.returncode == 0 and fetch_result.stdout.strip():
                    out.write(fetch_result.stdout)
                    if not fetch_result.stdout.endswith("\n"):
                        out.write("\n")
                    seq_count = fetch_result.stdout.count(">")
                    total_seqs += seq_count
                    print(f"OK ({seq_count} seqs)")
                else:
                    print(f"FAILED (efetch)")
                    failed.append(name)
            except Exception as e:
                print(f"ERROR: {e}")
                failed.append(name)

        # Be polite to NCBI
        time.sleep(1)

print()
print(f"Total contaminant sequences: {total_seqs}")
if failed:
    print(f"Failed downloads: {', '.join(failed)}", file=sys.stderr)
    sys.exit(1 if total_seqs == 0 else 0)
PYEOF
}

# Main
main() {
    local output_dir="."
    local config_file="$DEFAULT_CONFIG"
    local skip_viral=false
    local viral_fasta=""
    local package=false
    local dry_run=false

    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h) usage ;;
            --config) config_file="$2"; shift 2 ;;
            --skip-viral) skip_viral=true; shift ;;
            --viral-fasta) viral_fasta="$2"; shift 2 ;;
            --package) package=true; shift ;;
            --dry-run) dry_run=true; shift ;;
            -*) echo -e "${RED}Unknown option: $1${NC}" >&2; usage ;;
            *) output_dir="$1"; shift ;;
        esac
    done

    # Validate config
    if [[ ! -f "$config_file" ]]; then
        echo -e "${RED}Config file not found: $config_file${NC}" >&2
        echo "Expected contaminants.json at: $DEFAULT_CONFIG" >&2
        exit 1
    fi

    echo "========================================="
    echo "Building VICAST Contamination Database"
    echo "========================================="
    echo "Config:  $config_file"
    echo "Output:  ${output_dir}/${DB_NAME}"
    echo "Email:   $EMAIL"
    echo ""

    if $dry_run; then
        echo -e "${BLUE}DRY RUN - showing what would be downloaded:${NC}"
        echo ""
        echo "Viral source:"
        python3 -c "import json; c=json.load(open('$config_file')); print(f\"  {c['viral_source']['url']}\")"
        echo ""
        echo "Contaminants:"
        download_contaminants "$config_file" "/dev/null" "true"
        exit 0
    fi

    mkdir -p "$output_dir"
    WORK_DIR=$(mktemp -d)
    trap "rm -rf $WORK_DIR" EXIT
    cd "$WORK_DIR"

    # -------------------------------------------------------------------------
    # Step 1: Viral genomes
    # -------------------------------------------------------------------------
    echo "Step 1: Viral genomes"
    echo "---------------------"

    if [[ -n "$viral_fasta" ]]; then
        echo "Using pre-downloaded viral FASTA: $viral_fasta"
        cp "$viral_fasta" viral_genomes.fasta
    elif $skip_viral; then
        echo "Skipping viral download (--skip-viral)"
        if [[ ! -f viral_genomes.fasta ]]; then
            echo -e "${RED}Error: viral_genomes.fasta not found in work dir${NC}" >&2
            exit 1
        fi
    else
        local viral_url
        viral_url=$(python3 -c "import json; print(json.load(open('$config_file'))['viral_source']['url'])")
        echo "Downloading from: $viral_url"
        curl -fSL --progress-bar "$viral_url" -o viral.genomic.fna.gz || {
            echo -e "${RED}Failed to download viral genomes${NC}" >&2
            exit 1
        }
        gunzip viral.genomic.fna.gz
        mv viral.genomic.fna viral_genomes.fasta
    fi

    VIRAL_COUNT=$(grep -c "^>" viral_genomes.fasta)
    echo "Viral sequences: $VIRAL_COUNT"
    echo ""

    # -------------------------------------------------------------------------
    # Step 2: Contaminant genomes (from config)
    # -------------------------------------------------------------------------
    echo "Step 2: Contaminant genomes"
    echo "---------------------------"

    download_contaminants "$config_file" "contaminants.fasta" "false"

    CONTAM_COUNT=$(grep -c "^>" contaminants.fasta 2>/dev/null || echo "0")
    echo ""

    # -------------------------------------------------------------------------
    # Step 3: Build BLAST database
    # -------------------------------------------------------------------------
    echo "Step 3: Building BLAST database"
    echo "-------------------------------"

    cat viral_genomes.fasta contaminants.fasta > "${DB_NAME}.fasta"
    TOTAL=$(grep -c "^>" "${DB_NAME}.fasta")
    echo "Total sequences: $TOTAL"

    makeblastdb \
        -in "${DB_NAME}.fasta" \
        -dbtype nucl \
        -out "${output_dir}/${DB_NAME}" \
        -title "VICAST Contamination Screening Database" \
        -parse_seqids 2>&1

    # -------------------------------------------------------------------------
    # Step 4: Create info files
    # -------------------------------------------------------------------------

    # JSON info (machine-readable)
    python3 << PYEOF
import json
from datetime import datetime

config = json.load(open("$config_file"))
info = {
    "database_name": "${DB_NAME}",
    "created": datetime.now().isoformat(),
    "total_sequences": $TOTAL,
    "viral_sequences": $VIRAL_COUNT,
    "contaminant_sequences": $CONTAM_COUNT,
    "contaminants": config["contaminants"],
    "viral_source": config["viral_source"]["url"],
    "config_version": config["version"]
}
with open("${output_dir}/${DB_NAME}_info.json", "w") as f:
    json.dump(info, f, indent=2)
PYEOF

    # Text info (human-readable, backward compatible)
    cat > "${output_dir}/${DB_NAME}_info.txt" << EOF
VICAST Contamination Screening BLAST Database
=====================================================
Created: $(date)
Total sequences: $TOTAL
Viral sequences: $VIRAL_COUNT
Contaminant sequences: $CONTAM_COUNT

Organisms successfully included:
EOF

    python3 -c "
import json
config = json.load(open('$config_file'))
for c in config['contaminants']:
    print(f\"  - {c['name']} ({c['accession']}): [{c['category']}]\")
" >> "${output_dir}/${DB_NAME}_info.txt"

    cat >> "${output_dir}/${DB_NAME}_info.txt" << EOF

Database location: ${output_dir}/${DB_NAME}

To add new organisms: edit contaminants.json and re-run build_contamination_db.sh
To add user-specific sequences: use extend_blast_db.sh
EOF

    # Create alias
    echo ""
    echo "Creating vicast_combined alias..."
    if command -v blastdb_aliastool &> /dev/null; then
        (cd "$output_dir" && blastdb_aliastool -dblist "${DB_NAME}" \
            -dbtype nucl -out vicast_combined \
            -title "VICAST Combined BLAST DB") 2>/dev/null || true
    fi

    echo ""
    echo -e "${GREEN}Database built successfully!${NC}"
    echo "Location: ${output_dir}/${DB_NAME}"
    echo "Total sequences: $TOTAL ($VIRAL_COUNT viral + $CONTAM_COUNT contaminants)"

    # -------------------------------------------------------------------------
    # Step 5: Package for distribution (optional)
    # -------------------------------------------------------------------------
    if $package; then
        echo ""
        echo "Step 5: Packaging for distribution"
        echo "----------------------------------"
        local tarball="${output_dir}/${DB_NAME}.tar.gz"
        (cd "$output_dir" && tar czf "${DB_NAME}.tar.gz" \
            ${DB_NAME}.nhr ${DB_NAME}.nin ${DB_NAME}.nog \
            ${DB_NAME}.nsd ${DB_NAME}.nsi ${DB_NAME}.nsq \
            ${DB_NAME}_info.txt ${DB_NAME}_info.json 2>/dev/null) || \
        (cd "$output_dir" && tar czf "${DB_NAME}.tar.gz" ${DB_NAME}.*)
        local size
        size=$(ls -lh "$tarball" | awk '{print $5}')
        echo -e "${GREEN}Package created: $tarball ($size)${NC}"
        echo ""
        echo "To upload as GitHub Release:"
        echo "  gh release create blast-db-v1.0 $tarball \\"
        echo "      --repo mihinduk/VICAST \\"
        echo "      --title 'BLAST Contamination Database v1.0' \\"
        echo "      --notes 'Pre-built database: $TOTAL sequences'"
    fi
}

main "$@"
