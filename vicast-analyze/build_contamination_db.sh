#!/bin/bash

# =============================================================================
# Build VICAST Contamination Screening BLAST Database
# Downloads sequences from NCBI and builds a combined nucleotide database
# =============================================================================
#
# This script recreates the contamination database from NCBI sources.
# Use setup_blast_db.sh to install the pre-built version instead.
#
# Requirements: makeblastdb, datasets (NCBI), curl/wget
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

OUTPUT_DIR="${1:-.}"
DB_NAME="microbial_contaminants"
EMAIL="${NCBI_EMAIL:-vicast@example.com}"

echo "========================================="
echo "Building VICAST Contamination Database"
echo "========================================="
echo "Output: ${OUTPUT_DIR}/${DB_NAME}"
echo ""

mkdir -p "$OUTPUT_DIR"
WORK_DIR=$(mktemp -d)
trap "rm -rf $WORK_DIR" EXIT

cd "$WORK_DIR"

# -------------------------------------------------------------------------
# 1. Download NCBI RefSeq viral genomes
# -------------------------------------------------------------------------
echo "Step 1: Downloading NCBI RefSeq viral genomes..."

if command -v datasets &> /dev/null; then
    # Use NCBI datasets CLI (preferred)
    datasets download genome taxon "Viruses" \
        --reference \
        --include genome \
        --filename viral_genomes.zip 2>&1

    unzip -o viral_genomes.zip -d viral_download
    cat viral_download/ncbi_dataset/data/*/genomic.fna > viral_genomes.fasta
    rm -rf viral_download viral_genomes.zip
else
    # Fallback: download from NCBI FTP
    echo "  (datasets CLI not found, using FTP download)"
    curl -fSL "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz" \
        -o viral.1.1.genomic.fna.gz || {
        echo -e "${RED}Failed to download viral genomes${NC}" >&2
        exit 1
    }
    gunzip viral.1.1.genomic.fna.gz
    mv viral.1.1.genomic.fna viral_genomes.fasta
fi

VIRAL_COUNT=$(grep -c "^>" viral_genomes.fasta)
echo "  Downloaded $VIRAL_COUNT viral sequences"

# -------------------------------------------------------------------------
# 2. Download common lab contaminant genomes
# -------------------------------------------------------------------------
echo ""
echo "Step 2: Downloading common lab contaminant genomes..."

# Contaminants to include (accession: organism)
declare -A CONTAMINANTS=(
    ["GCF_000146045.2"]="Saccharomyces_cerevisiae_S288C"
    ["GCF_000182965.3"]="Candida_albicans_SC5314"
    ["GCF_000149245.2"]="Cryptococcus_neoformans_JEC21"
    ["GCF_000005845.2"]="Escherichia_coli_K12_MG1655"
    ["GCF_000006765.1"]="Pseudomonas_aeruginosa_PAO1"
    ["GCF_000013425.1"]="Staphylococcus_aureus_NCTC8325"
)

# Mycoplasma hyorhinis - use direct NCBI accession
MYCOPLASMA_ACC="NZ_CP038034.1"

> contaminants.fasta

for acc in "${!CONTAMINANTS[@]}"; do
    name="${CONTAMINANTS[$acc]}"
    echo "  Downloading $name ($acc)..."

    if command -v datasets &> /dev/null; then
        datasets download genome accession "$acc" --include genome --filename "${acc}.zip" 2>/dev/null
        unzip -o "${acc}.zip" -d "${acc}_dir" 2>/dev/null
        cat "${acc}_dir"/ncbi_dataset/data/*/genomic.fna >> contaminants.fasta 2>/dev/null || true
        rm -rf "${acc}.zip" "${acc}_dir"
    else
        # Fallback: use efetch
        python3 -c "
from Bio import Entrez, SeqIO
import sys
Entrez.email = '${EMAIL}'
try:
    handle = Entrez.esearch(db='assembly', term='${acc}[Assembly Accession]')
    results = Entrez.read(handle)
    handle.close()
    if results['IdList']:
        handle = Entrez.efetch(db='nucleotide', id='${acc}', rettype='fasta', retmode='text')
        with open('temp_contaminant.fasta', 'w') as f:
            f.write(handle.read())
        handle.close()
except Exception as e:
    print(f'Warning: Could not download ${name}: {e}', file=sys.stderr)
" 2>/dev/null
        if [[ -f temp_contaminant.fasta ]]; then
            cat temp_contaminant.fasta >> contaminants.fasta
            rm temp_contaminant.fasta
        fi
    fi

    count=$(grep -c "^>" contaminants.fasta 2>/dev/null || echo "0")
    echo "    $name: included ($count total contaminant sequences so far)"
done

# Download Mycoplasma hyorhinis
echo "  Downloading Mycoplasma_hyorhinis ($MYCOPLASMA_ACC)..."
python3 -c "
from Bio import Entrez, SeqIO
Entrez.email = '${EMAIL}'
handle = Entrez.efetch(db='nucleotide', id='${MYCOPLASMA_ACC}', rettype='fasta', retmode='text')
with open('mycoplasma.fasta', 'w') as f:
    f.write(handle.read())
handle.close()
" 2>/dev/null || echo "  Warning: Could not download Mycoplasma"

if [[ -f mycoplasma.fasta ]]; then
    cat mycoplasma.fasta >> contaminants.fasta
    rm mycoplasma.fasta
fi

CONTAM_COUNT=$(grep -c "^>" contaminants.fasta 2>/dev/null || echo "0")
echo "  Total contaminant sequences: $CONTAM_COUNT"

# -------------------------------------------------------------------------
# 3. Combine and build BLAST database
# -------------------------------------------------------------------------
echo ""
echo "Step 3: Building BLAST database..."

cat viral_genomes.fasta contaminants.fasta > "${DB_NAME}.fasta"
TOTAL=$(grep -c "^>" "${DB_NAME}.fasta")
echo "  Total sequences: $TOTAL"

makeblastdb \
    -in "${DB_NAME}.fasta" \
    -dbtype nucl \
    -out "${OUTPUT_DIR}/${DB_NAME}" \
    -title "Microbial Contaminants Database (No Mammals)" \
    -parse_seqids 2>&1

# -------------------------------------------------------------------------
# 4. Create info file
# -------------------------------------------------------------------------
cat > "${OUTPUT_DIR}/${DB_NAME}_info.txt" << EOF
Microbial Contamination BLAST Database (No Mammals)
=====================================================
Created: $(date)
Total sequences: $TOTAL

Organisms included:
EOF

for acc in "${!CONTAMINANTS[@]}"; do
    name="${CONTAMINANTS[$acc]}"
    count=$(grep -c "$name\|${acc}" contaminants.fasta 2>/dev/null || echo "?")
    echo "  - ${name}: ${count} sequences" >> "${OUTPUT_DIR}/${DB_NAME}_info.txt"
done
echo "  - Mycoplasma_hyorhinis: 1 sequences" >> "${OUTPUT_DIR}/${DB_NAME}_info.txt"
echo "  - viral genomes: ${VIRAL_COUNT} sequences" >> "${OUTPUT_DIR}/${DB_NAME}_info.txt"

cat >> "${OUTPUT_DIR}/${DB_NAME}_info.txt" << EOF

Database location: ${OUTPUT_DIR}/${DB_NAME}
Database name: ${DB_NAME}

Note: This database excludes mammalian genomes to reduce size and build time.

Usage:
blastn -query your_contigs.fa -db ${OUTPUT_DIR}/${DB_NAME} -outfmt 6 -max_target_seqs 5 -evalue 1e-10 -num_threads 4
EOF

echo ""
echo -e "${GREEN}Database built successfully!${NC}"
echo "Location: ${OUTPUT_DIR}/${DB_NAME}"
echo "Total sequences: $TOTAL"
echo ""
echo "To package for distribution:"
echo "  cd ${OUTPUT_DIR} && tar czf ${DB_NAME}.tar.gz ${DB_NAME}.*"
