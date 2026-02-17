#!/bin/bash
# =============================================================================
# test_blast_diagnostic.sh - Standalone BLAST test for viral_diagnostic.sh
#
# Runs ONLY the BLAST analysis + top hits parsing + HTML report generation
# on pre-assembled contigs, WITHOUT re-running mapping or assembly.
#
# Usage:
#   test_blast_diagnostic.sh <contigs_fasta> <accession> [threads]
#
# Example (Docker):
#   docker run --rm --user $(id -u):$(id -g) \
#       -v /mnt/pathogen2/kathie/vicast_analysis:/data \
#       -w /data vicast:latest \
#       bash /opt/vicast/vicast-analyze/test_blast_diagnostic.sh \
#         diagnostic_SRR5992153/SRR5992153_contigs_filtered.fa NC_001474.2 8
#
# This script is for iterative debugging of the BLAST portion of step 6.
# =============================================================================

set -euo pipefail

# Parse arguments
CONTIGS_FASTA="${1:?Usage: test_blast_diagnostic.sh <contigs_fasta> <accession> [threads]}"
ACCESSION="${2:?Usage: test_blast_diagnostic.sh <contigs_fasta> <accession> [threads]}"
THREADS="${3:-4}"

# Validate inputs
if [ ! -f "$CONTIGS_FASTA" ]; then
    echo "ERROR: Contigs FASTA not found: $CONTIGS_FASTA"
    exit 1
fi

# Derive sample name from contigs filename
# e.g., SRR5992153_contigs_filtered.fa -> SRR5992153
BASENAME=$(basename "$CONTIGS_FASTA")
SAMPLE_NAME="${BASENAME%%_contigs*}"
if [ "$SAMPLE_NAME" = "$BASENAME" ]; then
    # Fallback: strip extension
    SAMPLE_NAME="${BASENAME%.*}"
fi

# Resolve to absolute path before cd
CONTIGS_FASTA="$(readlink -f "$CONTIGS_FASTA")"

# Work in the same directory as the contigs file
WORK_DIR=$(dirname "$CONTIGS_FASTA")
cd "$WORK_DIR"

# Find pipeline directory (where this script and parse_blast_results.py live)
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "========================================="
echo "BLAST DIAGNOSTIC TEST"
echo "========================================="
echo "Contigs:   $CONTIGS_FASTA"
echo "Accession: $ACCESSION"
echo "Threads:   $THREADS"
echo "Work dir:  $WORK_DIR"
echo "Sample:    $SAMPLE_NAME"
echo "Pipeline:  $PIPELINE_DIR"
echo ""

# Count contigs
FILTERED_COUNT=$(grep -c ">" "$CONTIGS_FASTA" || echo 0)
echo "Contigs in FASTA: $FILTERED_COUNT"
echo ""

if [ "$FILTERED_COUNT" -eq 0 ]; then
    echo "ERROR: No contigs found in $CONTIGS_FASTA"
    exit 1
fi

# =============================================================================
# Step 1: BLAST search
# =============================================================================
echo "========================================="
echo "Step 1: Running BLAST"
echo "========================================="

# Check for local contamination database
LOCAL_DB_PATH="${BLAST_DB:-/opt/vicast/blast_db/microbial_contaminants}"

# Database diagnostics
echo "Database path: $LOCAL_DB_PATH"
if [ -f "${LOCAL_DB_PATH}.nhr" ] || [ -f "${LOCAL_DB_PATH}.00.nhr" ]; then
    echo "Database found. Info:"
    blastdbcmd -db "$LOCAL_DB_PATH" -info 2>&1 | head -5
    echo ""
else
    echo "WARNING: Local database NOT found at $LOCAL_DB_PATH"
    echo "Checking available files:"
    ls -la "$(dirname "$LOCAL_DB_PATH")"/ 2>/dev/null || echo "  Directory does not exist"
    echo ""
fi

# Add header to BLAST results
echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tMismatches\tGap_Opens\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tSubject_Title" > "${SAMPLE_NAME}_blast_all.tsv"

if [ -f "${LOCAL_DB_PATH}.nhr" ] || [ -f "${LOCAL_DB_PATH}.00.nhr" ]; then
    echo "Running local BLAST with $THREADS threads..."
    blastn -query "$CONTIGS_FASTA" \
           -db "$LOCAL_DB_PATH" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
           -max_target_seqs 5 \
           -max_hsps 1 \
           -evalue 1e-10 \
           -num_threads "$THREADS" \
           >> "${SAMPLE_NAME}_blast_all.tsv"
    BLAST_EXIT=$?
    echo "BLAST exit code: $BLAST_EXIT"
else
    echo "Running remote BLAST (limited to top 20 contigs)..."
    head -40 "$CONTIGS_FASTA" > "${SAMPLE_NAME}_top20_contigs.fa"
    blastn -query "${SAMPLE_NAME}_top20_contigs.fa" \
           -remote -db nt \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
           -max_target_seqs 5 \
           -max_hsps 1 \
           -evalue 1e-10 \
           >> "${SAMPLE_NAME}_blast_all.tsv"
    BLAST_EXIT=$?
    echo "BLAST exit code: $BLAST_EXIT"
    rm -f "${SAMPLE_NAME}_top20_contigs.fa"
fi

# Count results
TOTAL_HITS=$(tail -n +2 "${SAMPLE_NAME}_blast_all.tsv" | wc -l | tr -d ' ')
CONTIGS_WITH_HITS=$(tail -n +2 "${SAMPLE_NAME}_blast_all.tsv" | cut -f1 | sort -u | wc -l | tr -d ' ')
echo ""
echo "BLAST results: $TOTAL_HITS hits for $CONTIGS_WITH_HITS / $FILTERED_COUNT contigs"

# Show which contigs got hits vs. which didn't
echo ""
echo "Contigs WITH hits:"
tail -n +2 "${SAMPLE_NAME}_blast_all.tsv" | cut -f1 | sort -u | while read c; do
    echo "  $c"
done

echo ""
echo "Contigs WITHOUT hits:"
grep ">" "$CONTIGS_FASTA" | sed 's/>//;s/ .*//' | while read c; do
    if ! tail -n +2 "${SAMPLE_NAME}_blast_all.tsv" | cut -f1 | grep -q "^${c}$"; then
        echo "  $c"
    fi
done

# =============================================================================
# Step 2: Parse BLAST results with Python (tied hits + accession matching)
# =============================================================================
echo ""
echo "========================================="
echo "Step 2: Parsing BLAST results"
echo "========================================="

PARSE_SCRIPT="${PIPELINE_DIR}/parse_blast_results.py"
if [ ! -f "$PARSE_SCRIPT" ]; then
    echo "ERROR: parse_blast_results.py not found at: $PARSE_SCRIPT"
    exit 1
fi

python "$PARSE_SCRIPT" \
    "${SAMPLE_NAME}_blast_all.tsv" \
    "$CONTIGS_FASTA" \
    --accession "$ACCESSION" \
    --top-hits "${SAMPLE_NAME}_top_hits.tsv" \
    --viral-hits "${SAMPLE_NAME}_viral_blast.tsv" \
    --report "${SAMPLE_NAME}_blast_report_section.txt"

echo ""
echo "--- Top Hits TSV ---"
if [ -s "${SAMPLE_NAME}_top_hits.tsv" ]; then
    cat "${SAMPLE_NAME}_top_hits.tsv"
else
    echo "(empty)"
fi

echo ""
echo "--- Report Section ---"
if [ -s "${SAMPLE_NAME}_blast_report_section.txt" ]; then
    cat "${SAMPLE_NAME}_blast_report_section.txt"
else
    echo "(empty)"
fi

# =============================================================================
# Step 3: Generate HTML report (if qc_with_simple_plots.py available)
# =============================================================================
echo ""
echo "========================================="
echo "Step 3: HTML Report Generation"
echo "========================================="

QC_SCRIPT="${PIPELINE_DIR}/qc_with_simple_plots.py"

# The QC script expects a diagnostic_SAMPLE directory structure
# Check if we're already in one
CURRENT_DIR=$(basename "$WORK_DIR")
if [[ "$CURRENT_DIR" == diagnostic_* ]]; then
    # We're inside the diagnostic dir — run from parent
    DIAG_DIR_NAME="$CURRENT_DIR"
    PARENT_DIR=$(dirname "$WORK_DIR")

    if [ -f "$QC_SCRIPT" ]; then
        echo "Generating HTML report..."
        cd "$PARENT_DIR"
        python "$QC_SCRIPT" "$DIAG_DIR_NAME"
        QC_EXIT=$?
        cd "$WORK_DIR"
        if [ $QC_EXIT -eq 0 ]; then
            echo "HTML report generated successfully"
        else
            echo "WARNING: HTML report generation failed (exit code: $QC_EXIT)"
        fi
    else
        echo "QC script not found at: $QC_SCRIPT"
        echo "Skipping HTML report generation"
    fi
else
    echo "Not inside a diagnostic_* directory — skipping HTML report"
    echo "To generate HTML report, ensure contigs are in a diagnostic_SAMPLE/ directory"
fi

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "========================================="
echo "TEST COMPLETE"
echo "========================================="
echo ""
echo "Output files:"
echo "  ${WORK_DIR}/${SAMPLE_NAME}_blast_all.tsv      ($TOTAL_HITS total hits)"
echo "  ${WORK_DIR}/${SAMPLE_NAME}_top_hits.tsv        (top hits with ties)"
echo "  ${WORK_DIR}/${SAMPLE_NAME}_viral_blast.tsv     (viral hits only)"
echo "  ${WORK_DIR}/${SAMPLE_NAME}_blast_report_section.txt"
if [[ "$CURRENT_DIR" == diagnostic_* ]]; then
    echo "  ${WORK_DIR}/${CURRENT_DIR}_presentation_ready_report.html"
fi
echo ""
echo "To re-run with different parameters, just run this script again."
echo "No need to re-run mapping or assembly!"
