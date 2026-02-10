#!/bin/bash
# VICAST-Analyze: Post-Processing Workflow (Chunk 3)
# Runs: Parse mutations + Generate haplotype consensus
# REQUIRES: Annotation workflow must be run first (run_vicast_analyze_annotate_only.sh)

# Default parameter values
SAMPLE_NAME=""
ACCESSION=""
QUALITY=1000
DEPTH=200
FREQ_PARSE=0.01
FREQ_CONSENSUS=0.50

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --sample)
            SAMPLE_NAME="$2"
            shift 2
            ;;
        --accession)
            ACCESSION="$2"
            shift 2
            ;;
        --quality)
            QUALITY="$2"
            shift 2
            ;;
        --depth)
            DEPTH="$2"
            shift 2
            ;;
        --freq-parse)
            FREQ_PARSE="$2"
            shift 2
            ;;
        --freq-consensus)
            FREQ_CONSENSUS="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --sample <sample_name> --accession <accession> [options]"
            echo ""
            echo "Required arguments:"
            echo "  --sample <name>         Sample name (must match .snpEFF.ann.tsv filename)"
            echo "  --accession <accession> Virus accession number (e.g., AY532665.1)"
            echo ""
            echo "Optional arguments:"
            echo "  --quality <int>         Minimum quality score (default: 1000)"
            echo "  --depth <int>           Minimum coverage depth (default: 200)"
            echo "  --freq-parse <float>    Min allele frequency for parsing (default: 0.01)"
            echo "  --freq-consensus <float> Min allele frequency for consensus (default: 0.50)"
            echo ""
            echo "Examples:"
            echo "  # Use defaults (quality=1000, depth=200, freq_parse=0.01, freq_consensus=0.50)"
            echo "  $0 --sample WNV_sample --accession AY532665.1"
            echo ""
            echo "  # Custom quality and depth"
            echo "  $0 --sample WNV_sample --accession AY532665.1 --quality 500 --depth 100"
            echo ""
            echo "  # Detect rare quasispecies (0.5% frequency)"
            echo "  $0 --sample WNV_sample --accession AY532665.1 --freq-parse 0.005"
            echo ""
            echo "  # Very strict filtering"
            echo "  $0 --sample WNV_sample --accession AY532665.1 --quality 5000 --depth 500 --freq-parse 0.05 --freq-consensus 0.80"
            echo ""
            echo "  # Lower quality sample - relaxed thresholds"
            echo "  $0 --sample WNV_sample --accession AY532665.1 --quality 100 --depth 50 --freq-parse 0.02"
            echo ""
            echo "This script runs post-processing (Chunk 3):"
            echo "  1. Parse & filter mutations (parse_snpeff_tsv.py)"
            echo "  2. Generate realistic haplotypes & consensus (generate_realistic_haplotype_consensus.py)"
            echo ""
            echo "PREREQUISITES:"
            echo "  Run annotation workflow first: ./run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

# Check required arguments
if [ -z "$SAMPLE_NAME" ] || [ -z "$ACCESSION" ]; then
    echo "❌ Error: Missing required arguments"
    echo ""
    echo "Usage: $0 --sample <sample_name> --accession <accession> [options]"
    echo "Use --help for full usage information"
    exit 1
fi

# Source configuration file if it exists
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"

if [ -f "$CONFIG_FILE" ]; then
    source "$CONFIG_FILE"
else
    MAMBA_CMD="conda run -n viral_genomics"
    SNPEFF_DIR="${SNPEFF_DIR:-/ref/sahlab/software/snpEff}"  # Use env var or default
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found in PATH"
    echo "Please activate conda first:"
    echo "  source $CONDA_BASE/bin/activate (or your conda installation)"
    exit 1
fi

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate vicast_analyze

echo "============================================="
echo "VICAST-ANALYZE: POST-PROCESSING (Chunk 3)"
echo "============================================="
echo "  Sample: $SAMPLE_NAME"
echo "  Accession: $ACCESSION"
echo "  Quality threshold: $QUALITY"
echo "  Depth threshold: $DEPTH"
echo "  Frequency (parse): $FREQ_PARSE"
echo "  Frequency (consensus): $FREQ_CONSENSUS"
echo "============================================="
echo ""

# Verify prerequisite files exist
echo "Verifying prerequisite files..."
TSV_FILE="cleaned_seqs/variants/${SAMPLE_NAME}.snpEFF.ann.tsv"
REF_FILE="cleaned_seqs/${ACCESSION}.fasta"

if [ ! -f "$TSV_FILE" ]; then
    echo ""
    echo "❌ ERROR: Annotation TSV not found: $TSV_FILE"
    echo ""
    echo "Please run annotation workflow first:"
    echo "  ./run_vicast_analyze_annotate_only.sh <R1> <R2> $ACCESSION"
    exit 1
fi

if [ ! -f "$REF_FILE" ]; then
    echo ""
    echo "❌ ERROR: Reference FASTA not found: $REF_FILE"
    echo ""
    echo "Please run QC or annotation workflow first to download reference"
    exit 1
fi

echo "✓ Found annotation TSV: $TSV_FILE"
echo "✓ Found reference FASTA: $REF_FILE"
echo ""

# Create results directory if it doesn't exist
RESULTS_DIR="${SAMPLE_NAME}_results"
mkdir -p "$RESULTS_DIR"

# ============================================================================
# STEP 1: Parse & Filter Mutations
# ============================================================================
echo "============================================="
echo "STEP 1: Parse & Filter Mutations"
echo "============================================="
echo "Parameters:"
echo "  Quality: >= $QUALITY"
echo "  Depth: >= $DEPTH"
echo "  Allele Frequency: >= $FREQ_PARSE"
echo ""

FILTERED_TSV="${RESULTS_DIR}/${SAMPLE_NAME}_filtered_mutations.tsv"

python ${SCRIPT_DIR}/parse_snpeff_tsv.py \
  "$TSV_FILE" \
  "$FILTERED_TSV" \
  --quality $QUALITY \
  --depth $DEPTH \
  --freq $FREQ_PARSE

PARSE_EXIT_CODE=$?

if [ $PARSE_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "❌ Parse mutations failed with exit code $PARSE_EXIT_CODE"
    exit $PARSE_EXIT_CODE
fi

echo ""
echo "✓ Mutation parsing complete"
echo "  Output: $FILTERED_TSV"
echo ""

# ============================================================================
# STEP 2: Generate Realistic Haplotype Consensus
# ============================================================================
echo "============================================="
echo "STEP 2: Generate Haplotype Consensus"
echo "============================================="
echo "Parameters:"
echo "  Quality: >= $QUALITY"
echo "  Depth: >= $DEPTH"
echo "  Allele Frequency: >= $FREQ_CONSENSUS (majority rule)"
echo ""

# Check if known_viruses.json needs updating
echo "Checking virus configuration..."
python ${SCRIPT_DIR}/auto_add_virus_to_json.py \
  "$ACCESSION" \
  "$SNPEFF_DIR" \
  "${SCRIPT_DIR}/known_viruses.json"

AUTO_ADD_EXIT_CODE=$?

if [ $AUTO_ADD_EXIT_CODE -eq 2 ]; then
    echo "⚠️  Warning: Failed to auto-add virus configuration"
    echo "   Continuing anyway - will generate polyprotein only"
elif [ $AUTO_ADD_EXIT_CODE -eq 0 ]; then
    echo "✓ Virus configuration ready"
fi

echo ""

OUTPUT_PREFIX="${RESULTS_DIR}/${SAMPLE_NAME}_consensus"

python ${SCRIPT_DIR}/generate_realistic_haplotype_consensus.py \
  --vcf "$FILTERED_TSV" \
  --reference "$REF_FILE" \
  --accession "$ACCESSION" \
  --quality $QUALITY \
  --depth $DEPTH \
  --freq $FREQ_CONSENSUS \
  --output-prefix "$OUTPUT_PREFIX"

CONSENSUS_EXIT_CODE=$?

if [ $CONSENSUS_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "❌ Consensus generation failed with exit code $CONSENSUS_EXIT_CODE"
    exit $CONSENSUS_EXIT_CODE
fi

echo ""
echo "============================================="
echo "✓ POST-PROCESSING COMPLETE (Chunk 3)"
echo "============================================="
echo ""
echo "OUTPUT FILES:"
echo "  1. Filtered mutations:"
echo "     ${FILTERED_TSV}"
echo ""
echo "  2. Consensus genome:"
echo "     ${OUTPUT_PREFIX}_consensus.fasta"
echo ""
echo "  3. Haplotype proteins:"
echo "     ${OUTPUT_PREFIX}_proteins.fasta"
echo ""
echo "  4. Summary report:"
echo "     ${OUTPUT_PREFIX}_summary_report.txt"
echo ""
echo "NEXT STEPS:"
echo "  Review summary report:"
echo "    cat ${OUTPUT_PREFIX}_summary_report.txt"
echo ""
echo "  Check filtered mutations:"
echo "    head -20 ${FILTERED_TSV}"
echo ""
echo "  View haplotype proteins:"
echo "    grep '>' ${OUTPUT_PREFIX}_proteins.fasta"
echo ""
echo "============================================="
