#!/bin/bash
# VICAST-Analyze: Post-Processing Workflow
# Two-tier pipeline:
#   Step 1: Parse & filter mutations (parse_snpeff_tsv.py)
#   Step 2: Generate consensus genome + proteins (generate_consensus_genome.py)
#   Step 3: BAM co-occurrence analysis (check_read_cooccurrence.py)
#   Step 4: Generate variant genomes + proteins (generate_variant_genomes.py)
#
# REQUIRES: Annotation workflow must be run first (run_vicast_analyze_annotate_only.sh)

# Default parameter values
SAMPLE_NAME=""
ACCESSION=""
BAM_FILE=""
QUALITY=1000
DEPTH=200
FREQ_PARSE=0.01
CONSENSUS_AF=""
GENOME_TYPE=""
MIN_MINOR_AF=0.005
MAX_MINOR_AF=0.05

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
        --bam)
            BAM_FILE="$2"
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
        --consensus-af)
            CONSENSUS_AF="$2"
            shift 2
            ;;
        --genome-type)
            GENOME_TYPE="$2"
            shift 2
            ;;
        --min-minor-af)
            MIN_MINOR_AF="$2"
            shift 2
            ;;
        --max-minor-af)
            MAX_MINOR_AF="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --sample <name> --accession <accession> --bam <bam_file> [options]"
            echo ""
            echo "Required arguments:"
            echo "  --sample <name>           Sample name (must match .snpEFF.ann.tsv filename)"
            echo "  --accession <accession>   Virus accession number (e.g., AY532665.1)"
            echo "  --bam <file>              BAM file from QC step (sorted, indexed)"
            echo ""
            echo "Optional arguments:"
            echo "  --quality <int>           Minimum quality score (default: 1000)"
            echo "  --depth <int>             Minimum coverage depth (default: 200)"
            echo "  --freq-parse <float>      Min allele frequency for parsing (default: 0.01)"
            echo "  --consensus-af <float>    Override consensus AF threshold"
            echo "                            (default: auto — 0.95 for ssRNA/ssDNA, 0.45 for dsRNA)"
            echo "  --genome-type <str>       Override genome type: ss|ds|ssRNA|dsRNA"
            echo "                            (default: auto from known_viruses.json)"
            echo "  --min-minor-af <float>    Min AF for minor variant genomes (default: 0.005)"
            echo "  --max-minor-af <float>    Max AF for minor variant genomes (default: 0.05)"
            echo ""
            echo "Pipeline steps:"
            echo "  1. Parse & filter mutations        (parse_snpeff_tsv.py)"
            echo "  2. Generate consensus genome       (generate_consensus_genome.py)"
            echo "  3. BAM co-occurrence analysis      (check_read_cooccurrence.py)"
            echo "  4. Generate minor variant genomes  (generate_variant_genomes.py)"
            echo ""
            echo "Examples:"
            echo "  # Standard ssRNA virus (AF >= 0.95 for consensus)"
            echo "  $0 --sample DRR878516 --accession NC_045512.2 \\"
            echo "     --bam cleaned_seqs/variants/DRR878516_aligned_sorted.bam"
            echo ""
            echo "  # Override consensus threshold for dsRNA"
            echo "  $0 --sample MY_SAMPLE --accession NC_XXXXXX.X \\"
            echo "     --bam my_bam.bam --genome-type ds"
            echo ""
            echo "  # Custom minor variant range (1%-10%)"
            echo "  $0 --sample MY_SAMPLE --accession NC_XXXXXX.X \\"
            echo "     --bam my_bam.bam --min-minor-af 0.01 --max-minor-af 0.10"
            echo ""
            echo "PREREQUISITES:"
            echo "  Run annotation workflow first: run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>"
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
if [ -z "$SAMPLE_NAME" ] || [ -z "$ACCESSION" ] || [ -z "$BAM_FILE" ]; then
    echo "Error: Missing required arguments"
    echo ""
    echo "Usage: $0 --sample <name> --accession <accession> --bam <bam_file> [options]"
    echo "Use --help for full usage information"
    exit 1
fi

# =============================================================================
# Configuration Loading
# =============================================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
VICAST_HOME="$( cd "${SCRIPT_DIR}/.." && pwd )"

# Source user configuration if it exists
if [ -f "$HOME/.vicast/config.sh" ]; then
    source "$HOME/.vicast/config.sh"
elif [ -f "${VICAST_HOME}/vicast_config.template.sh" ]; then
    source "${VICAST_HOME}/vicast_config.template.sh"
fi

# Source local pipeline config if it exists
if [ -f "${SCRIPT_DIR}/pipeline_config.sh" ]; then
    source "${SCRIPT_DIR}/pipeline_config.sh"
fi

# =============================================================================
# Environment Detection
# =============================================================================

# Detect SnpEff paths (needed for auto_add_virus_to_json.py)
if [ -z "$SNPEFF_DIR" ]; then
    if [ -n "$SNPEFF_HOME" ]; then
        SNPEFF_DIR="${SNPEFF_HOME}"
    fi
fi

# Detect conda environment
VICAST_CONDA_ENV="${VICAST_CONDA_ENV:-vicast_analyze}"

# Check if conda/micromamba is available (optional in Docker)
CONDA_AVAILABLE=false
if command -v conda &> /dev/null; then
    CONDA_AVAILABLE=true
    CONDA_CMD="conda"
elif command -v micromamba &> /dev/null; then
    CONDA_AVAILABLE=true
    CONDA_CMD="micromamba"
    VICAST_CONDA_ENV="base"
fi

# Activate conda/micromamba environment if available
if [ "$CONDA_AVAILABLE" = true ]; then
    if [ "$CONDA_CMD" = "conda" ]; then
        eval "$(conda shell.bash hook)"
        if conda activate "$VICAST_CONDA_ENV" 2>/dev/null; then
            echo "Activated conda environment: $VICAST_CONDA_ENV"
        else
            echo "Warning: Could not activate conda environment '$VICAST_CONDA_ENV'"
        fi
    elif [ "$CONDA_CMD" = "micromamba" ]; then
        echo "Detected micromamba (Docker environment) - using current environment"
    fi
else
    echo "Warning: conda/micromamba not found - assuming tools are in PATH"
fi

echo "============================================="
echo "VICAST-ANALYZE: POST-PROCESSING"
echo "============================================="
echo "  Sample:              $SAMPLE_NAME"
echo "  Accession:           $ACCESSION"
echo "  BAM file:            $BAM_FILE"
echo "  Quality threshold:   $QUALITY"
echo "  Depth threshold:     $DEPTH"
echo "  Frequency (parse):   $FREQ_PARSE"
if [ -n "$CONSENSUS_AF" ]; then
    echo "  Consensus AF:        $CONSENSUS_AF (user override)"
else
    echo "  Consensus AF:        auto (genome-type aware)"
fi
if [ -n "$GENOME_TYPE" ]; then
    echo "  Genome type:         $GENOME_TYPE (user override)"
else
    echo "  Genome type:         auto (from known_viruses.json)"
fi
echo "  Minor variant range: $MIN_MINOR_AF - $MAX_MINOR_AF"
echo "============================================="
echo ""

# Verify prerequisite files exist
echo "Verifying prerequisite files..."
TSV_FILE="cleaned_seqs/variants/${SAMPLE_NAME}.snpEFF.ann.tsv"
REF_FILE="cleaned_seqs/${ACCESSION}.fasta"

if [ ! -f "$TSV_FILE" ]; then
    echo ""
    echo "ERROR: Annotation TSV not found: $TSV_FILE"
    echo ""
    echo "Please run annotation workflow first:"
    echo "  run_vicast_analyze_annotate_only.sh <R1> <R2> $ACCESSION"
    exit 1
fi

if [ ! -f "$REF_FILE" ]; then
    echo ""
    echo "ERROR: Reference FASTA not found: $REF_FILE"
    echo ""
    echo "Please run QC or annotation workflow first to download reference"
    exit 1
fi

if [ ! -f "$BAM_FILE" ]; then
    echo ""
    echo "ERROR: BAM file not found: $BAM_FILE"
    exit 1
fi

echo "  Found annotation TSV: $TSV_FILE"
echo "  Found reference FASTA: $REF_FILE"
echo "  Found BAM file: $BAM_FILE"
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
    echo "Parse mutations failed with exit code $PARSE_EXIT_CODE"
    exit $PARSE_EXIT_CODE
fi

echo ""
echo "Mutation parsing complete"
echo "  Output: $FILTERED_TSV"
echo ""

# ============================================================================
# STEP 2: Generate Consensus Genome (Tier 1)
# ============================================================================
echo "============================================="
echo "STEP 2: Generate Consensus Genome (Tier 1)"
echo "============================================="

# Check if known_viruses.json needs updating
echo "Checking virus configuration..."
python ${SCRIPT_DIR}/auto_add_virus_to_json.py \
  "$ACCESSION" \
  "$SNPEFF_DIR" \
  "${SCRIPT_DIR}/known_viruses.json"

AUTO_ADD_EXIT_CODE=$?

if [ $AUTO_ADD_EXIT_CODE -eq 2 ]; then
    echo "Warning: Failed to auto-add virus configuration"
    echo "  Continuing anyway - will generate polyprotein only"
elif [ $AUTO_ADD_EXIT_CODE -eq 0 ]; then
    echo "Virus configuration ready"
fi
echo ""

OUTPUT_PREFIX="${RESULTS_DIR}/${SAMPLE_NAME}"

# Build consensus-genome args
CONSENSUS_ARGS="--vcf $FILTERED_TSV --reference $REF_FILE --accession $ACCESSION --output-prefix $OUTPUT_PREFIX"
if [ -n "$CONSENSUS_AF" ]; then
    CONSENSUS_ARGS="$CONSENSUS_ARGS --consensus-af $CONSENSUS_AF"
fi
if [ -n "$GENOME_TYPE" ]; then
    CONSENSUS_ARGS="$CONSENSUS_ARGS --genome-type $GENOME_TYPE"
fi

python ${SCRIPT_DIR}/generate_consensus_genome.py $CONSENSUS_ARGS

CONSENSUS_EXIT_CODE=$?

if [ $CONSENSUS_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "Consensus generation failed with exit code $CONSENSUS_EXIT_CODE"
    exit $CONSENSUS_EXIT_CODE
fi

echo ""

# ============================================================================
# STEP 3: BAM Co-occurrence Analysis
# ============================================================================
echo "============================================="
echo "STEP 3: BAM Co-occurrence Analysis"
echo "============================================="

COOCCURRENCE_TSV="${RESULTS_DIR}/${SAMPLE_NAME}_cooccurrence.tsv"

python ${SCRIPT_DIR}/check_read_cooccurrence.py \
  --bam "$BAM_FILE" \
  --vcf "$FILTERED_TSV" \
  --output "$COOCCURRENCE_TSV" \
  --min-freq "$MIN_MINOR_AF"

COOCCURRENCE_EXIT_CODE=$?

if [ $COOCCURRENCE_EXIT_CODE -ne 0 ]; then
    echo ""
    echo "Warning: Co-occurrence analysis exited with code $COOCCURRENCE_EXIT_CODE"
    echo "  This is non-fatal — variant genomes will treat all variants as unlinked"
    # Create empty co-occurrence file so Step 4 can proceed
    echo "variant1_pos	variant2_pos	distance	evidence_level	total_spanning_reads	informative_reads	both_alt	variant1_only	variant2_only	both_ref	undetermined	cooccurrence_rate	linkage_proven" > "$COOCCURRENCE_TSV"
fi

echo ""

# ============================================================================
# STEP 4: Generate Variant Genomes (Tier 2)
# ============================================================================
echo "============================================="
echo "STEP 4: Generate Variant Genomes (Tier 2)"
echo "============================================="

CONSENSUS_FASTA="${OUTPUT_PREFIX}_consensus.fasta"

if [ ! -f "$CONSENSUS_FASTA" ]; then
    echo "Warning: Consensus FASTA not found: $CONSENSUS_FASTA"
    echo "  Skipping variant genome generation"
else
    python ${SCRIPT_DIR}/generate_variant_genomes.py \
      --vcf "$FILTERED_TSV" \
      --cooccurrence "$COOCCURRENCE_TSV" \
      --consensus "$CONSENSUS_FASTA" \
      --accession "$ACCESSION" \
      --min-af "$MIN_MINOR_AF" \
      --max-af "$MAX_MINOR_AF" \
      --output-prefix "$OUTPUT_PREFIX"

    VARIANT_EXIT_CODE=$?

    if [ $VARIANT_EXIT_CODE -ne 0 ]; then
        echo ""
        echo "Warning: Variant genome generation exited with code $VARIANT_EXIT_CODE"
        echo "  Consensus genome was generated successfully"
    fi
fi

echo ""

# ============================================================================
# SUMMARY
# ============================================================================
echo "============================================="
echo "POST-PROCESSING COMPLETE"
echo "============================================="
echo ""
echo "OUTPUT FILES:"
echo ""
echo "  Filtered mutations:"
echo "    ${FILTERED_TSV}"
echo ""
echo "  Tier 1 - Consensus genome:"
echo "    ${OUTPUT_PREFIX}_consensus.fasta"
echo "    ${OUTPUT_PREFIX}_consensus_proteins.fasta"
echo "    ${OUTPUT_PREFIX}_consensus_report.txt"
echo ""
echo "  Co-occurrence analysis:"
echo "    ${COOCCURRENCE_TSV}"
echo ""
echo "  Tier 2 - Variant genomes:"
echo "    ${OUTPUT_PREFIX}_variant_genomes.fasta"
echo "    ${OUTPUT_PREFIX}_variant_proteins.fasta"
echo "    ${OUTPUT_PREFIX}_variant_report.txt"
echo ""
echo "NEXT STEPS:"
echo "  Review consensus report:"
echo "    cat ${OUTPUT_PREFIX}_consensus_report.txt"
echo ""
echo "  Review variant report:"
echo "    cat ${OUTPUT_PREFIX}_variant_report.txt"
echo ""
echo "  Check co-occurrence linkage:"
echo "    cat ${COOCCURRENCE_TSV}"
echo ""
echo "============================================="
