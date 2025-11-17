#!/bin/bash
# Enhanced Module 1: Viral Pipeline with Automatic Virus Configuration
# This version automatically adds new viruses to known_viruses.json during assembly

# Function to display usage
usage() {
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> <threads> [--large-files]"
    echo ""
    echo "Arguments:"
    echo "  R1_fastq      : Path to R1 FASTQ file"
    echo "  R2_fastq      : Path to R2 FASTQ file"
    echo "  accession     : GenBank accession number (e.g., NC_001477.1)"
    echo "  threads       : Number of threads to use"
    echo "  --large-files : Optional flag for large file handling"
    echo ""
    echo "Prerequisites:"
    echo "  Activate conda before running: source /ref/sahlab/software/anaconda3/bin/activate"
    echo ""
    echo "Configuration:"
    echo "  Create pipeline_config.sh from pipeline_config.template.sh to customize paths"
    echo ""
    echo "This enhanced version automatically:"
    echo "  - Downloads virus reference and annotation"
    echo "  - Adds new viruses to configuration database"
    echo "  - Verifies proper gene annotation exists"
    echo "  - Runs complete assembly and variant calling"
    exit 1
}

# Check if correct number of arguments provided
if [ $# -lt 4 ]; then
    usage
fi

# Assign arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=$4
LARGE_FILES_FLAG=""

# Check for optional flags
if [ $# -eq 5 ] && [ "$5" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
fi

# Extract sample name from R1 filename
SAMPLE_NAME=$(basename "$R1" | sed 's/_R[12]\..*//')

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Source configuration file if it exists
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"
if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "⚠️  Configuration file not found: $CONFIG_FILE"
    echo "   Using hardcoded paths (consider creating pipeline_config.sh from template)"
    # Fallback to hardcoded paths for backward compatibility
    MAMBA_CMD="conda run -n viral_genomics"
    PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
fi

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found in PATH"
    echo ""
    echo "Please activate conda first:"
    echo "  source /ref/sahlab/software/anaconda3/bin/activate"
    echo ""
    echo "Then re-run this script."
    exit 1
fi

# Set up HTCF environment
echo "========================================="
echo "Enhanced Viral Pipeline - Module 1"
echo "========================================="
echo "Sample: $SAMPLE_NAME"
echo "Accession: $ACCESSION"
echo "Working directory: $(pwd)"
echo "Time: $(date)"
echo "========================================="

# Step 1: Check virus configuration before running pipeline
echo ""
echo "Step 1: Checking virus configuration..."
KNOWN_VIRUSES_FILE="${PIPELINE_BASE}/viral_pipeline/visualization/known_viruses.json"

VIRUS_CONFIGURED=false
if [ -f "$KNOWN_VIRUSES_FILE" ]; then
    if grep -q "\"$ACCESSION\"" "$KNOWN_VIRUSES_FILE"; then
        echo "✅ Virus $ACCESSION found in configuration"
        VIRUS_CONFIGURED=true
    else
        echo "⚠️  Virus $ACCESSION not in configuration database"
    fi
else
    echo "⚠️  Configuration file not found, will create new one"
    mkdir -p "$(dirname "$KNOWN_VIRUSES_FILE")"
    echo "{}" > "$KNOWN_VIRUSES_FILE"
fi

# Step 2: If virus not configured, download GenBank and add it
if [ "$VIRUS_CONFIGURED" = false ]; then
    echo ""
    echo "Step 2: Downloading virus information..."
    
    # Create temporary directory for downloads
    TEMP_DIR="./temp_virus_setup_$$"
    mkdir -p "$TEMP_DIR"
    
    # Download GenBank file
    echo "Downloading GenBank file for $ACCESSION..."
    $MAMBA_CMD efetch -db nucleotide -id "$ACCESSION" -format gb > "$TEMP_DIR/${ACCESSION}.gb"
    
    if [ $? -ne 0 ] || [ ! -s "$TEMP_DIR/${ACCESSION}.gb" ]; then
        echo "❌ Error: Failed to download GenBank file for $ACCESSION"
        rm -rf "$TEMP_DIR"
        exit 1
    fi
    
    echo "GenBank file downloaded successfully"
    
    # Add virus to configuration
    echo "Adding virus to configuration database..."
    $MAMBA_CMD python3 "${PIPELINE_BASE}/viral_pipeline/utils/add_virus_to_config.py" \
        "$TEMP_DIR/${ACCESSION}.gb" \
        --config "$KNOWN_VIRUSES_FILE"
    
    if [ $? -eq 0 ]; then
        echo "✅ Successfully added $ACCESSION to virus configuration!"
        
        # Verify it has proper gene annotation (not just polyprotein)
        GENE_COUNT=$($MAMBA_CMD python3 -c "
import json
with open('$KNOWN_VIRUSES_FILE', 'r') as f:
    data = json.load(f)
    if '$ACCESSION' in data:
        print(len(data['$ACCESSION'].get('gene_coords', {})))
    else:
        print(0)
")
        
        if [ "$GENE_COUNT" -gt 1 ]; then
            echo "✅ Found $GENE_COUNT genes in annotation"
        else
            echo "⚠️  Warning: Only found $GENE_COUNT gene(s). Manual configuration may be needed."
            echo "   The GenBank file may use polyprotein annotation."
            echo "   Consensus generation will use polyprotein fallback."
        fi
    else
        echo "⚠️  Failed to add virus to configuration. Continuing anyway..."
    fi
    
    # Clean up
    rm -rf "$TEMP_DIR"
else
    echo "Step 2: Skipped (virus already configured)"
fi

# Step 3: Run the main pipeline
echo ""
echo "Step 3: Running viral assembly pipeline..."
echo "========================================="

# Run the consolidated pipeline
bash "${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh" "$R1" "$R2" "$ACCESSION" "$THREADS" $LARGE_FILES_FLAG

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "========================================="
    echo "✅ Pipeline completed successfully!"
    echo "========================================="
    
    # Copy known_viruses.json to results directory for reference
    if [ -f "$KNOWN_VIRUSES_FILE" ]; then
        cp "$KNOWN_VIRUSES_FILE" "./cleaned_seqs/variants/"
        echo "Virus configuration copied to: ./cleaned_seqs/variants/known_viruses.json"
    fi
    
    # Create next steps file
    cat > "next_steps_${SAMPLE_NAME}.txt" << EOF
# Next steps for ${SAMPLE_NAME} (${ACCESSION})

# Module 2: Generate depth file
${MAMBA_CMD} bash -c '
  mkdir -p ${SAMPLE_NAME}_results &&
  echo -e "chrom\\tposition\\tdepth" > ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt &&
  samtools depth cleaned_seqs/mapping/${SAMPLE_NAME}.lofreq.final.bam >> ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt'

# Module 3: Parse mutations
${MAMBA_CMD} \\
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/parse_snpeff_tsv.py \\
  "cleaned_seqs/variants/${SAMPLE_NAME}.snpEFF.ann.tsv" \\
  "${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv" \\
  --quality 1000 --depth 200 --freq 0.01

# Modules 4 & 5 (visualization) removed - focus on data generation
# Visualization broke with new genomes due to hardcoded assumptions
# Use generated TSV files for custom visualization as needed

# Module 7: Diagnostic Report (Optional)
sbatch ${PIPELINE_BASE}/viral_pipeline/analysis/submit_viral_diagnostic.sh \\
  "${R1}" \\
  "${R2}" \\
  ${ACCESSION} \\
  diagnostic_${SAMPLE_NAME} \\
  4

# Module 8: Generate consensus
${MAMBA_CMD} \\
  python3 ${PIPELINE_BASE}/viral_pipeline/utils/generate_filtered_consensus.py \\
  --vcf ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \\
  --reference cleaned_seqs/${ACCESSION}.fasta \\
  --accession ${ACCESSION} \\
  --quality 1000 --depth 200 --freq 0.50 \\
  --output-prefix ${SAMPLE_NAME}_results/${SAMPLE_NAME}_consensus
EOF
    
    echo ""
    echo "📝 Visualization commands saved to: next_steps_${SAMPLE_NAME}.txt"
    echo ""
else
    echo ""
    echo "❌ Pipeline failed with exit code: $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi