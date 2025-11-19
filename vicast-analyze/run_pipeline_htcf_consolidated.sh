#!/bin/bash
# Consolidated viral pipeline for HTCF with automatic snpEff database building

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads] [--large-files|--extremely-large-files]"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --large-files"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 4 --extremely-large-files"
    echo ""
    echo "Prerequisites:"
    echo "  Activate conda before running: source /ref/sahlab/software/anaconda3/bin/activate"
    echo ""
    echo "Configuration:"
    echo "  Create pipeline_config.sh from pipeline_config.template.sh to customize paths"
    echo ""
    echo "This script will:"
    echo "  1. Check if snpEff database exists for the accession"
    echo "  2. Build the database if needed (viral-friendly settings)"
    echo "  3. Run the complete viral genomics pipeline"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3
THREADS=${4:-4}  # Default to 4 threads if not specified
LARGE_FILES_FLAG=""

# Check if memory flags are provided (can be in position 4 or 5)
if [ $# -ge 4 ] && [ "$4" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
    THREADS=4  # Reset to default since flag was in threads position
elif [ $# -ge 4 ] && [ "$4" == "--extremely-large-files" ]; then
    LARGE_FILES_FLAG="--extremely-large-files"
    THREADS=4  # Reset to default since flag was in threads position
elif [ $# -ge 5 ] && [ "$5" == "--large-files" ]; then
    LARGE_FILES_FLAG="--large-files"
elif [ $# -ge 5 ] && [ "$5" == "--extremely-large-files" ]; then
    LARGE_FILES_FLAG="--extremely-large-files"
fi

# Source configuration file if it exists, with fallbacks
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"

if [ -f "$CONFIG_FILE" ]; then
    echo "Loading configuration from: $CONFIG_FILE"
    source "$CONFIG_FILE"
else
    echo "⚠️  No configuration file found. Using hardcoded paths."
    # Fallback defaults for backward compatibility
    MAMBA_CMD="conda run -n viral_genomics"
    SNPEFF_DIR="/ref/sahlab/software/snpEff"
    SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
    JAVA_PATH="java"
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

# Set up environment
echo "Setting up HTCF environment..."

# Set paths
PIPELINE_DIR="$(cd "$(dirname "$0")"; pwd)"

# Validate input files exist
if [ ! -f "$R1" ]; then
    echo "Error: R1 file not found: $R1"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: R2 file not found: $R2"
    exit 1
fi

echo "Running viral pipeline with:"
echo "  R1: $R1"
echo "  R2: $R2"
echo "  Accession: $ACCESSION"
echo "  Threads: $THREADS"
echo "  Memory mode: $(if [ "$LARGE_FILES_FLAG" == "--extremely-large-files" ]; then echo "extreme (256GB)"; elif [ "$LARGE_FILES_FLAG" == "--large-files" ]; then echo "large (64GB)"; else echo "standard (32GB)"; fi)"
echo "  Working directory: $(pwd)"

# Function to build snpEff database
build_snpeff_database() {
    local acc=$1
    echo "Building snpEff database for $acc..."
    
    # Create database directory
    mkdir -p "$SNPEFF_DIR/data/$acc"
    
    # Download GenBank file to variants folder (for collaborators)
    mkdir -p ./cleaned_seqs/variants
    if [ ! -f "./cleaned_seqs/variants/$acc.gb" ]; then
        echo "Downloading GenBank file for $acc..."
        $MAMBA_CMD efetch -db nucleotide -id "$acc" -format gb > "./cleaned_seqs/variants/$acc.gb"
        if [ $? -ne 0 ] || [ ! -s "./cleaned_seqs/variants/$acc.gb" ]; then
            echo "Error: Failed to download GenBank file for $acc"
            return 1
        fi
        echo "GenBank file downloaded successfully ($(wc -l < "./cleaned_seqs/variants/$acc.gb") lines)"
    fi
    
    # Copy GenBank file to snpEff database directory
    cp "./cleaned_seqs/variants/$acc.gb" "$SNPEFF_DIR/data/$acc/genes.gbk"
    
    # Download FASTA file if not exists
    if [ ! -f "$SNPEFF_DIR/data/$acc/sequences.fa" ]; then
        echo "Downloading FASTA file for $acc..."
        $MAMBA_CMD efetch -db nucleotide -id "$acc" -format fasta > "$SNPEFF_DIR/data/$acc/sequences.fa"
        if [ $? -ne 0 ] || [ ! -s "$SNPEFF_DIR/data/$acc/sequences.fa" ]; then
            echo "Error: Failed to download FASTA file for $acc"
            return 1
        fi
        echo "FASTA file downloaded successfully"
    fi
    
    # Add to snpEff config BEFORE building - proper syntax
    if ! grep -q "$acc.genome" "$SNPEFF_DIR/snpEff.config"; then
        echo "Adding $acc to snpEff configuration..."
        echo "" >> "$SNPEFF_DIR/snpEff.config"
        echo "# $acc" >> "$SNPEFF_DIR/snpEff.config"
        echo "$acc.genome : $acc" >> "$SNPEFF_DIR/snpEff.config"
    fi
    
    # Verify config was updated
    echo "Verifying config file contains $acc..."
    if grep -q "$acc.genome" "$SNPEFF_DIR/snpEff.config"; then
        echo "Config file updated successfully"
    else
        echo "Error: Failed to update config file"
        return 1
    fi
    
    # Build database with viral-friendly settings using mamba run
    echo "Building snpEff database with viral-friendly settings..."
    cd "$SNPEFF_DIR"
    $MAMBA_CMD java -jar snpEff.jar build -v -noCheckProtein -noCheckCds -noLog "$acc"
    local build_result=$?
    
    # Return to original directory
    cd - > /dev/null
    
    if [ $build_result -eq 0 ]; then
        echo "Successfully built snpEff database for $acc"
        return 0
    else
        echo "Error: Failed to build snpEff database for $acc"
        return 1
    fi
}

# Check if snpEff database exists using mamba run
echo "Checking snpEff database for $ACCESSION..."
if $MAMBA_CMD java -jar "$SNPEFF_JAR" databases 2>/dev/null | grep -q "$ACCESSION"; then
    echo "snpEff database found for $ACCESSION"
    
    # Check if the database is actually built (has snpEffectPredictor.bin)
    if [ ! -f "$SNPEFF_DIR/data/$ACCESSION/snpEffectPredictor.bin" ]; then
        echo "Database exists but not built. Building now..."
        build_snpeff_database "$ACCESSION"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to build snpEff database"
            exit 1
        fi
    fi
else
    echo "snpEff database not found for $ACCESSION. Building now..."
    build_snpeff_database "$ACCESSION"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to build snpEff database"
        exit 1
    fi
fi

echo "Running the viral pipeline..."

# Run the pipeline without --add-to-snpeff since we handled it above
if [ -n "$LARGE_FILES_FLAG" ]; then
    echo "Using $LARGE_FILES_FLAG for increased memory allocation"
    $MAMBA_CMD python ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH" \
        --large-files
else
    $MAMBA_CMD python ${PIPELINE_DIR}/viral_pipeline.py \
        --r1 "$R1" \
        --r2 "$R2" \
        --accession "$ACCESSION" \
        --threads $THREADS \
        --snpeff-jar "$SNPEFF_JAR" \
        --java-path "$JAVA_PATH"
fi

PIPELINE_EXIT_CODE=$?

if [ $PIPELINE_EXIT_CODE -eq 0 ]; then
    echo ""
    echo "Pipeline completed successfully!"
    echo ""
    
    # Extract sample name from R1 filename to only show files for current sample
    SAMPLE_NAME=$(basename "$R1" | sed 's/_R1\.fastq\.gz$//' | sed 's/_R1\.qc\.fastq\.gz$//')
    
    echo "Output files:"
    find ./cleaned_seqs/variants -name "*${SAMPLE_NAME}*.tsv" -o -name "*${SAMPLE_NAME}*.vcf" 2>/dev/null | sort
    echo ""
    echo "Key results:"
    if ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv >/dev/null 2>&1; then
        echo "  Annotation TSV: $(ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv | head -1)"
        echo "  Total annotated variants: $(wc -l < $(ls ./cleaned_seqs/variants/*${SAMPLE_NAME}*.snpEFF.ann.tsv | head -1))"
    fi
else
    echo "Pipeline failed with exit code $PIPELINE_EXIT_CODE"
    exit $PIPELINE_EXIT_CODE
fi