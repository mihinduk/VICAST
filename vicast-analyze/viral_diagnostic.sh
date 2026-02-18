#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --job-name=viral_diagnostic
#SBATCH --output=%x_%j.out

# Suppress micromamba shell initialization warnings (Docker environment)
export MAMBA_NO_BANNER=1
export MAMBA_ROOT_PREFIX=/opt/conda 2>/dev/null || true

# Viral Contamination Diagnostic Module
# Runs mapping check, assembly, and viral BLAST for sample quality assessment
# Usage: ./viral_diagnostic.sh <R1> <R2> <accession> [threads]
#    or: ./viral_diagnostic.sh <R1> <R2> <accession> <sample_name> [threads]
# Sample name is auto-derived from R1 filename if not provided.

# Check arguments
if [ $# -lt 3 ]; then
    echo "Usage: $0 <R1_fastq> <R2_fastq> <accession> [threads]"
    echo "   or: $0 <R1_fastq> <R2_fastq> <accession> <sample_name> [threads]"
    echo ""
    echo "Example: $0 SRR5992153_1.fastq.gz SRR5992153_2.fastq.gz NC_001474.2 8"
    echo "Example: $0 sample_R1.fastq.gz sample_R2.fastq.gz GQ433359.1 SMS_14 4"
    echo ""
    echo "If sample_name is omitted, it is derived from the R1 filename."
    echo ""
    echo "This diagnostic script will:"
    echo "  1. Check mapping statistics to reference"
    echo "  2. Perform de novo assembly with MEGAHIT"
    echo "  3. BLAST contigs against viral nt database"
    echo "  4. Generate contamination assessment report"
    exit 1
fi

# Parse arguments
R1=$1
R2=$2
ACCESSION=$3

# Smart argument parsing: detect whether $4 is sample_name or threads
# If $4 is purely numeric and there is no $5, treat $4 as threads
if [ $# -eq 3 ]; then
    # Only 3 args: auto-derive sample name, default threads
    SAMPLE_NAME=""
    THREADS=4
elif [ $# -eq 4 ] && [[ "$4" =~ ^[0-9]+$ ]]; then
    # 4 args with numeric $4: treat as threads, auto-derive sample name
    SAMPLE_NAME=""
    THREADS=$4
elif [ $# -eq 4 ]; then
    # 4 args with non-numeric $4: treat as sample name
    SAMPLE_NAME=$4
    THREADS=4
else
    # 5+ args: explicit sample name and threads
    SAMPLE_NAME=$4
    THREADS=${5:-4}
fi

# Auto-derive sample name from R1 if not provided
if [ -z "$SAMPLE_NAME" ]; then
    R1_BASE=$(basename "$R1")
    # Remove _R1/_R2/_1/_2 suffixes, _001, and .fastq.gz extension
    # Matches viral_pipeline.py logic: re.sub(r'(_R?[12])?(_001)?\.fastq\.gz$', '', r1_base)
    SAMPLE_NAME=$(echo "$R1_BASE" | sed -E 's/(_R?[12])?(_001)?\.fastq\.gz$//')
fi

EXTREME_MEMORY_FLAG=""

# Check for extreme memory flag (can appear in any position)
if [[ "$*" == *"--extremely-large-files"* ]]; then
    EXTREME_MEMORY_FLAG="--extremely-large-files"
    echo "EXTREME MEMORY MODE: Using high memory settings for large files"
fi

# Resolve pipeline directory (where this script and helper scripts live)
# Must handle SLURM path resolution (SLURM copies scripts to temp locations)
if [ -n "${SLURM_JOB_ID:-}" ]; then
    ORIGINAL_SCRIPT=$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K[^ ]+' || echo "")
    if [ -n "$ORIGINAL_SCRIPT" ]; then
        PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT")"
    else
        PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi
else
    PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Set up environment and paths
echo "========================================="
echo "VIRAL CONTAMINATION DIAGNOSTIC MODULE"
echo "========================================="
echo "Sample: $SAMPLE_NAME"
echo "R1: $R1"
echo "R2: $R2"
echo "Reference: $ACCESSION"
echo "Threads: $THREADS"
echo "Started at: $(date)"
echo "========================================="

# Create output directory (remove if exists to avoid conflicts)
DIAGNOSTIC_DIR="./diagnostic_${SAMPLE_NAME}"
if [ -d "$DIAGNOSTIC_DIR" ]; then
    echo "Removing existing diagnostic directory..."
    rm -rf "$DIAGNOSTIC_DIR"
fi
mkdir -p "$DIAGNOSTIC_DIR"
cd "$DIAGNOSTIC_DIR"

# Set up conda/micromamba environment commands
echo "Setting up environment commands..."

# Detect conda/micromamba (Docker vs native)
# IMPORTANT: Check micromamba FIRST because it provides a conda symlink
CONDA_AVAILABLE=false
CONDA_CMD=""
VICAST_ENV="vicast_analyze"

if command -v micromamba &> /dev/null; then
    CONDA_AVAILABLE=true
    CONDA_CMD="micromamba"
    VICAST_ENV="base"
    echo "Detected micromamba (Docker environment)"
elif command -v conda &> /dev/null; then
    CONDA_AVAILABLE=true
    CONDA_CMD="conda"
    echo "Detected conda"
fi

# Activate conda environment if available
if [ "$CONDA_AVAILABLE" = true ]; then
    if [ "$CONDA_CMD" = "micromamba" ]; then
        # Docker: tools already in base environment, no activation needed
        echo "Using micromamba base environment (tools already available)"
        GENOMICS_CMD=""  # Direct execution
    elif [ "$CONDA_CMD" = "conda" ]; then
        # Native conda: try to activate environment
        eval "$(conda shell.bash hook)"
        if conda activate "$VICAST_ENV" 2>/dev/null; then
            echo "Activated conda environment: $VICAST_ENV"
            GENOMICS_CMD=""  # Direct execution (already activated)
        else
            echo "Could not activate conda environment, using conda run"
            GENOMICS_CMD="conda run -n $VICAST_ENV"
        fi
    fi
else
    echo "Warning: conda/micromamba not found - assuming tools are in PATH"
    GENOMICS_CMD=""  # Direct execution
fi

# Initialize base filenames for consistent use throughout script
# Extract sample name (strip _1/_2 or _R1/_R2 suffix) - same pattern as wrapper scripts
SAMPLE_BASE=$(basename "$R1" | sed -E "s/_R?[12](_[0-9]+)?..*$//")
# Construct base filenames with _R1/_R2 suffix to match QC output files
R1_BASE="${SAMPLE_BASE}_R1"
R2_BASE="${SAMPLE_BASE}_R2"
echo "Base filenames: $R1_BASE and $R2_BASE"

# Step 1: Quick mapping check
echo ""
echo "Step 1: Mapping Statistics Check"
echo "================================"

# Get reference genome (reuse from main pipeline if available)
if [ ! -f "${ACCESSION}.fasta" ]; then
    # Check if reference exists in main pipeline output
    if [ -f "../cleaned_seqs/${ACCESSION}.fasta" ]; then
        echo "Reusing reference genome from main pipeline: ../cleaned_seqs/${ACCESSION}.fasta"
        ln -s "../cleaned_seqs/${ACCESSION}.fasta" "${ACCESSION}.fasta"
    elif [ -f "../${ACCESSION}.fasta" ]; then
        echo "Reusing reference genome from parent directory: ../${ACCESSION}.fasta"
        ln -s "../${ACCESSION}.fasta" "${ACCESSION}.fasta"
    else
        echo "Downloading reference genome: $ACCESSION"
        # Use Python BioPython (more reliable than efetch command)
        python -c "
from Bio import Entrez
import os
Entrez.email = os.environ.get('NCBI_EMAIL', 'vicast_docker@example.com')
try:
    handle = Entrez.efetch(db='nucleotide', id='${ACCESSION}', rettype='fasta', retmode='text')
    with open('${ACCESSION}.fasta', 'w') as f:
        f.write(handle.read())
    handle.close()
    print('Reference genome downloaded successfully')
except Exception as e:
    print(f'ERROR: Failed to download reference: {e}')
    exit(1)
"
        if [ ! -s "${ACCESSION}.fasta" ]; then
            echo "ERROR: Failed to download reference genome"
            exit 1
        fi
    fi
fi

# Index reference
echo "Indexing reference genome..."
if [ -z "$GENOMICS_CMD" ]; then
    bwa index "${ACCESSION}.fasta" 2>&1
    INDEX_EXIT=$?
else
    $GENOMICS_CMD bwa index "${ACCESSION}.fasta"
    INDEX_EXIT=$?
fi

if [ $INDEX_EXIT -ne 0 ]; then
    echo "ERROR: BWA indexing failed with exit code $INDEX_EXIT"
    echo "Reference: ${ACCESSION}.fasta"
    ls -lh "${ACCESSION}.fasta" 2>&1
    exit 1
fi

# Verify index files were created
if [ ! -f "${ACCESSION}.fasta.bwt" ]; then
    echo "ERROR: BWA index files not created"
    echo "Expected: ${ACCESSION}.fasta.{amb,ann,bwt,pac,sa}"
    ls -lh "${ACCESSION}.fasta"* 2>&1 || echo "No index files found"
    exit 1
fi
echo "✓ BWA index complete"

# Quick alignment for mapping stats using cleaned reads
echo "Performing quick alignment for mapping statistics using cleaned reads..."
# Look for cleaned reads in the parent cleaned_seqs directory
CLEANED_R1="../cleaned_seqs/${R1_BASE}.qc.fastq.gz"
CLEANED_R2="../cleaned_seqs/${R2_BASE}.qc.fastq.gz"

if [ -f "$CLEANED_R1" ] && [ -f "$CLEANED_R2" ]; then
    echo "  Using cleaned reads: $CLEANED_R1 and $CLEANED_R2"
    echo "  Reference: ${ACCESSION}.fasta"
    echo "  Output BAM: ${SAMPLE_NAME}_quick.bam"

    if [ -z "$GENOMICS_CMD" ]; then
        # Direct execution (Docker/micromamba base environment)
        echo "  Running BWA mem + samtools (direct execution)..."
        echo "  This may take several minutes for large datasets..."

        # Log files for debugging
        BWA_LOG="${SAMPLE_NAME}_bwa.log"
        SAMTOOLS_LOG="${SAMPLE_NAME}_samtools.log"

        set -o pipefail  # Fail if any command in pipe fails
        bwa mem -t $THREADS "${ACCESSION}.fasta" "$CLEANED_R1" "$CLEANED_R2" 2>"$BWA_LOG" | \
            samtools view -@ $THREADS -bS - 2>>"$SAMTOOLS_LOG" | \
            samtools sort -@ $THREADS -o "${SAMPLE_NAME}_quick.bam" - 2>>"$SAMTOOLS_LOG"
        BWA_EXIT=$?
        set +o pipefail

        if [ $BWA_EXIT -ne 0 ]; then
            echo "  ERROR: Mapping pipeline failed with exit code $BWA_EXIT"
            echo "  BWA log:"
            tail -10 "$BWA_LOG"
            echo "  Samtools log:"
            tail -10 "$SAMTOOLS_LOG"
        else
            echo "  ✓ Mapping complete"
            # Show summary from logs
            if [ -s "$BWA_LOG" ]; then
                echo "  BWA summary:"
                grep -E "Processed|aligned" "$BWA_LOG" | tail -3
            fi
        fi
    else
        # Conda run execution
        $GENOMICS_CMD bash -c "bwa mem -t $THREADS ${ACCESSION}.fasta $CLEANED_R1 $CLEANED_R2 | \
            samtools view -@ $THREADS -bS - | \
            samtools sort -@ $THREADS -o ${SAMPLE_NAME}_quick.bam -"
    fi
else
    echo "  Warning: Cleaned reads not found, using original reads"
    echo "  Looking for: $CLEANED_R1 and $CLEANED_R2"
    if [ -z "$GENOMICS_CMD" ]; then
        # Direct execution (Docker/micromamba base environment)
        bwa mem -t $THREADS "${ACCESSION}.fasta" "../$R1" "../$R2" 2>&1 | \
            samtools view -@ $THREADS -bS - 2>&1 | \
            samtools sort -@ $THREADS -o "${SAMPLE_NAME}_quick.bam" - 2>&1
    else
        # Conda run execution
        $GENOMICS_CMD bash -c "bwa mem -t $THREADS ${ACCESSION}.fasta ../$R1 ../$R2 | \
            samtools view -@ $THREADS -bS - | \
            samtools sort -@ $THREADS -o ${SAMPLE_NAME}_quick.bam -"
    fi
fi

# Check if BAM was created
if [ ! -f "${SAMPLE_NAME}_quick.bam" ]; then
    echo "ERROR: Failed to create BAM file from mapping"
    echo "Skipping mapping statistics..."
    TOTAL_READS=0
    MAPPED_READS=0
else
    if [ -z "$GENOMICS_CMD" ]; then
        samtools index "${SAMPLE_NAME}_quick.bam" 2>&1
    else
        $GENOMICS_CMD samtools index "${SAMPLE_NAME}_quick.bam"
    fi
fi

# Generate mapping statistics (only if BAM exists)
if [ -f "${SAMPLE_NAME}_quick.bam" ]; then
    echo "Generating mapping statistics..."
    if [ -z "$GENOMICS_CMD" ]; then
        samtools flagstat "${SAMPLE_NAME}_quick.bam" > "${SAMPLE_NAME}_mapping_stats.txt" 2>&1
        samtools idxstats "${SAMPLE_NAME}_quick.bam" > "${SAMPLE_NAME}_idxstats.txt" 2>&1
        TOTAL_READS=$(samtools view -c "${SAMPLE_NAME}_quick.bam" 2>/dev/null || echo "0")
        MAPPED_READS=$(samtools view -c -F 4 "${SAMPLE_NAME}_quick.bam" 2>/dev/null || echo "0")
    else
        $GENOMICS_CMD samtools flagstat "${SAMPLE_NAME}_quick.bam" > "${SAMPLE_NAME}_mapping_stats.txt"
        $GENOMICS_CMD samtools idxstats "${SAMPLE_NAME}_quick.bam" > "${SAMPLE_NAME}_idxstats.txt"
        TOTAL_READS=$($GENOMICS_CMD samtools view -c "${SAMPLE_NAME}_quick.bam")
        MAPPED_READS=$($GENOMICS_CMD samtools view -c -F 4 "${SAMPLE_NAME}_quick.bam")
    fi
fi

# Extract duplication rate from fastp report
# First try to get it from the original main pipeline fastp report
# Remove the _R1 suffix to get the sample base name
SAMPLE_BASE=$(echo "$R1_BASE" | sed 's/_R1$//')
ORIGINAL_FASTP_JSON="../cleaned_seqs/${SAMPLE_BASE}_fastp_report.json"
echo "  Looking for fastp report: $ORIGINAL_FASTP_JSON"
if [ -f "$ORIGINAL_FASTP_JSON" ]; then
    # Extract duplication rate from main pipeline fastp JSON
    echo "  Found fastp report, extracting duplication rate..."
    # Extract duplication rate from the "duplication": {"rate": 0.67605} structure
    RAW_DUP_RATE=$(grep -A2 '"duplication"' "$ORIGINAL_FASTP_JSON" | grep -o '"rate": *[0-9.]*' | grep -o '[0-9.]*')
    echo "  Raw duplication rate extracted: '$RAW_DUP_RATE'"
    
    if [ -n "$RAW_DUP_RATE" ] && [ "$RAW_DUP_RATE" != "0" ]; then
        DUPLICATION_RATE=$(awk "BEGIN {printf \"%.2f\", $RAW_DUP_RATE * 100}")
        echo "  Using duplication rate from main pipeline: ${DUPLICATION_RATE}%"
    else
        DUPLICATION_RATE="0.00"
        echo "  Could not extract duplication rate, using ${DUPLICATION_RATE}%"
    fi
elif [ -f "${SAMPLE_NAME}_fastp_report.json" ] && [ -s "${SAMPLE_NAME}_fastp_report.json" ]; then
    # Fallback: try to extract from local fastp JSON (if it's not empty)
    DUPLICATION_RATE=$(grep -o '"dup_rate":[0-9.]*' "${SAMPLE_NAME}_fastp_report.json" 2>/dev/null | cut -d: -f2 || echo "0")
    DUPLICATION_RATE=$(awk "BEGIN {printf \"%.2f\", $DUPLICATION_RATE * 100}")
    echo "  Using duplication rate from local fastp: ${DUPLICATION_RATE}%"
else
    # Default fallback
    DUPLICATION_RATE="0.00"
    echo "  Warning: Could not extract duplication rate, using ${DUPLICATION_RATE}%"
fi

# Calculate duplication metrics based on fastp duplication rate
if [ "$TOTAL_READS" -gt 0 ]; then
    # Use awk for calculations
    MAPPING_PERCENT=$(awk "BEGIN {printf \"%.2f\", $MAPPED_READS * 100 / $TOTAL_READS}")
    
    # Calculate unique reads based on fastp duplication rate
    DUPLICATE_READS=$(awk "BEGIN {printf \"%.0f\", $TOTAL_READS * $DUPLICATION_RATE / 100}")
    UNIQUE_READS=$((TOTAL_READS - DUPLICATE_READS))
    
    if [ "$UNIQUE_READS" -gt 0 ]; then
        # Calculate unique mapped reads correctly
        # Since BWA doesn't mark duplicates, we estimate based on the duplication rate
        MAPPED_UNIQUE_READS=$(awk "BEGIN {printf \"%.0f\", $MAPPED_READS * (100 - $DUPLICATION_RATE) / 100}")
        # Calculate the actual percentage of unique reads that mapped
        DEDUPLICATED_MAPPING_PERCENT=$(awk "BEGIN {printf \"%.2f\", $MAPPED_UNIQUE_READS * 100 / $UNIQUE_READS}")
    else
        MAPPED_UNIQUE_READS=0
        DEDUPLICATED_MAPPING_PERCENT=0
    fi
else
    MAPPING_PERCENT=0
    DUPLICATE_READS=0
    UNIQUE_READS=0
    MAPPED_UNIQUE_READS=0
    DEDUPLICATED_MAPPING_PERCENT=0
fi

echo "Mapping Results:"
echo "  Total reads: $TOTAL_READS"
echo "  Duplicate reads: $DUPLICATE_READS (${DUPLICATION_RATE}%)"
echo "  Unique reads: $UNIQUE_READS"
echo "  Mapped reads (all): $MAPPED_READS (${MAPPING_PERCENT}%)"
echo "  Mapped unique reads: $MAPPED_UNIQUE_READS (${DEDUPLICATED_MAPPING_PERCENT}% of unique)"

# Step 2: De novo assembly
echo ""
echo "Step 2: De Novo Assembly"
echo "========================"

# Check for existing cleaned reads first (using pipeline naming convention)
# Use the base filenames initialized at script start
EXISTING_R1="../cleaned_seqs/${R1_BASE}.qc.fastq.gz"
EXISTING_R2="../cleaned_seqs/${R2_BASE}.qc.fastq.gz"

if [ -f "$EXISTING_R1" ] && [ -f "$EXISTING_R2" ]; then
    echo "Found existing cleaned reads from main pipeline, reusing them..."
    echo "  Using: $EXISTING_R1"
    echo "  Using: $EXISTING_R2"
    
    # Create symlinks using the same base filenames for consistency
    ln -sf "$EXISTING_R1" "${SAMPLE_NAME}_R1.qc.fastq.gz"
    ln -sf "$EXISTING_R2" "${SAMPLE_NAME}_R2.qc.fastq.gz"
    
    # Create a placeholder fastp report
    echo "Reused existing cleaned reads from main pipeline" > "${SAMPLE_NAME}_fastp_reused.txt"
    touch "${SAMPLE_NAME}_fastp_report.html"
    touch "${SAMPLE_NAME}_fastp_report.json"
    
else
    echo "No existing cleaned reads found, cleaning reads with fastp..."
    echo "  Looking for: $EXISTING_R1"
    echo "  Looking for: $EXISTING_R2"
    
    # Clean reads with identical settings to main pipeline
    $GENOMICS_CMD fastp -i "../$R1" -I "../$R2" \
          -o "${SAMPLE_NAME}_R1.qc.fastq.gz" \
          -O "${SAMPLE_NAME}_R2.qc.fastq.gz" \
          -h "${SAMPLE_NAME}_fastp_report.html" \
          -j "${SAMPLE_NAME}_fastp_report.json" \
          --thread "$THREADS" \
          --qualified_quality_phred 15 \
          --unqualified_percent_limit 40 \
          --length_required 36
fi

# Run MEGAHIT assembly
echo "Running MEGAHIT assembly..."

# Remove existing assembly directory if it exists
if [ -d "assembly_${SAMPLE_NAME}" ]; then
    echo "Removing existing assembly directory..."
    rm -rf "assembly_${SAMPLE_NAME}"
fi

# MEGAHIT is available in current environment (no separate activation needed in Docker)
# Build MEGAHIT command with conditional extreme memory settings
MEGAHIT_CMD="megahit -1 \"${SAMPLE_NAME}_R1.qc.fastq.gz\" -2 \"${SAMPLE_NAME}_R2.qc.fastq.gz\" -o \"assembly_${SAMPLE_NAME}\" --presets meta-sensitive --min-contig-len 500 --memory 0.7 -t \"$THREADS\""

if [ -n "$EXTREME_MEMORY_FLAG" ]; then
    echo "Adding extreme memory settings for MEGAHIT..."
    MEGAHIT_CMD="$MEGAHIT_CMD --memory 0.9"  # Use 90% of available memory
fi

echo "MEGAHIT command: $MEGAHIT_CMD"
eval $MEGAHIT_CMD

# Check if assembly succeeded
if [ ! -f "assembly_${SAMPLE_NAME}/final.contigs.fa" ]; then
    echo "ERROR: Assembly failed - no contigs produced"
    exit 1
fi

# Get assembly statistics
CONTIG_COUNT=$(grep -c ">" "assembly_${SAMPLE_NAME}/final.contigs.fa")
echo "Assembly completed: $CONTIG_COUNT contigs generated"

# Step 3: Viral BLAST analysis
echo ""
echo "Step 3: Viral Contamination BLAST"
echo "================================="

# Filter and sort contigs by length (>1000bp, largest first)
echo "Filtering contigs >1000bp and sorting by size for BLAST analysis..."
seqkit seq -m 1000 "assembly_${SAMPLE_NAME}/final.contigs.fa" | \
seqkit sort --by-length --reverse > "${SAMPLE_NAME}_contigs_filtered.fa"

FILTERED_COUNT=$(grep -c ">" "${SAMPLE_NAME}_contigs_filtered.fa")
echo "Filtered contigs for BLAST: $FILTERED_COUNT"

if [ "$FILTERED_COUNT" -gt 0 ]; then
    # Check for local contamination database first
    # Prefer vicast_combined alias (supports user extensions), fallback to base DB
    LOCAL_DB_PATH="${BLAST_DB:-/opt/vicast/blast_db/vicast_combined}"
    if [[ ! -f "${LOCAL_DB_PATH}.nal" ]] && [[ ! -f "${LOCAL_DB_PATH}.nhr" ]]; then
        LOCAL_DB_PATH="${BLAST_DB_DIR:-/opt/vicast/blast_db}/microbial_contaminants"
    fi
    
    # Add header to BLAST results
    echo -e "Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\tMismatches\tGap_Opens\tQuery_Start\tQuery_End\tSubject_Start\tSubject_End\tE_value\tBit_Score\tSubject_Title" > "${SAMPLE_NAME}_blast_all.tsv"

    # BLAST is available in current environment (no separate activation needed in Docker)
    echo "Running BLAST analysis..."
    
    if [ -n "$LOCAL_DB_PATH" ] && { [ -f "${LOCAL_DB_PATH}.nal" ] || [ -f "${LOCAL_DB_PATH}.nhr" ] || [ -f "${LOCAL_DB_PATH}.00.nhr" ]; }; then
        echo "Using local contamination database for fast BLAST..."
        echo "Note: Contigs are sorted by size (largest first) for priority analysis"
        
        blastn -query "${SAMPLE_NAME}_contigs_filtered.fa" \
               -db "$LOCAL_DB_PATH" \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
               -max_target_seqs 5 \
               -max_hsps 1 \
               -evalue 1e-10 \
               -num_threads "$THREADS" \
               >> "${SAMPLE_NAME}_blast_all.tsv"
    else
        echo "Local database not found, using remote BLAST (slower, may not process all contigs)..."
        echo "Note: Contigs are sorted by size (largest first) for priority analysis"
        echo "WARNING: Remote BLAST may timeout with many contigs. Consider building local database."
        
        # For remote BLAST, only process top 20 contigs to avoid timeout
        head -40 "${SAMPLE_NAME}_contigs_filtered.fa" > "${SAMPLE_NAME}_top20_contigs.fa"
        CONTIG_COUNT_TOP20=$(grep -c ">" "${SAMPLE_NAME}_top20_contigs.fa")
        echo "Processing top $CONTIG_COUNT_TOP20 contigs only (to avoid timeout)..."
        
        blastn -query "${SAMPLE_NAME}_top20_contigs.fa" \
               -remote -db nt \
               -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
               -max_target_seqs 5 \
               -max_hsps 1 \
               -evalue 1e-10 \
               >> "${SAMPLE_NAME}_blast_all.tsv"
        
        rm -f "${SAMPLE_NAME}_top20_contigs.fa"
    fi
    
    # Parse BLAST results with Python (handles tied hits, accession matching, kingdom classification)
    echo "Parsing BLAST results..."
    if [ -s "${SAMPLE_NAME}_blast_all.tsv" ] && [ $(wc -l < "${SAMPLE_NAME}_blast_all.tsv") -gt 1 ]; then
        python "${PIPELINE_DIR}/parse_blast_results.py" \
            "${SAMPLE_NAME}_blast_all.tsv" \
            "${SAMPLE_NAME}_contigs_filtered.fa" \
            --accession "$ACCESSION" \
            --top-hits "${SAMPLE_NAME}_top_hits.tsv" \
            --viral-hits "${SAMPLE_NAME}_viral_blast.tsv" \
            --report "${SAMPLE_NAME}_blast_report_section.txt"

        TOP_HIT_COUNT=$(tail -n +2 "${SAMPLE_NAME}_top_hits.tsv" 2>/dev/null | wc -l | tr -d ' ')
        VIRAL_HIT_COUNT=$(tail -n +2 "${SAMPLE_NAME}_viral_blast.tsv" 2>/dev/null | wc -l | tr -d ' ')
        echo "Top hits: $TOP_HIT_COUNT entries (best hits per contig)"
        echo "Viral hits: $VIRAL_HIT_COUNT"
    else
        echo -e "Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\tAlignment_Length\tQuery_Coverage\tE_value\tBit_Score\tSubject_Title" > "${SAMPLE_NAME}_viral_blast.tsv"
        echo -e "Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\tAlignment_Length\tQuery_Coverage\tE_value\tBit_Score\tKingdom/Type\tSubject_Title" > "${SAMPLE_NAME}_top_hits.tsv"
        echo "No BLAST results obtained"
    fi
    
    # Alternative: Local viral database if available
    # blastn -query "${SAMPLE_NAME}_contigs_filtered.fa" \
    #        -db /ref/common/blastdb/viral_nt \
    #        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" \
    #        -max_target_seqs 5 \
    #        -max_hsps 1 \
    #        -evalue 1e-5 \
    #        > "${SAMPLE_NAME}_viral_blast.tsv"
    
    echo "BLAST analysis completed"
else
    echo "No contigs >1000bp found - skipping BLAST analysis"
    touch "${SAMPLE_NAME}_viral_blast.tsv"
fi

# Step 4: Generate diagnostic report
echo ""
echo "Step 4: Generating Diagnostic Report"
echo "===================================="

cat > "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
========================================
VIRAL CONTAMINATION DIAGNOSTIC REPORT
========================================
Sample: $SAMPLE_NAME
Reference: $ACCESSION
Analysis Date: $(date)
Working Directory: $(pwd)

========================================
MAPPING STATISTICS
========================================
Total Reads: $TOTAL_READS
Duplicate Reads: $DUPLICATE_READS (${DUPLICATION_RATE}%)
Unique Reads: $UNIQUE_READS

Raw Mapping: $MAPPED_READS reads (${MAPPING_PERCENT}% of total)
Deduplicated Mapping: $MAPPED_UNIQUE_READS reads (${DEDUPLICATED_MAPPING_PERCENT}% of unique)

Interpretation (based on deduplicated mapping):
- >70%: Good mapping to reference, likely correct organism
- 30-70%: Moderate mapping, possible mixed infection or variant
- <30%: Poor mapping, likely wrong reference or contamination

Note: High duplication rates (>50%) are expected in deep sequencing
of cell cultures for mutation detection and do not indicate quality issues.

========================================
ASSEMBLY STATISTICS
========================================
Total Contigs: $CONTIG_COUNT
Contigs >1000bp: $FILTERED_COUNT

========================================
BLAST RESULTS SUMMARY
========================================
EOF

# Add top hits summary to report (from parse_blast_results.py --report output)
if [ -s "${SAMPLE_NAME}_blast_report_section.txt" ]; then
    cat "${SAMPLE_NAME}_blast_report_section.txt" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
elif [ -s "${SAMPLE_NAME}_top_hits.tsv" ] && [ $(wc -l < "${SAMPLE_NAME}_top_hits.tsv") -gt 1 ]; then
    # Fallback: simple top hits table
    echo "TOP HITS BY CONTIG:" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    tail -n +2 "${SAMPLE_NAME}_top_hits.tsv" | head -20 >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
fi

# Viral contamination analysis section (already covered by parse_blast_results.py report)
if [ ! -s "${SAMPLE_NAME}_viral_blast.tsv" ] || [ $(tail -n +2 "${SAMPLE_NAME}_viral_blast.tsv" 2>/dev/null | wc -l | tr -d ' ') -eq 0 ]; then
    echo "VIRAL CONTAMINATION ANALYSIS:" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "No viral contamination detected in assembled contigs." >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "This indicates:" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "  - Clean viral culture (expected result)" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "  - No detectable co-infection with other viruses" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "  - Assembly contigs match target reference genome" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
fi

# Add file references with full paths
DIAG_DIR=$(pwd)
cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF

DETAILED RESULTS FILES:
- Top hits summary: ${DIAG_DIR}/${SAMPLE_NAME}_top_hits.tsv
- All BLAST results: ${DIAG_DIR}/${SAMPLE_NAME}_blast_all.tsv
- Viral BLAST results: ${DIAG_DIR}/${SAMPLE_NAME}_viral_blast.tsv
- Assembly contigs: ${DIAG_DIR}/assembly_${SAMPLE_NAME}/final.contigs.fa
- HTML report: ${DIAG_DIR}/diagnostic_${SAMPLE_NAME}_presentation_ready_report.html
EOF

if [ ! -s "${SAMPLE_NAME}_blast_all.tsv" ] || [ $(wc -l < "${SAMPLE_NAME}_blast_all.tsv") -le 1 ]; then
    echo "" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "No significant BLAST hits found (E-value < 1e-10)" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "This may indicate:" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "  - Clean viral culture (no bacterial/fungal contamination)" >> "${SAMPLE_NAME}_diagnostic_report.txt"
    echo "  - Low contig complexity or very small contigs" >> "${SAMPLE_NAME}_diagnostic_report.txt"
fi

# Add recommendations
cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF

========================================
RECOMMENDATIONS
========================================
EOF

# Generate recommendations based on deduplicated mapping percentage
if (( $(awk "BEGIN {print ($DEDUPLICATED_MAPPING_PERCENT < 30)}") )); then
    cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
WARNING: Low deduplicated mapping percentage (${DEDUPLICATED_MAPPING_PERCENT}%)
- Check BLAST results for actual organism identity
- Consider using different reference genome
- Investigate potential sample contamination or mislabeling
- Raw mapping: ${MAPPING_PERCENT}% (includes ${DUPLICATION_RATE}% duplicates)
EOF
elif (( $(awk "BEGIN {print ($DEDUPLICATED_MAPPING_PERCENT < 70)}") )); then
    cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
CAUTION: Moderate deduplicated mapping percentage (${DEDUPLICATED_MAPPING_PERCENT}%)
- Review BLAST results for mixed infections
- Check for strain variants or recombinants
- Consider manual review of mapping results
- Raw mapping: ${MAPPING_PERCENT}% (includes ${DUPLICATION_RATE}% duplicates)
EOF
else
    cat >> "${SAMPLE_NAME}_diagnostic_report.txt" << EOF
GOOD: High deduplicated mapping percentage (${DEDUPLICATED_MAPPING_PERCENT}%)
- Mapping suggests correct reference genome
- Proceed with standard analysis pipeline
- Monitor for any unusual variant patterns
- Raw mapping: ${MAPPING_PERCENT}% (includes ${DUPLICATION_RATE}% duplicates)
EOF
fi

# Step 5: Cleanup and final summary
echo ""
echo "Cleaning up temporary files..."
rm -f "${ACCESSION}.fasta."* # Remove BWA index files
# Keep BAM files for QC review - they're useful for debugging
# rm -f "${SAMPLE_NAME}_quick.bam" "${SAMPLE_NAME}_quick.bam.bai"

echo ""
echo "========================================="
echo "GENERATING QC VISUALIZATION REPORT"
echo "========================================="

# Generate presentation-ready QC report
# PIPELINE_DIR already resolved at top of script
QC_SCRIPT_PATH="${PIPELINE_DIR}/qc_with_simple_plots.py"

if [ -f "$QC_SCRIPT_PATH" ]; then
    echo "Creating QC visualization report..."
    # Run QC script from parent directory so it can find diagnostic_${SAMPLE_NAME}
    cd ..
    python "$QC_SCRIPT_PATH" "diagnostic_${SAMPLE_NAME}"
    QC_EXIT_CODE=$?
    cd "diagnostic_${SAMPLE_NAME}"
    
    if [ $QC_EXIT_CODE -eq 0 ]; then
        echo "✅ QC report generated: diagnostic_${SAMPLE_NAME}/diagnostic_${SAMPLE_NAME}_presentation_ready_report.html"
        echo "🌐 Open this file in a web browser for presentation-ready results"
    else
        echo "⚠️  QC report generation failed - check diagnostic output files"
    fi
else
    echo "⚠️  QC script not found at: $QC_SCRIPT_PATH"
    echo "   Manual QC generation: python ${PIPELINE_DIR}/qc_with_simple_plots.py diagnostic_${SAMPLE_NAME}"
fi

# Map Docker /data path to host path if HOST_PWD is set
DISPLAY_DIR="${DIAGNOSTIC_DIR}"
if [[ -n "${HOST_PWD:-}" ]]; then
    ABS_DIAG_DIR="$(cd "${DIAGNOSTIC_DIR}" 2>/dev/null && pwd || echo "${DIAGNOSTIC_DIR}")"
    if [[ "$ABS_DIAG_DIR" == /data* ]]; then
        DISPLAY_DIR="${HOST_PWD}${ABS_DIAG_DIR#/data}"
    fi
fi

echo ""
echo "========================================="
echo "DIAGNOSTIC ANALYSIS COMPLETE"
echo "========================================="
echo ""
echo "📊 REVIEW QC OUTPUTS:"
echo ""
echo "1. HTML Report (open in browser):"
echo "   ${DISPLAY_DIR}/diagnostic_${SAMPLE_NAME}_presentation_ready_report.html"
echo ""
echo "2. Text Report:"
echo "   ${DISPLAY_DIR}/${SAMPLE_NAME}_diagnostic_report.txt"
echo ""
echo "3. BLAST Results:"
echo "   ${DISPLAY_DIR}/${SAMPLE_NAME}_viral_blast.tsv"
echo "   ${DISPLAY_DIR}/${SAMPLE_NAME}_top_hits.tsv"
echo ""
echo "4. Mapping Quality:"
echo "   Deduplicated mapping: ${DEDUPLICATED_MAPPING_PERCENT}% (target: >70%)"
echo "   Total contigs >1kb: ${FILTERED_COUNT}"
echo ""
echo "========================================="
echo "📋 INTERPRETATION GUIDANCE"
echo "========================================="
echo ""
if (( $(awk "BEGIN {print ($DEDUPLICATED_MAPPING_PERCENT >= 70)}") )); then
    echo "✅ GOOD QUALITY - Proceed with variant analysis"
    echo ""
    echo "Next step: Review HTML report and continue to annotation (Steps 7-9)"
elif (( $(awk "BEGIN {print ($DEDUPLICATED_MAPPING_PERCENT >= 30)}") )); then
    echo "⚠️  MODERATE QUALITY - Review recommended"
    echo ""
    echo "Action: Check BLAST results for mixed infection or contamination"
    echo "        Consider manual review before proceeding"
else
    echo "❌ LOW QUALITY - Investigation required"
    echo ""
    echo "Action: Check BLAST results for organism identity"
    echo "        May need different reference genome"
    echo "        Review diagnostic report for details"
fi
echo ""
echo "========================================="
echo "🔬 CONTINUE TO ANNOTATION (Steps 7-9)"
echo "========================================="
echo ""
echo "After reviewing QC results, run annotation workflow:"
echo ""
echo "Command:"
echo "  cd .."
echo "  bash ${PIPELINE_DIR}/viral_diagnostic.sh \\"
echo "    ${R1} ${R2} ${ACCESSION} ${SAMPLE_NAME} ${THREADS}"
echo ""
echo "Or from parent directory:"
echo "  run_vicast_analyze_annotate_only.sh ${R1} ${R2} ${ACCESSION}"
echo ""
echo "Variant filtering parameters (modify in viral_pipeline.py if needed):"
echo "  Step 7 (high confidence): freq ≥1%, depth ≥200×, qual ≥1000"
echo "  Step 8 (low frequency):   freq 0.5-5%, depth ≥200×, qual ≥1000"
echo ""
echo "========================================="
echo ""
echo "Diagnostic analysis completed at: $(date)"
echo ""
echo "  - assembly_${SAMPLE_NAME}/final.contigs.fa (assembled contigs - sorted by size)"
