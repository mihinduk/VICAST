#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --job-name=segmented_virus
#SBATCH --output=%x_%j.out

# Process Segmented Virus Pipeline
# Runs variant calling and consensus generation for all segments
# Usage: sbatch process_segmented_virus.sh <R1> <R2> <virus_id>

if [ $# -lt 3 ]; then
    cat << USAGE
Usage: $0 <R1_fastq> <R2_fastq> <virus_id>

Example:
  sbatch process_segmented_virus.sh \\
    influenza_R1.fastq.gz \\
    influenza_R2.fastq.gz \\
    influenza_pr8

The virus_id must be in known_viruses.json with 'segment_accessions' field.

Output will be created in: [sample_name]_[virus_id]_results/
  - One subdirectory per segment
  - Consolidated proteins and genomes
  - Summary report
USAGE
    exit 1
fi

R1=$1
R2=$2
VIRUS_ID=$3
THREADS=${SLURM_CPUS_PER_TASK:-8}

# Extract sample name from R1
SAMPLE_NAME=$(basename "$R1" | sed -E 's/_R?[12](_[0-9]+)?\..*$//')
OUTPUT_DIR="${SAMPLE_NAME}_${VIRUS_ID}_results"

echo "========================================="
echo "SEGMENTED VIRUS PIPELINE"
echo "========================================="
echo "Virus: $VIRUS_ID"
echo "Sample: $SAMPLE_NAME"
echo "R1: $R1"
echo "R2: $R2"
echo "Output: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo "Started: $(date)"
echo "========================================="

# Setup paths
if [ -n "$SLURM_JOB_ID" ]; then
    ORIGINAL_SCRIPT=$(scontrol show job $SLURM_JOB_ID | grep -oP 'Command=\K[^ ]+' || echo "")
    if [ -n "$ORIGINAL_SCRIPT" ]; then
        PIPELINE_DIR="$(dirname "$ORIGINAL_SCRIPT")"
    else
        PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi
else
    PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

KNOWN_VIRUSES="${PIPELINE_DIR}/known_viruses.json"

# Activate environment
source /ref/sahlab/software/anaconda3/bin/activate
conda activate vicast_analyze

# Extract segment information
python3 << PYEOF
import json
import sys

with open("$KNOWN_VIRUSES", 'r') as f:
    viruses = json.load(f)

if "$VIRUS_ID" not in viruses:
    print(f"ERROR: $VIRUS_ID not in known_viruses.json", file=sys.stderr)
    sys.exit(1)

virus_data = viruses["$VIRUS_ID"]

if "segment_accessions" not in virus_data:
    print(f"ERROR: $VIRUS_ID doesn't have segment_accessions", file=sys.stderr)
    print("Use standard pipeline for non-segmented viruses", file=sys.stderr)
    sys.exit(1)

segment_accessions = virus_data["segment_accessions"]
print(f"Found {len(segment_accessions)} segments:", list(segment_accessions.keys()))

# Write segment file
with open("/tmp/segments_$VIRUS_ID.txt", 'w') as out:
    for seg_name, accession in segment_accessions.items():
        out.write(f"{seg_name}\t{accession}\n")
PYEOF

if [ $? -ne 0 ]; then
    echo "Failed to read segment information"
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit 1

# Step 1: Quality control (once for all segments)
echo ""
echo "Step 1: Quality Control (shared across segments)"
echo "========================================="

QC_DIR="../cleaned_seqs"
mkdir -p "$QC_DIR"

QC_R1="${QC_DIR}/${SAMPLE_NAME}_R1.qc.fastq.gz"
QC_R2="${QC_DIR}/${SAMPLE_NAME}_R2.qc.fastq.gz"

if [ -f "$QC_R1" ] && [ -f "$QC_R2" ]; then
    echo "Found existing cleaned reads, reusing..."
else
    echo "Cleaning reads with fastp..."
    fastp -i "../$R1" -I "../$R2" \
          -o "$QC_R1" -O "$QC_R2" \
          -h "${SAMPLE_NAME}_fastp.html" \
          -j "${SAMPLE_NAME}_fastp.json" \
          --detect_adapter_for_pe \
          --correction \
          --thread "$THREADS"
fi

# Step 2: Process each segment
echo ""
echo "Step 2: Processing Segments"
echo "========================================="

SEGMENT_NUM=0
FAILED_SEGMENTS=()

while IFS=$'\t' read -r seg_name accession; do
    SEGMENT_NUM=$((SEGMENT_NUM + 1))
    
    echo ""
    echo "[$SEGMENT_NUM] Segment: $seg_name ($accession)"
    echo "----------------------------------------"
    
    mkdir -p "${seg_name}"
    
    # Download reference
    REF_PATH="${QC_DIR}/${accession}.fasta"
    if [ ! -f "$REF_PATH" ]; then
        echo "  Downloading reference..."
        python3 << PYEOF
from Bio import Entrez, SeqIO
Entrez.email = "vicast@example.com"
try:
    handle = Entrez.efetch(db="nucleotide", id="${accession}", rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    record.id = "${accession}"
    record.description = ""
    SeqIO.write(record, "${REF_PATH}", "fasta")
    handle.close()
    print(f"Downloaded {len(record.seq)} bp")
except Exception as e:
    print(f"ERROR downloading ${accession}: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
        if [ $? -ne 0 ]; then
            echo "  ✗ Failed to download reference"
            FAILED_SEGMENTS+=("$seg_name")
            continue
        fi
    fi
    
    # Align to reference
    echo "  Aligning reads..."
    bwa index "$REF_PATH" 2>/dev/null
    
    bwa mem -t "$THREADS" "$REF_PATH" "$QC_R1" "$QC_R2" 2>/dev/null | \
    samtools view -@ "$THREADS" -bS - | \
    samtools sort -@ "$THREADS" -o "${seg_name}/${seg_name}.bam" -
    
    samtools index "${seg_name}/${seg_name}.bam"
    
    # Variant calling
    echo "  Calling variants..."
    lofreq call-parallel --pp-threads "$THREADS" \
        -f "$REF_PATH" \
        -o "${seg_name}/${seg_name}.vcf" \
        "${seg_name}/${seg_name}.bam" 2>/dev/null
    
    # Annotate with SnpEff
    echo "  Annotating variants..."
    snpEff -noStats "${accession}" "${seg_name}/${seg_name}.vcf" \
        > "${seg_name}/${seg_name}.ann.vcf" 2>/dev/null
    
    # Convert to TSV
    python3 "${PIPELINE_DIR}/viral_pipeline/visualization/vcf_to_tsv_annotated.py" \
        "${seg_name}/${seg_name}.ann.vcf" \
        "${seg_name}/${seg_name}.ann.tsv" 2>/dev/null
    
    # Parse and filter mutations
    echo "  Filtering mutations..."
    python3 "${PIPELINE_DIR}/viral_pipeline/visualization/parse_snpeff_tsv.py" \
        "${seg_name}/${seg_name}.ann.tsv" \
        "${seg_name}/${seg_name}_filtered_mutations.tsv" \
        --quality 1000 --depth 20 --freq 0.01 2>/dev/null
    
    # Generate consensus
    echo "  Generating consensus..."
    python3 "${PIPELINE_DIR}/generate_filtered_consensus.py" \
        --vcf "${seg_name}/${seg_name}_filtered_mutations.tsv" \
        --reference "$REF_PATH" \
        --accession "${accession}" \
        --quality 1000 --depth 20 --freq 0.50 \
        --output-prefix "${seg_name}/${seg_name}_consensus" 2>&1 | grep -E "(✅|⚠️|mutations|proteins)"
    
    if [ -f "${seg_name}/${seg_name}_consensus.fasta" ]; then
        echo "  ✓ Consensus generated"
    else
        echo "  ✗ Consensus generation failed"
        FAILED_SEGMENTS+=("$seg_name")
    fi
    
done < /tmp/segments_${VIRUS_ID}.txt

# Step 3: Consolidate results
echo ""
echo "Step 3: Consolidating Results"
echo "========================================="

# All segments concatenated
cat /dev/null > "${VIRUS_ID}_all_segments.fasta"
for seg_dir in */; do
    [ -d "$seg_dir" ] || continue
    seg_name="${seg_dir%/}"
    if [ -f "${seg_dir}${seg_name}_consensus.fasta" ]; then
        cat "${seg_dir}${seg_name}_consensus.fasta" >> "${VIRUS_ID}_all_segments.fasta"
    fi
done

# All proteins combined  
cat /dev/null > "${VIRUS_ID}_all_proteins.fasta"
for seg_dir in */; do
    [ -d "$seg_dir" ] || continue
    seg_name="${seg_dir%/}"
    if [ -f "${seg_dir}${seg_name}_consensus_proteins.fasta" ]; then
        cat "${seg_dir}${seg_name}_consensus_proteins.fasta" >> "${VIRUS_ID}_all_proteins.fasta"
    fi
done

TOTAL_PROTEINS=$(grep -c "^>" "${VIRUS_ID}_all_proteins.fasta" 2>/dev/null || echo "0")

# Summary report
cat > "${VIRUS_ID}_summary_report.txt" << REPORT
========================================
SEGMENTED VIRUS ANALYSIS SUMMARY
========================================
Virus: $VIRUS_ID
Sample: $SAMPLE_NAME
Analysis Date: $(date)
Total Segments: $SEGMENT_NUM
Total Proteins: $TOTAL_PROTEINS

SEGMENT RESULTS:
========================================
REPORT

for seg_dir in */; do
    [ -d "$seg_dir" ] || continue
    seg_name="${seg_dir%/}"
    
    if [ -f "${seg_dir}${seg_name}_consensus_summary_report.txt" ]; then
        mut_count=$(grep -oP "Total mutations passing filters: \K\d+" "${seg_dir}${seg_name}_consensus_summary_report.txt" || echo "0")
        prot_count=$(grep -c "^>" "${seg_dir}${seg_name}_consensus_proteins.fasta" 2>/dev/null || echo "0")
        
        echo "✓ $seg_name: $mut_count mutations, $prot_count proteins" >> "${VIRUS_ID}_summary_report.txt"
        
        # Add detailed mutations
        if [ -f "${seg_dir}${seg_name}_consensus_summary_report.txt" ]; then
            echo "" >> "${VIRUS_ID}_summary_report.txt"
            grep -A 30 "MUTATIONS PER PROTEIN" "${seg_dir}${seg_name}_consensus_summary_report.txt" | head -25 >> "${VIRUS_ID}_summary_report.txt"
        fi
    else
        echo "✗ $seg_name: FAILED" >> "${VIRUS_ID}_summary_report.txt"
    fi
done

cat >> "${VIRUS_ID}_summary_report.txt" << REPORT

OUTPUT FILES:
========================================
- ${VIRUS_ID}_all_segments.fasta (all 8 consensus genomes)
- ${VIRUS_ID}_all_proteins.fasta ($TOTAL_PROTEINS proteins from all segments)
- Individual segment results in subdirectories (PB2/, PB1/, etc.)
REPORT

echo ""
echo "========================================="
echo "SEGMENTED VIRUS PIPELINE COMPLETE"
echo "========================================="
echo "Virus: $VIRUS_ID ($SEGMENT_NUM segments)"
echo "Sample: $SAMPLE_NAME"
echo "Total Proteins: $TOTAL_PROTEINS"

if [ ${#FAILED_SEGMENTS[@]} -gt 0 ]; then
    echo "⚠️  Failed segments: ${FAILED_SEGMENTS[*]}"
fi

echo ""
echo "Key Output Files:"
echo "  - ${VIRUS_ID}_all_segments.fasta"
echo "  - ${VIRUS_ID}_all_proteins.fasta"
echo "  - ${VIRUS_ID}_summary_report.txt"
echo ""
echo "Completed: $(date)"
echo "========================================="
