# VICAST Contamination Screening Guide

**Module:** viral_diagnostic.sh
**Version:** 2.2.0
**Last Updated:** 2026-02-05

---

## Table of Contents

1. [Overview](#overview)
2. [When and Why to Use](#when-and-why-to-use)
3. [De Novo Assembly Approach](#de-novo-assembly-approach)
4. [BLAST Database Setup](#blast-database-setup)
5. [Running the Pipeline](#running-the-pipeline)
6. [Interpreting Results](#interpreting-results)
7. [Publication Reporting](#publication-reporting)
8. [Troubleshooting](#troubleshooting)

---

## Overview

Contamination screening is a critical quality control step for viral genomics studies, especially for **publication-quality analysis**. VICAST's `viral_diagnostic.sh` module provides comprehensive contamination detection through:

- De novo assembly with MEGAHIT
- BLAST-based taxonomic classification
- Viral, bacterial, and fungal screening
- Automated reporting with coverage metrics

### What is Detected?

| Contamination Type | Detection Method | Reported In |
|-------------------|------------------|-------------|
| **Other viruses** | BLAST against viral database | viral_hits.txt |
| **Bacteria** | BLAST against microbial database | bacterial_hits.txt |
| **Fungi** | BLAST against microbial database | fungal_hits.txt |
| **Human/host** | BLAST against microbial database | top_hits.tsv |

### Why This Matters

**Scientific integrity:**
- Distinguishes true viral sequences from contamination
- Identifies co-infections that could confound results
- Validates that sequencing targeted the intended virus

**Publication requirements:**
- Many journals require contamination screening
- Reviewers often ask: "How do you know it's not contamination?"
- Strengthens claims about viral evolution and adaptation

---

## When and Why to Use

### Always Use For

✅ **Publications** - Essential for peer review
✅ **Novel findings** - Unexpected mutations or phenotypes
✅ **Clinical samples** - Higher contamination risk
✅ **New viruses** - Validate identity
✅ **Failed experiments** - Understand what went wrong

### Optional For

⚠️ **Well-characterized lab strains** - Lower risk, but still recommended
⚠️ **Routine passage tracking** - If previous passages were clean
⚠️ **Internal pilot studies** - Can defer until publication

### Critical Scenarios

**Scenario 1: Unexpected Phenotype**
```
Problem: Virus suddenly shows increased pathogenicity at passage 5
Question: Is this due to viral evolution or bacterial contamination?
Solution: Run contamination screening on passage 5 sample
```

**Scenario 2: Low Mapping Rate**
```
Problem: Only 30% of reads map to expected virus reference
Question: What are the other 70% of reads?
Solution: De novo assembly + BLAST reveals co-infection
```

**Scenario 3: Publication Submission**
```
Problem: Reviewer asks "How do you know this is authentic?"
Question: Can you prove no contamination?
Solution: Include contamination screening in supplementary materials
```

---

## De Novo Assembly Approach

### Why De Novo Assembly?

**Reference-independent detection:**
- Doesn't require knowing what contaminant to look for
- Detects unexpected organisms
- Assembles contigs long enough for confident BLAST hits

**Sensitivity:**
- Even minor contaminants (1-5% of reads) produce assembled contigs
- Longer contigs (>1kb) give unambiguous taxonomy

### Assembly Strategy

**MEGAHIT meta-sensitive mode:**
```bash
megahit \
    -1 R1.qc.fastq.gz \
    -2 R2.qc.fastq.gz \
    -o assembly_sample \
    --presets meta-sensitive \
    --min-contig-len 500 \
    --memory 0.7 \
    -t 4
```

**Parameters:**
- `--presets meta-sensitive` - Optimized for mixed samples
- `--min-contig-len 500` - Only keep substantial contigs
- `--memory 0.7` - Use 70% of available RAM (adjustable)

**Why meta-sensitive?**
- Designed for metagenomes (mixed organisms)
- More sensitive than single-genome mode
- Better at assembling minor contaminants

### What Gets Assembled

**Expected outputs:**

| Sample Quality | Expected Contigs | Dominant Organism | Minor Contaminants |
|----------------|------------------|-------------------|-------------------|
| Clean virus | 1-10 | Reference virus | None or trace host |
| Bacterial contamination | 100-1000 | Mix of virus + bacteria | Multiple species |
| Co-infection | 10-100 | Multiple viruses | Depends on abundance |
| Failed prep | 1000+ | Host genome | Multiple contaminants |

---

## BLAST Database Setup

### Option 1: Local Database (Recommended)

**Advantages:**
- Fast (minutes instead of hours)
- No internet dependency
- Can process unlimited samples
- Complete results (no timeout)

**Setup:**

```bash
# 1. Download or create microbial database
# Example: NCBI microbial genomes
mkdir -p /ref/sahlab/data/microbes_db

cd /ref/sahlab/data/microbes_db

# Download microbial genomes (bacteria, fungi, viruses)
# Option A: From NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_viruses_rep_genomes.tar.gz
tar -xzf ref_viruses_rep_genomes.tar.gz

# Option B: Build custom database
cat viruses.fasta bacteria.fasta fungi.fasta > microbial_contaminants.fasta

makeblastdb -in microbial_contaminants.fasta \
            -dbtype nucl \
            -out microbial_contaminants \
            -title "Microbial Contaminants Database"

# 2. Set environment variable
export BLAST_DB=/ref/sahlab/data/microbes_db/microbial_contaminants

# 3. Make permanent
echo 'export BLAST_DB=/ref/sahlab/data/microbes_db/microbial_contaminants' >> ~/.bashrc
```

**Verify:**
```bash
# Check database
blastdbcmd -db $BLAST_DB -info

# Should show: number of sequences, database size, etc.
```

### Option 2: Remote BLAST (Fallback)

**Used automatically if local database not found:**
- Searches against NCBI nt database
- **Slow** - can take hours for many contigs
- **Limited** - processes only top 20 contigs (timeout protection)
- Internet connection required

**When to use:**
- Quick check on single sample
- No local database available
- Small number of contigs (<20)

---

## Running the Pipeline

### Basic Usage

```bash
cd /path/to/your/data

/path/to/VICAST/vicast-analyze/viral_diagnostic.sh \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    NC_001477 \
    sample_name \
    4  # threads
```

**Arguments:**
1. **R1_fastq** - Forward reads (paired-end)
2. **R2_fastq** - Reverse reads
3. **accession** - Expected virus reference (e.g., NC_001477)
4. **sample_name** - Sample identifier for outputs
5. **threads** - Number of CPU cores (default: 4)

### As Part of QC Workflow

**The diagnostic module runs automatically in:**

```bash
# Runs as Step 6
./run_vicast_analyze_qc_only.sh R1.fq.gz R2.fq.gz NC_001477 4
```

**Or run standalone:**

```bash
# Independent contamination check
./viral_diagnostic.sh R1.fq.gz R2.fq.gz NC_001477 my_sample 8
```

### For Large Files

**If you have very high coverage (>10GB FASTQ):**

```bash
./viral_diagnostic.sh \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    NC_001477 \
    sample_name \
    8 \
    --extremely-large-files  # Increases memory allocation
```

**Effect:**
- MEGAHIT uses 90% memory instead of 70%
- Allows assembly of very deep sequencing data

### On HPC/SLURM

**Submit as job:**

```bash
#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --job-name=contamination_check

# Load conda
source activate vicast_analyze

# Run diagnostic
/path/to/VICAST/vicast-analyze/viral_diagnostic.sh \
    $R1 $R2 $ACCESSION $SAMPLE 8
```

**Submit:**
```bash
sbatch contamination_check.slurm
```

---

## Interpreting Results

### Output Directory Structure

```
diagnostic_sample_name/
├── sample_name_diagnostic_report.txt          # Main summary report
├── diagnostic_sample_name_presentation_ready_report.html  # HTML report
│
├── Mapping Statistics
├── sample_name_mapping_stats.txt              # Samtools flagstat
├── sample_name_idxstats.txt                   # Per-contig mapping
│
├── Assembly
├── assembly_sample_name/final.contigs.fa      # All assembled contigs
├── sample_name_contigs_filtered.fa            # Contigs >1kb
│
├── BLAST Results
├── sample_name_blast_all.tsv                  # All BLAST hits
├── sample_name_viral_blast.tsv                # Viral hits only
├── sample_name_top_hits.tsv                   # Best hit per contig
│
└── blast_results/
    ├── viral_hits.txt                         # Classified viral
    ├── bacterial_hits.txt                     # Classified bacterial
    └── fungal_hits.txt                        # Classified fungal
```

### Main Diagnostic Report

**File:** `sample_name_diagnostic_report.txt`

**Key sections:**

#### 1. Mapping Statistics

```
========================================
MAPPING STATISTICS
========================================
Total Reads: 1234567
Duplicate Reads: 450000 (36.5%)
Unique Reads: 784567

Raw Mapping: 1200000 reads (97.2% of total)
Deduplicated Mapping: 750000 reads (95.6% of unique)

Interpretation (based on deduplicated mapping):
- >70%: Good mapping to reference, likely correct organism
- 30-70%: Moderate mapping, possible mixed infection or variant
- <30%: Poor mapping, likely wrong reference or contamination
```

**Interpretation:**

| Deduplicated Mapping | Interpretation |
|---------------------|----------------|
| **>90%** | Excellent - clean sample, correct reference |
| **70-90%** | Good - minor contamination or strain variation |
| **30-70%** | Moderate - check BLAST for co-infection |
| **<30%** | Poor - wrong virus or major contamination |

**Note on duplication:**
- High duplication (>50%) is **normal** for deep viral sequencing
- Doesn't indicate quality problems
- Focus on **deduplicated mapping rate**

#### 2. Assembly Statistics

```
========================================
ASSEMBLY STATISTICS
========================================
Total Contigs: 245
Contigs >1000bp: 12
```

**Interpretation:**

| Total Contigs | >1kb Contigs | Interpretation |
|---------------|--------------|----------------|
| 1-20 | 1-5 | Clean - likely single virus |
| 20-100 | 5-20 | Moderate - check BLAST results |
| 100-1000 | 20-100 | Contaminated - multiple organisms |
| >1000 | >100 | Severe contamination or host genome |

#### 3. Viral Contamination Analysis

**Coverage-based classification:**

```
CONFIRMED VIRAL CONTAMINANTS (≥80% query coverage):
========================================================

West Nile virus (2 contigs):
    k141_1      89.5%  98.2%
    k141_5      85.3%  97.8%

POTENTIAL VIRAL CONTAMINANTS (50-80% query coverage):
========================================================

Dengue virus type 2 (1 contig):
    k141_3      65.2%  85.3%
```

**Columns:**
- **Contig ID** - Assembly contig identifier
- **Query Coverage** - % of contig matching BLAST hit
- **Percent Identity** - % nucleotide identity

**Interpretation:**

| Query Coverage | Identity | Interpretation |
|----------------|----------|----------------|
| **≥80%** | ≥95% | CONFIRMED - Same species/strain |
| **≥80%** | 90-95% | CONFIRMED - Related strain |
| **≥80%** | 80-90% | CONFIRMED - Related species |
| **50-80%** | ≥95% | POTENTIAL - Partial match or chimera |
| **50-80%** | <95% | POTENTIAL - Conserved region only |
| **<50%** | Any | EXCLUDED - Non-specific hit |

**Why exclude <50% coverage?**
- Often conserved protein domains
- Could be host genome with viral-like regions
- Not sufficient evidence for contamination

### Top Hits Summary

**File:** `sample_name_top_hits.tsv`

**Format:**
```
Contig_ID  Length  Subject_ID  Percent_ID  Alignment_Length  Query_Coverage  E_value  Kingdom  Subject_Title
k141_1     12809   NC_001477   98.5%       12500            97.5%           0.0      Virus    West Nile virus, complete genome
k141_2     5234    NZ_CP012345 95.2%       4800             91.7%           0.0      Bacteria Escherichia coli strain...
```

**Best practices:**

1. **Sort by contig length** (largest first)
   - Longer contigs are more reliable
   - <1kb contigs may be assembly artifacts

2. **Focus on high coverage** (>80%)
   - Full-length hits are confident
   - Partial hits (<50%) may be non-specific

3. **Check E-values** (should be near 0)
   - E-value > 1e-10 is suspicious
   - Very high E-value = poor match

### HTML Report

**File:** `diagnostic_sample_name_presentation_ready_report.html`

**View in browser:**
```bash
firefox diagnostic_sample_name_presentation_ready_report.html
```

**Contains:**
- Visual summary with plots
- Coverage distribution
- Variant frequency histogram
- Top contamination hits
- Ready for presentations/meetings

---

## Publication Reporting

### Methods Section Language

**Example text:**

```
Quality Control and Contamination Screening

To assess sample purity and identify potential contamination, we performed
de novo assembly and BLAST-based contamination screening. Quality-controlled
reads were assembled de novo using MEGAHIT v1.2.9 (meta-sensitive preset,
minimum contig length 500 bp). Assembled contigs longer than 1 kb were
searched against a custom microbial database using BLASTn (E-value < 1e-10,
top 5 hits per contig). Contamination was assessed based on BLAST query
coverage: contigs with ≥80% coverage to non-target organisms were classified
as confirmed contamination, while 50-80% coverage were considered potential
contamination. Contigs with <50% coverage were excluded as non-specific matches.

Mapping statistics were calculated using BWA-MEM v0.7.17 and SAMtools v1.15.
Deduplicated mapping rates (percentage of unique reads mapping to reference)
were used to assess sample quality, with >70% considered acceptable for analysis.
```

### Results Section Language

**Clean sample:**

```
Contamination screening revealed [X]% of deduplicated reads mapped to the
expected virus reference (NC_001477). De novo assembly produced [X] contigs
>1 kb, all of which showed >95% identity to [Virus Name] with no evidence of
viral, bacterial, or fungal contamination (Table S1).
```

**Minor contamination:**

```
Contamination screening showed [X]% deduplicated mapping to the target virus.
De novo assembly identified [X] contigs matching [target virus] and [X] minor
contigs (<2% of total assembly) matching [contaminant species]. These minor
contaminants were removed from downstream analysis. All reported variants were
confirmed to map exclusively to the target viral genome.
```

**Co-infection:**

```
Contamination screening revealed evidence of co-infection with [Virus A] and
[Virus B]. Mapping to [Virus A] reference (NC_XXXXXX) showed [X]% of reads,
while [X]% mapped to [Virus B] (NC_YYYYYY). De novo assembly produced distinct
contig sets matching each virus (Table S2). Subsequent analyses were performed
separately for each virus.
```

### Supplementary Materials

**Table S1: Contamination Screening Results**

| Sample | Total Reads | Dedup Mapping (%) | Contigs >1kb | Confirmed Contaminants | Potential Contaminants |
|--------|-------------|-------------------|--------------|------------------------|------------------------|
| P0_rep1 | 1,234,567 | 96.8% | 5 | None | None |
| P0_rep2 | 1,456,789 | 95.3% | 6 | None | None |
| P5_rep1 | 2,345,678 | 94.1% | 8 | None | 1 (E. coli, 1.2%) |

**Figure S1: Contamination Report**
- Include HTML report screenshot
- Or export key plots from HTML

### Addressing Reviewer Concerns

**Common reviewer questions:**

**Q: "How do you know this is not contamination?"**
```
A: We performed comprehensive contamination screening using de novo assembly
and BLAST analysis (Methods, Contamination Screening). All samples showed
>90% deduplicated mapping to the expected virus with no confirmed viral,
bacterial, or fungal contaminants (Table S1).
```

**Q: "What about bacterial contamination?"**
```
A: BLAST screening of assembled contigs against a microbial database revealed
no bacterial hits with >50% query coverage in any sample. Minor hits (<1% of
reads) to Escherichia coli were observed, likely from lab environment, but
did not affect variant calling as they did not map to the viral reference.
```

**Q: "Could this be a co-infection?"**
```
A: Contamination screening explicitly tested for co-infection by BLAST
screening against a comprehensive viral database. Only [target virus] was
detected with >80% query coverage. No other viral species exceeded 50%
coverage threshold (Table S1).
```

---

## Troubleshooting

### Issue: No Contigs Assembled

**Cause:** Very low read count or all reads are duplicates

**Solution:**
```bash
# Check read count
zcat R1.fastq.gz | wc -l
# Should be >4000 lines (1000 reads minimum)

# Check fastp report
cat cleaned_seqs/sample_fastp_report.json | grep total_reads

# If <100k reads, may need more sequencing
```

### Issue: Thousands of Contigs

**Cause:** Severe contamination or host genome presence

**Solution:**
```bash
# Check top hits
head -20 sample_top_hits.tsv

# If mostly host genome:
# - Improve DNase treatment
# - Enrich for viral particles before sequencing

# If mostly bacteria:
# - Contamination during culture
# - Use antibiotics in culture medium
```

### Issue: BLAST Takes Forever

**Cause:** Remote BLAST with many contigs

**Solution:**
```bash
# Cancel job
# Set up local BLAST database (see BLAST Database Setup)
# Re-run with local database

export BLAST_DB=/path/to/local/database
./viral_diagnostic.sh R1.fq.gz R2.fq.gz NC_001477 sample 4
```

### Issue: "No significant BLAST hits found"

**Cause:** Novel virus or wrong database

**Solution:**
```bash
# Check if contigs exist
ls -lh assembly_sample/final.contigs.fa

# Try different database
# Or use NCBI remote BLAST:
# (automatically used if local database not available)

# Inspect contigs manually
cat sample_contigs_filtered.fa
```

### Issue: Mapping Stats Show 100% Duplication

**Cause:** Calculation error or actual PCR artifact

**Solution:**
```bash
# Check raw fastp output
cat cleaned_seqs/sample_fastp_report.json | grep dup_rate

# If truly 100%, likely issue with:
# - Over-amplification during library prep
# - Very low input material
# - Need to re-sequence with less PCR cycles
```

---

## Best Practices

### When to Run

1. **Always for publications** - Non-negotiable
2. **First sample of new virus** - Validate identity
3. **After protocol changes** - Verify no new contamination introduced
4. **Unexpected results** - Rule out contamination as cause
5. **Low mapping rates** - Identify what else is in sample

### Quality Standards

**Acceptable for publication:**
- Deduplicated mapping ≥70%
- No confirmed viral contaminants (≥80% coverage)
- Bacterial/fungal contaminants <5% of reads
- Top 10 contigs all match target virus

**Requires attention:**
- Deduplicated mapping 50-70%
- Potential viral contaminants (50-80% coverage)
- Bacterial contaminants 5-20%

**Reject sample:**
- Deduplicated mapping <50%
- Confirmed co-infection (unless that's the study question)
- Bacterial contamination >20%
- Most contigs don't match target virus

### Preventing Contamination

**Lab practices:**
1. Use dedicated pipettes for viral work
2. UV treat biosafety cabinets between samples
3. Use filter tips
4. Negative controls in every sequencing batch
5. Sequence controls alongside samples

**Computational:**
1. Always check for cross-contamination between samples
2. Run controls through same pipeline
3. Look for sample swaps (if multiple viruses)

---

## Advanced Topics

### Custom BLAST Databases

**Create virus-specific database:**

```bash
# Download all viruses of interest
cat virus1.fasta virus2.fasta virus3.fasta > my_viruses.fasta

makeblastdb -in my_viruses.fasta \
            -dbtype nucl \
            -parse_seqids \
            -out my_viruses_db \
            -title "Lab Virus Database"

# Use in diagnostic
export BLAST_DB=/path/to/my_viruses_db
```

**Benefits:**
- Faster BLAST (smaller database)
- More relevant hits
- Easier interpretation

### Comparing Across Samples

**Batch analysis:**

```bash
#!/bin/bash
# compare_contamination.sh

SAMPLES="sample1 sample2 sample3"

for SAMPLE in $SAMPLES; do
    ./viral_diagnostic.sh \
        ${SAMPLE}_R1.fastq.gz \
        ${SAMPLE}_R2.fastq.gz \
        NC_001477 \
        $SAMPLE \
        8
done

# Summarize results
echo "Sample,Mapping%,Contigs,Viral_Hits" > contamination_summary.csv

for SAMPLE in $SAMPLES; do
    MAPPING=$(grep "Deduplicated Mapping:" diagnostic_${SAMPLE}/${SAMPLE}_diagnostic_report.txt | awk '{print $3}')
    CONTIGS=$(grep "Contigs >1000bp:" diagnostic_${SAMPLE}/${SAMPLE}_diagnostic_report.txt | awk '{print $3}')
    VIRAL=$(wc -l < diagnostic_${SAMPLE}/${SAMPLE}_viral_blast.tsv)

    echo "$SAMPLE,$MAPPING,$CONTIGS,$VIRAL" >> contamination_summary.csv
done
```

### Integration with Main Pipeline

**Automatic contamination flagging:**

```bash
# In your analysis script
MAPPING=$(grep "Deduplicated Mapping:" diagnostic_${SAMPLE}/${SAMPLE}_diagnostic_report.txt | \
          awk '{print $3}' | sed 's/%//')

if (( $(echo "$MAPPING < 70" | bc -l) )); then
    echo "WARNING: Low mapping rate for $SAMPLE - check contamination report!"
    # Flag for manual review
fi
```

---

## Quick Reference

### Essential Commands

```bash
# Basic contamination check
./viral_diagnostic.sh R1.fq.gz R2.fq.gz NC_001477 sample_name 4

# With large files
./viral_diagnostic.sh R1.fq.gz R2.fq.gz NC_001477 sample_name 8 --extremely-large-files

# Check results
cat diagnostic_sample_name/sample_name_diagnostic_report.txt

# View HTML report
firefox diagnostic_sample_name/diagnostic_sample_name_presentation_ready_report.html
```

### Interpretation Cheat Sheet

| Metric | Good | Caution | Bad |
|--------|------|---------|-----|
| Dedup Mapping | >90% | 70-90% | <70% |
| Total Contigs | <50 | 50-100 | >100 |
| Contigs >1kb | <10 | 10-20 | >20 |
| Viral Contaminants | 0 | 1-2 minor | >2 or major |

### Output Files

| File | Purpose |
|------|---------|
| `*_diagnostic_report.txt` | Main summary - read this first |
| `*_presentation_ready_report.html` | Visual report for presentations |
| `*_top_hits.tsv` | Best BLAST hit per contig |
| `*_viral_blast.tsv` | All viral BLAST hits |
| `final.contigs.fa` | All assembled contigs |

---

## Next Steps

After contamination screening:
- **Clean samples:** Proceed with [Variant Calling Guide](VICAST_ANALYZE_GUIDE.md)
- **Contaminated samples:** Troubleshoot lab protocol or re-sequence
- **Co-infections:** Analyze each virus separately

---

**Last Updated:** 2026-02-05
**VICAST Version:** 2.2.0
