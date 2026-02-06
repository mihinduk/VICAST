# VICAST-Analyze: Variant Calling Pipeline Guide

**Module:** vicast-analyze
**Version:** 2.2.0
**Last Updated:** 2026-02-05

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Complete Workflow](#complete-workflow)
4. [Quality Control Parameters](#quality-control-parameters)
5. [Alignment and Mapping](#alignment-and-mapping)
6. [Variant Calling](#variant-calling)
7. [Filtering Strategies](#filtering-strategies)
8. [Output Interpretation](#output-interpretation)
9. [Advanced Topics](#advanced-topics)

---

## Overview

VICAST-analyze is a comprehensive variant calling pipeline optimized for **cultured virus passage studies**, designed to detect and track low-frequency variants across viral populations.

### Key Features

- ✅ **Quality control** with fastp (trimming, filtering, deduplication stats)
- ✅ **Read alignment** with bwa mem (optimized for viral genomes)
- ✅ **Low-frequency variant calling** with lofreq (down to 1% frequency)
- ✅ **Variant annotation** with SnpEff (functional effects)
- ✅ **Contamination detection** via de novo assembly + BLAST
- ✅ **Multi-tier analysis** (high/medium/low frequency variants)
- ✅ **Consensus generation** from dominant variants
- ✅ **Quasispecies support** with individual variant proteins

### Use Cases

| Use Case | Typical Parameters | Key Output |
|----------|-------------------|------------|
| **Passage tracking** | freq=0.01, depth=200 | Mutation emergence/fixation |
| **Consensus genome** | freq=0.50, depth=200 | Dominant sequence |
| **Quasispecies** | freq=0.01, quality=1000 | Population diversity |
| **SNP calling** | freq=0.20, quality=1000 | Major variants only |
| **Contamination check** | viral_diagnostic.sh | Assembly + BLAST |

---

## Pipeline Architecture

### Three-Chunk Workflow

VICAST-analyze uses a **three-chunk architecture** with two manual decision points:

```
┌──────────────────────────────────────────┐
│ CHUNK 1: QC & Diagnostics (Automated)   │
│ Script: run_vicast_analyze_qc_only.sh   │
│ Steps: 1-6                               │
│   1. Prepare reference genome            │
│   2. Calculate read statistics           │
│   3. Clean reads (fastp QC)              │
│   4. Map reads (bwa) & call variants     │
│   5. Generate depth file                 │
│   6. Run diagnostic report               │
└──────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────┐
│ 🔍 DECISION POINT 1: Quality Assessment │
│ Review: Coverage, contamination          │
│ Decide: Proceed? Reject? Re-sequence?   │
└──────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────┐
│ CHUNK 2: Annotation (Automated)         │
│ Script: run_vicast_analyze_annotate_only │
│ Steps: 7-9                               │
│   7. Filter variants                     │
│   8. Annotate with SnpEff                │
│   9. Parse to TSV                        │
└──────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────┐
│ 🔍 DECISION POINT 2: Parameter Tuning   │
│ Review: Variant distribution, quality    │
│ Decide: Adjust thresholds?               │
└──────────────────────────────────────────┘
                    ↓
┌──────────────────────────────────────────┐
│ CHUNK 3: Analysis (Manual)              │
│ Scripts: parse_snpeff_tsv.py            │
│          generate_realistic_haplotype    │
│   - Filter mutations (custom params)     │
│   - Generate consensus genomes           │
│   - Create variant proteins              │
└──────────────────────────────────────────┘
```

### Why This Architecture?

**Benefits:**
- **Quality gates** - Stop early if sample is bad
- **Parameter optimization** - Adjust based on data characteristics
- **Efficiency** - Resume from VCF (skip hours of re-processing)
- **Flexibility** - Customize filtering for each analysis goal

---

## Complete Workflow

### Prerequisites

```bash
# Activate environment
conda activate vicast_analyze

# Verify tools
which bwa samtools lofreq fastp

# Verify SnpEff database exists
java -jar $SNPEFF_JAR databases | grep NC_001477
# If not, create with VICAST-annotate first!
```

### Input Files

**Required:**
- **R1.fastq.gz** - Forward reads (Illumina paired-end)
- **R2.fastq.gz** - Reverse reads
- **Accession** - Reference genome (must be in SnpEff database)

**Optional:**
- **Threads** - Number of CPU cores (default: 4)

### Workflow Option 1: Three-Chunk (Recommended)

**Chunk 1: QC & Diagnostics**

```bash
cd /path/to/your/data

/path/to/VICAST/vicast-analyze/run_vicast_analyze_qc_only.sh \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    NC_001477 \
    8  # threads
```

**Outputs:**
```
sample_results/sample_depth.txt           # Coverage per position
diagnostic_sample/sample_diagnostic_report.txt  # Text summary
diagnostic_sample/diagnostic_sample_presentation_ready_report.html  # HTML report
diagnostic_sample/blast_results/viral_hits.txt   # Viral contamination
cleaned_seqs/variants/sample_vars.vcf     # Raw variants
```

**Decision Point 1: Review QC**

```bash
# Check coverage
awk 'NR>1 {sum+=$3; count++} END {print "Mean depth:", sum/count}' \
    sample_results/sample_depth.txt

# Check contamination
cat diagnostic_sample/sample_diagnostic_report.txt

# Check mapping rate
grep "mapped" diagnostic_sample/sample_mapping_stats.txt
```

**If QC passes, continue:**

**Chunk 2: Annotation**

```bash
/path/to/VICAST/vicast-analyze/run_vicast_analyze_annotate_only.sh \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    NC_001477
```

**Outputs:**
```
cleaned_seqs/variants/sample.snpEFF.ann.vcf   # Annotated VCF
cleaned_seqs/variants/sample.snpEFF.ann.tsv   # Annotated TSV (easier to read)
```

**Decision Point 2: Evaluate Parameters**

```bash
# Count variants
grep -v "^#" cleaned_seqs/variants/sample.snpEFF.ann.tsv | wc -l

# Check frequency distribution
grep -v "^#" cleaned_seqs/variants/sample.snpEFF.ann.tsv | \
    awk -F'\t' '{print $8}' | grep -oP 'AF=[0-9.]+' | cut -d= -f2 | sort -n
```

**Chunk 3: Custom Analysis**

```bash
# Filter mutations (quasispecies detection)
python /path/to/VICAST/vicast-analyze/parse_snpeff_tsv.py \
    cleaned_seqs/variants/sample.snpEFF.ann.tsv \
    sample_results/sample_filtered_mutations.tsv \
    --quality 1000 \
    --depth 200 \
    --freq 0.01

# Generate consensus (majority variants)
python /path/to/VICAST/vicast-analyze/generate_realistic_haplotype_consensus.py \
    --vcf sample_results/sample_filtered_mutations.tsv \
    --reference cleaned_seqs/NC_001477.fasta \
    --accession NC_001477 \
    --quality 1000 \
    --depth 200 \
    --freq 0.50 \
    --output-prefix sample_results/sample_consensus
```

### Workflow Option 2: Full Pipeline (All at Once)

**For experienced users with known parameters:**

```bash
/path/to/VICAST/vicast-analyze/run_vicast_analyze_full.sh \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    NC_001477 \
    8
```

**Runs all steps 1-9 without stopping.**

**Use when:**
- You trust your data quality
- Parameters are well-established for your virus
- Running in batch mode

---

## Quality Control Parameters

### Step 3: fastp QC

**Default settings (optimized for viral sequencing):**

```bash
fastp \
    --qualified_quality_phred 15 \      # Q15 minimum
    --unqualified_percent_limit 40 \    # Max 40% low-quality bases
    --length_required 36 \              # Minimum read length 36bp
    --thread 4
```

**Adjust for different sequencing quality:**

| Sequencing Quality | --qualified_quality_phred | --length_required |
|--------------------|--------------------------|-------------------|
| High (NovaSeq) | 20 | 50 |
| Standard (MiSeq) | 15 (default) | 36 (default) |
| Low (older platforms) | 10 | 30 |

**Customize in viral_pipeline.py:**

```python
# Line ~200
fastp_cmd = [
    "fastp",
    "-i", r1, "-I", r2,
    "-o", r1_qc, "-O", r2_qc,
    "--qualified_quality_phred", "20",  # Change here
    "--length_required", "50",          # Change here
]
```

### Quality Metrics

**Good QC indicators:**
- **Total reads:** >1M (for 10kb genome, >100x coverage)
- **% passing filter:** >80%
- **Duplication rate:** Variable (high is OK for deep sequencing)
- **GC content:** Should match expected viral genome

**Review fastp report:**

```bash
# Open HTML report in browser
firefox cleaned_seqs/sample_fastp_report.html

# Or check JSON
grep "total_reads" cleaned_seqs/sample_fastp_report.json
```

---

## Alignment and Mapping

### Step 4: bwa mem Alignment

**Default bwa mem parameters:**

```bash
bwa mem -t 4 reference.fasta R1.qc.fastq.gz R2.qc.fastq.gz
```

**Why bwa mem for viruses?**
- Fast and accurate
- Handles high coverage well
- Good for short genomes (<50kb)
- Tolerates some mismatches (good for divergent strains)

### Mapping Quality Assessment

**Key metrics:**

```bash
# View mapping stats
samtools flagstat sample_aligned.bam
```

**Output:**
```
1234567 + 0 in total (QC-passed reads + QC-failed reads)
1200000 + 0 mapped (97.20% : N/A)
1234567 + 0 paired in sequencing
617283 + 0 read1
617284 + 0 read2
```

**Interpretation:**

| Mapping Rate | Interpretation | Action |
|--------------|----------------|--------|
| >90% | Excellent - correct reference | Proceed |
| 70-90% | Good - possible strain variation | Check variants |
| 30-70% | Moderate - wrong strain or contamination | Check BLAST results |
| <30% | Poor - wrong reference or no virus | Use different reference |

### Depth of Coverage

**Check depth distribution:**

```bash
# Summary statistics
awk 'NR>1 {sum+=$3; count++; if($3<50) low++}
     END {print "Mean:", sum/count, "Low (<50x):", low}' \
     sample_results/sample_depth.txt
```

**Recommended coverage:**

| Analysis Goal | Minimum Depth | Ideal Depth |
|---------------|---------------|-------------|
| Consensus genome | 20x | 100x |
| Major SNPs | 50x | 200x |
| Quasispecies (1%) | 200x | 1000x |
| Rare variants (0.5%) | 500x | 2000x |

---

## Variant Calling

### Step 4 (continued): lofreq Variant Calling

**Why lofreq?**
- Designed for **low-frequency variants**
- No need to specify ploidy (unlike GATK, freebayes)
- Calculates quality scores accounting for sequencing errors
- Good for viral quasispecies

**Default lofreq parameters:**

```bash
lofreq call \
    --call-indels \         # Include insertions/deletions
    -f reference.fasta \
    -o variants.vcf \
    sample_aligned.bam
```

**Key lofreq features:**

1. **No frequency cutoff** - Calls all variants with statistical support
2. **Quality-aware** - Considers base quality, mapping quality
3. **Strand bias filter** - Removes artifacts
4. **Indel calling** - Detects insertions and deletions

### Variant Quality Scores

**VCF QUAL field interpretation:**

| QUAL Score | Meaning | Phred Error Probability |
|------------|---------|------------------------|
| 10 | Low confidence | 10% chance of error |
| 100 | Moderate | 0.00001% error |
| 1000 | High | ~0% error |
| 10000+ | Very high | Essentially certain |

**Default filtering: QUAL >= 1000** (very high confidence)

---

## Filtering Strategies

### Two-Stage Filtering

**Stage 1: Parse variants (quasispecies detection)**

```bash
python parse_snpeff_tsv.py \
    input.snpEFF.ann.tsv \
    output_filtered.tsv \
    --quality 1000 \
    --depth 200 \
    --freq 0.01
```

**Stage 2: Consensus generation (majority rule)**

```bash
python generate_realistic_haplotype_consensus.py \
    --vcf output_filtered.tsv \
    --reference ref.fasta \
    --accession NC_001477 \
    --quality 1000 \
    --depth 200 \
    --freq 0.50 \
    --output-prefix sample_consensus
```

### Parameter Selection Guide

#### Quality Score (--quality)

**Recommended values:**

| Sample Quality | --quality | Rationale |
|----------------|-----------|-----------|
| High (>500x, clean) | 1000-5000 | Strict - only confident calls |
| Medium (200-500x) | 500-1000 | Moderate |
| Low (50-200x) | 100-500 | Relaxed - accept more calls |

**Example:**
```bash
--quality 1000  # Default: High confidence only
--quality 500   # Medium: More permissive
--quality 100   # Low: Very permissive (use with caution)
```

#### Depth Cutoff (--depth)

**Recommended values:**

| Analysis Goal | --depth | Rationale |
|---------------|---------|-----------|
| Consensus only | 20-50 | Basic coverage |
| Major variants | 100-200 | Good support |
| Quasispecies (1%) | 200-500 | At least 2-5 reads per 1% variant |
| Rare variants (0.5%) | 500-1000 | High depth needed |

**Calculate minimum depth for frequency:**
```
min_depth = 100 / min_frequency
# For 1% variants: 100 / 1 = 100x minimum
# For 0.5%: 100 / 0.5 = 200x minimum
```

#### Frequency Cutoff (--freq)

**Stage 1 (parse_snpeff_tsv.py): Broad detection**

| Goal | --freq | Description |
|------|--------|-------------|
| Quasispecies | 0.01 (1%) | Detect minority variants |
| Major variants | 0.05 (5%) | Moderate frequency |
| Dominant only | 0.20 (20%) | High frequency only |

**Stage 2 (consensus): Majority rule**

| Goal | --freq | Description |
|------|--------|-------------|
| Consensus genome | 0.50 (50%) | Majority allele |
| Near-fixation | 0.80 (80%) | Almost fixed |
| Fixed only | 0.95 (95%) | Essentially unanimous |

### Example Parameter Sets

**Conservative (publication-quality):**
```bash
# Stage 1: Detect quasispecies
--quality 1000 --depth 200 --freq 0.01

# Stage 2: Consensus
--quality 1000 --depth 200 --freq 0.50
```

**Exploratory (passage tracking):**
```bash
# Stage 1: Sensitive detection
--quality 500 --depth 100 --freq 0.01

# Stage 2: Dominant variants
--quality 1000 --depth 200 --freq 0.50
```

**Strict (low coverage sample):**
```bash
# Stage 1: Require more reads per variant
--quality 1000 --depth 100 --freq 0.05

# Stage 2: Majority rule
--quality 1000 --depth 100 --freq 0.50
```

---

### Manual Checkpoint 2: Variant Review and Filtering

After automated filtering (Step 8A), you should manually review variants before final annotation. **You have two options:**

---

#### 🔀 Choose Your Pathway

| Criterion | Pathway A: Frequency-Based | Pathway B: BAM Co-Occurrence |
|-----------|---------------------------|------------------------------|
| **Time required** | 30-60 min | 2-5 hours (1-3h compute + 1-2h review) |
| **Expertise needed** | Basic virology | Basic virology + read interpretation |
| **Best for** | Routine analysis, passage tracking | Publications, novel findings |
| **Evidence type** | Frequency + quality scores | Read-level linkage proof |
| **Variant confidence** | Conservative thresholds | Evidence-based decisions |
| **False positive rate** | Higher (conservative removal) | Lower (validated variants) |
| **Biological claims** | "Variants detected at X%" | "Variants co-occur on same reads" |
| **Reviewer questions** | "How do you know it's real?" | Direct evidence provided |

**Decision Guide:**

**Use Pathway A (Frequency-Based) when:**
- ✅ Routine passage studies with well-characterized viruses
- ✅ Time constraints (need results quickly)
- ✅ Previous samples showed clean, straightforward patterns
- ✅ Variants are clearly high-frequency (>50%) or clearly low (<5%)
- ✅ Internal exploratory work (not for publication yet)

**Use Pathway B (BAM Co-Occurrence) when:**
- ✅ Preparing for publication submission
- ✅ Novel or unexpected findings (need to prove they're real)
- ✅ Multiple mid-frequency variants (20-50%) - need to know if linked
- ✅ Co-infection suspected - need to distinguish haplotypes
- ✅ Borderline variants from Pathway A need validation
- ✅ Reviewers likely to ask "How do you know these are linked?"

**Can you switch pathways mid-analysis?**
- ✅ YES! Start with Pathway A for quick review
- If you find borderline/uncertain cases → Run Pathway B for those specific regions
- BAM analysis can be run anytime on the same .bam and .vcf files

---

#### Pathway A: Frequency-Based Review (Standard)

**Use when:**
- Routine passage tracking
- Time constraints
- Well-characterized viral systems
- Previous samples showed clean variant patterns

**Decision criteria based on:**
- Variant frequency
- Quality scores
- Coverage depth
- Genomic context
- Biological plausibility

##### Standard Review Process

```bash
# Review filtered VCF
less cleaned_seqs/variants/sample_vars.filt.vcf

# Or export to TSV for easier viewing
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\t%INFO/AF\n' \
    cleaned_seqs/variants/sample_vars.filt.vcf > variants_summary.tsv
```

##### Frequency-Based Decision Rules

**Rule 1: High-frequency variants (≥50%)**

| Frequency | Quality | Decision | Rationale |
|-----------|---------|----------|-----------|
| ≥80% | ≥1000 | ✅ KEEP | Dominant variant, very reliable |
| 50-80% | ≥1000 | ✅ KEEP | Major variant, reliable |
| 50-80% | <1000 | ⚠️ CHECK | Review genomic context |

**Action for low quality:**
```bash
# Check alignment at position
samtools view cleaned_seqs/mapping/sample.bam NC_001477.1:POS-50-POS+50 | less

# Look for:
# - Homopolymer runs (AAAA, TTTT)
# - Low-complexity regions
# - Alignment artifacts
```

---

**Rule 2: Mid-frequency variants (20-50%)**

**Consider:**
- Is this near other mid-frequency variants? (possible mixed population)
- Is this in a variable region? (biological expectation)
- Is the quality score high? (>1500)

| Context | Decision | Action |
|---------|----------|--------|
| Multiple nearby mid-freq variants | ✅ KEEP | Possible co-infection or mixed population |
| Known variable region (e.g., Spike RBD) | ✅ KEEP | Biologically expected |
| Isolated + high quality | ✅ KEEP | Mark as uncertain, report conservatively |
| Isolated + low quality | ❌ REMOVE | Likely artifact |
| In homopolymer region | ❌ REMOVE | Indel artifact |

---

**Rule 3: Low-frequency variants (5-20%)**

**Conservative approach (recommended without BAM evidence):**

| Frequency | Depth | Quality | Decision |
|-----------|-------|---------|----------|
| 10-20% | ≥500X | ≥2000 | ✅ KEEP |
| 10-20% | 200-500X | ≥1500 | ⚠️ CHECK | Review context |
| 10-20% | <200X | Any | ❌ REMOVE | Insufficient evidence |
| 5-10% | ≥500X | ≥2000 | ⚠️ UNCERTAIN | Report with caveats |
| 5-10% | <500X | Any | ❌ REMOVE | Too unreliable |
| <5% | Any | Any | ❌ REMOVE | Below detection threshold |

**Rationale:** Without read-level evidence, we cannot prove these variants are real vs. artifacts

---

**Rule 4: Known hotspots and biological context**

✅ **KEEP even if borderline, when:**
- In known variable regions (literature support)
- In epitope regions (immune pressure)
- In drug target sites (resistance)
- Consistent with passage history

❌ **REMOVE even if frequency is OK, when:**
- In conserved essential regions (unlikely to be real)
- Premature stop codons in essential genes (lethal)
- Multiple variants at same codon (likely artifact)

---

##### Simplified Review Template

Use this for frequency-based filtering:

```bash
# Create review spreadsheet
cat > variant_review_frequency_based.tsv << 'EOF'
Position	Frequency	Quality	Depth	Gene	Effect	Decision	Confidence	Rationale	Notes
# Fill in for each variant that needs review
# Decision: KEEP/REMOVE/UNCERTAIN
# Confidence: HIGH/MEDIUM/LOW
# Rationale: Why you made this decision
EOF
```

**Example entries:**

```tsv
Position	Freq	Quality	Depth	Gene	Effect	Decision	Confidence	Rationale	Notes
2450	0.92	3500	8000	E	D614G	KEEP	HIGH	High freq, high quality, spike RBD	Known variable site
5000	0.45	1200	3000	NS3	L234F	REMOVE	LOW	Mid-freq orphan, homopolymer region	AAAAAA run
7500	0.18	2200	600	NS5	A456V	UNCERTAIN	MEDIUM	Borderline freq, depth too low	Report with caveat
9000	0.08	1800	450	NS5	K789R	REMOVE	LOW	Below 10%, insufficient depth	Below threshold
```

---

##### Quality Control Checklist (Frequency-Based)

```
┌────────────────────────────────────────────────────┐
│ VARIANT REVIEW CHECKLIST - Frequency-Based        │
├────────────────────────────────────────────────────┤
│ Total variants after automated filtering: ___     │
│                                                    │
│ HIGH-FREQUENCY (≥50%):                             │
│   - Total: ___                                     │
│   - KEEP: ___                                      │
│   - REMOVE: ___                                    │
│   - UNCERTAIN: ___                                 │
│                                                    │
│ MID-FREQUENCY (20-50%):                            │
│   - Total: ___                                     │
│   - KEEP: ___                                      │
│   - REMOVE: ___                                    │
│   - Flagged for biological review: ___            │
│                                                    │
│ LOW-FREQUENCY (5-20%):                             │
│   - Total: ___                                     │
│   - KEEP (high confidence): ___                    │
│   - UNCERTAIN (report with caveats): ___          │
│   - REMOVE: ___                                    │
│                                                    │
│ VERY LOW (<5%):                                    │
│   - All removed: YES/NO                            │
│   - If NO, explain: ___                            │
│                                                    │
│ DECISION:                                          │
│   [ ] Proceed to annotation with current set      │
│   [ ] Need BAM co-occurrence for borderline cases │
│   [ ] Consult literature for biological context   │
└────────────────────────────────────────────────────┘
```

---

##### Conservative Filtering Recommendation

**Without BAM evidence, use conservative thresholds:**

```bash
# More stringent than with BAM evidence
# Err on the side of removing uncertain variants

# High-frequency only (safest)
--quality 1500 --depth 200 --freq 0.50

# Include moderate-frequency (standard)
--quality 1000 --depth 200 --freq 0.20

# Include low-frequency (requires high depth)
--quality 2000 --depth 500 --freq 0.10
```

**Publication language (frequency-based only):**
> "Variants were filtered based on quality score (≥1000), depth (≥200X), and frequency (≥20%). Manual review removed variants in low-complexity regions and those with biological implausibility. Without read-level linkage evidence, we applied conservative criteria to minimize false positives."

---

### Pathway B: Evidence-Based Review with BAM Co-occurrence (Enhanced)

**Use when:**
- Publication-quality analyses
- Novel or unexpected findings
- Complex samples (co-infections, mixed populations)
- Need to validate variant linkage for biological claims
- Borderline variants from Pathway A need resolution

#### Step 8B: BAM Co-occurrence Analysis (Optional but Recommended)

**Enhanced Variant Review with Read-Level Evidence**

After automated filtering (Step 8A), you can generate read-level co-occurrence data to make evidence-based decisions instead of relying solely on frequency thresholds.

#### When to Use This Step

✅ **Highly recommended for:**
- Publication-quality analyses
- Studies with unexpected findings
- Samples with complex haplotype structure
- When you need to distinguish real variants from artifacts

⚠️ **Optional for:**
- Routine passage tracking
- Well-characterized systems
- Time-constrained exploratory work

#### Running Co-occurrence Analysis

```bash
# After Step 8A: Automated filtering completes
# Generate read-level co-occurrence evidence

python ../scripts/check_read_cooccurrence.py \
    --bam cleaned_seqs/mapping/sample.lofreq.realign.bam \
    --vcf cleaned_seqs/variants/sample_vars.filt.vcf \
    --output cleaned_seqs/variants/sample_cooccurrence.tsv \
    --max-distance 500
```

**Parameters:**
- `--bam`: Aligned BAM file (from Step 4)
- `--vcf`: Filtered VCF file (from Step 8A)
- `--output`: Output TSV file with co-occurrence metrics
- `--max-distance`: Maximum bp between variants to check (default: 500)

**Runtime:** ~1-3 hours for typical viral genome with 100-300 variants

#### Output Format

The co-occurrence TSV contains:

```tsv
variant1_pos  variant2_pos  distance  cooccurrence_rate  informative_reads  linkage_proven
22673         22674         1         0.9974             69616              True
22673         22679         6         0.9996             68697              True
5000          5100          100       0.12               42                 True
```

**Key columns:**
- **cooccurrence_rate**: Fraction of reads carrying BOTH alternative alleles
- **informative_reads**: Number of reads spanning both positions
- **linkage_proven**: Whether we have sufficient evidence (>10 reads)

#### Manual Review with Co-occurrence Data

**MANUAL CHECKPOINT 2: Enhanced Variant Filtering**

Use the co-occurrence data to make evidence-based decisions:

**Decision Rule 1: High-frequency clusters (≥80%)**
- Look for variants with >95% co-occurrence
- These are on the same haplotype → KEEP ALL

**Decision Rule 2: Mid-frequency orphans (30-70%)**
- Check if variant co-occurs with any cluster
- No linkage found → FLAG for manual inspection (possible artifact)

**Decision Rule 3: Low-frequency variants (<20%)**
- Check if linked (>80%) to major haplotype → KEEP (emerging variant)
- Isolated singleton with <10% co-occurrence → REMOVE (artifact)

**Decision Rule 4: Mutually exclusive variants**
- Two variants with 0% co-occurrence but each links to different clusters
- These represent distinct haplotypes → KEEP BOTH

#### Review Template

Download and use the variant review checklist:

```bash
# Copy template to your working directory
cp ../docs/user_guides/variant_review_template.tsv variant_review.tsv

# Edit in Excel or text editor
# Fill in decisions based on co-occurrence patterns
```

**Template columns:**
- Position, Frequency, Variant Type
- Co-occurrence Pattern (from TSV)
- Decision (KEEP/REMOVE/FLAG/UNCERTAIN)
- Confidence (HIGH/MEDIUM/LOW/NONE)
- Notes (biological context, reviewer comments)

#### Example Review Session

```bash
# 1. Load both files for comparison
less cleaned_seqs/variants/sample_vars.filt.vcf
less cleaned_seqs/variants/sample_cooccurrence.tsv

# 2. Identify high-frequency clusters
awk -F'\t' '$12 >= 0.95' sample_cooccurrence.tsv > high_linkage.tsv

# 3. Flag suspicious singletons
# (Variants with moderate frequency but no high linkage)

# 4. Apply decisions to create final VCF
# (Manual editing or scripted filtering)

# 5. Document review in checklist
vim variant_review.tsv
```

#### Integration with Workflow Scripts

**Option A: Run manually between steps**

```bash
# Run QC and get filtered VCF
./run_vicast_analyze_qc_only.sh sample_R1.fq sample_R2.fq NC_001477.1 4

# Generate co-occurrence data
python ../scripts/check_read_cooccurrence.py \
    --bam cleaned_seqs/mapping/sample.lofreq.realign.bam \
    --vcf cleaned_seqs/variants/sample_vars.filt.vcf \
    --output cleaned_seqs/variants/sample_cooccurrence.tsv

# MANUAL REVIEW HERE
# Review co-occurrence patterns and make filtering decisions

# Continue with annotation
./run_vicast_analyze_annotate_only.sh sample_R1.fq sample_R2.fq NC_001477.1
```

**Option B: Modify workflow script to include automatic pause**

Add to `run_vicast_analyze_annotate_only.sh` after filtering:

```bash
# Step 8B: BAM co-occurrence (optional)
if [ -f "${PIPELINE_BASE}/scripts/check_read_cooccurrence.py" ]; then
    echo "Generating co-occurrence evidence..."
    python "${PIPELINE_BASE}/scripts/check_read_cooccurrence.py" \
        --bam "${OUTPUT_DIR}/mapping/${SAMPLE}.lofreq.realign.bam" \
        --vcf "${OUTPUT_DIR}/variants/${SAMPLE}_vars.filt.vcf" \
        --output "${OUTPUT_DIR}/variants/${SAMPLE}_cooccurrence.tsv"

    echo ""
    echo "=========================================="
    echo "MANUAL CHECKPOINT 2: Variant Review"
    echo "=========================================="
    echo "Co-occurrence results: ${OUTPUT_DIR}/variants/${SAMPLE}_cooccurrence.tsv"
    echo "Filtered VCF: ${OUTPUT_DIR}/variants/${SAMPLE}_vars.filt.vcf"
    echo ""
    echo "Review variants using co-occurrence evidence before proceeding."
    echo "See: docs/user_guides/BAM_COOCCURRENCE_GUIDE.md"
    echo ""
    echo "Press ENTER to continue to SnpEff annotation..."
    echo "Or press Ctrl+C to stop and review manually"
    read -r
fi
```

#### Benefits for Publication

Including BAM co-occurrence analysis strengthens your publication by:

✅ **Validation**: Direct read-level evidence supports variant calls
✅ **Rigor**: Demonstrates thorough quality control beyond frequency thresholds
✅ **Transparency**: Provides reviewers with linkage evidence
✅ **Novelty**: Distinguishes your analysis from frequency-only approaches

**Example Methods text:**
> "Variant co-occurrence was assessed at the read level using pysam (v0.22). For variant pairs within 500 bp, reads spanning both positions were examined to determine linkage. Variants showing >95% co-occurrence were considered to be on the same haplotype. This read-level evidence informed manual filtering decisions (see Supplementary Table X)."

#### See Also

- **Detailed decision rules**: [BAM Co-Occurrence Guide](BAM_COOCCURRENCE_GUIDE.md#using-co-occurrence-data-for-manual-variant-review)
- **Template spreadsheet**: `docs/user_guides/variant_review_template.tsv`
- **Example analysis**: BAM test results summary

---

## Output Interpretation

### Key Output Files

```
working_directory/
├── cleaned_seqs/
│   ├── sample_R1.qc.fastq.gz          # QC'd reads
│   ├── sample_R2.qc.fastq.gz
│   ├── sample_fastp_report.html       # QC report
│   ├── sample_aligned.bam             # Aligned reads
│   ├── NC_001477.fasta                # Reference genome
│   └── variants/
│       ├── sample_vars.vcf            # Raw variants (lofreq)
│       ├── sample.snpEFF.ann.vcf      # Annotated VCF
│       └── sample.snpEFF.ann.tsv      # Annotated TSV (easier to read)
│
├── sample_results/
│   ├── sample_depth.txt               # Per-position depth
│   ├── sample_filtered_mutations.tsv  # Filtered variants
│   ├── sample_consensus_consensus.fasta    # Consensus genome
│   ├── sample_consensus_proteins.fasta     # Variant proteins
│   └── sample_consensus_summary_report.txt # Summary
│
└── diagnostic_sample/
    ├── sample_diagnostic_report.txt        # QC summary
    ├── diagnostic_sample_presentation_ready_report.html  # HTML
    ├── assembly_sample/final.contigs.fa    # De novo assembly
    └── blast_results/
        ├── viral_hits.txt                  # Viral contamination
        ├── bacterial_hits.txt              # Bacterial contamination
        └── fungal_hits.txt                 # Fungal contamination
```

### Understanding TSV Output

**Annotated TSV columns:**

| Column | Description | Example |
|--------|-------------|---------|
| CHROM | Chromosome/Reference ID | NC_001477 |
| POS | Position in genome | 1523 |
| REF | Reference allele | A |
| ALT | Alternate allele | G |
| QUAL | Quality score | 15234.0 |
| Total_Depth | Coverage at position | 856 |
| Allele_Frequency | AF (0-1) | 0.12 (12%) |
| EFFECT | Variant effect | missense_variant |
| PUTATIVE_IMPACT | Impact level | MODERATE |
| GENE_NAME | Gene name | Protein_E |
| HGVSc | Nucleotide change | c.1523A>G |
| HGVSp | Protein change | p.Lys508Arg |

**Filter in Excel/R:**

```r
# R example
library(tidyverse)

variants <- read_tsv("sample.snpEFF.ann.tsv", comment = "##")

# High-impact variants only
high_impact <- variants %>%
    filter(PUTATIVE_IMPACT == "HIGH")

# Missense mutations in envelope protein
env_missense <- variants %>%
    filter(GENE_NAME == "Protein_E",
           EFFECT == "missense_variant")
```

### Consensus Genome Output

**File: sample_consensus_consensus.fasta**

```fasta
>sample_filtered_consensus Consensus with 15 mutations (Q>=1000, D>=200, F>=0.50)
ATGAAGAACGTTCGCGTGTGGAGTATCGCCTTGTCGCTGATCATCGGGAGC...
```

**Use for:**
- Phylogenetic analysis
- Submission to GenBank
- Comparison across passages

### Variant Proteins Output

**File: sample_consensus_proteins.fasta**

**Individual proteins for each mutation:**

```fasta
>sample_ProteinE ProteinE protein [L287I] (from 862T>A)
MKTLILGAVILGVATAAQITAGIALHQSFSISKR...

>sample_ProteinE ProteinE protein [K310R] (from 929A>G)
MKTLILGAVILGVATAAQITAGIALHQSFSISKR...

>sample_ProteinNS5 ProteinNS5 protein
MDPWRKGERNGSMKLTYKCDH... (reference, no mutations)
```

**Use for:**
- Protein structure prediction (AlphaFold)
- Functional analysis of each variant
- Understanding quasispecies diversity

### Summary Report

**File: sample_consensus_summary_report.txt**

```
================================================================================
INDIVIDUAL VARIANT PROTEIN GENERATION SUMMARY REPORT
================================================================================

Total mutations passing filters: 15
Multi-allelic sites found: 2
Unique protein sequences: 18

--------------------------------------------------------------------------------
MUTATIONS PER PROTEIN
--------------------------------------------------------------------------------

ProteinE:
  - L287I: 862T>A (missense_variant, 12.50%)
  - K310R: 929A>G (missense_variant, 8.30%)

ProteinM:
  - V194I: 580G>A (missense_variant, 15.20%)
```

---

## Advanced Topics

### Large Files (High Coverage)

**For very deep sequencing (>10GB FASTQ):**

```bash
# Use extreme memory flag
./run_vicast_analyze_full.sh \
    sample_R1.fastq.gz \
    sample_R2.fastq.gz \
    NC_001477 \
    16 \
    --extremely-large-files
```

**Adjusts:**
- MEGAHIT assembly memory (90% instead of 70%)
- Samtools memory limits
- Tmp file locations

### Batch Processing

**Process multiple samples:**

```bash
#!/bin/bash
# batch_process.sh

SAMPLES="sample1 sample2 sample3 sample4"
REFERENCE="NC_001477"
THREADS=8

for SAMPLE in $SAMPLES; do
    echo "Processing $SAMPLE..."

    ./run_vicast_analyze_full.sh \
        ${SAMPLE}_R1.fastq.gz \
        ${SAMPLE}_R2.fastq.gz \
        $REFERENCE \
        $THREADS

    # Generate consensus
    python generate_realistic_haplotype_consensus.py \
        --vcf ${SAMPLE}_results/${SAMPLE}_filtered_mutations.tsv \
        --reference cleaned_seqs/${REFERENCE}.fasta \
        --accession $REFERENCE \
        --output-prefix ${SAMPLE}_results/${SAMPLE}_consensus
done
```

### HPC/SLURM Submission

**Submit as job array:**

```bash
#!/bin/bash
#SBATCH --array=1-10
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=4:00:00

# Sample list file: samples.txt
# sample1_R1.fastq.gz sample1_R2.fastq.gz
# sample2_R1.fastq.gz sample2_R2.fastq.gz
# ...

# Read sample from list
R1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt | awk '{print $1}')
R2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt | awk '{print $2}')

# Activate environment
source activate vicast_analyze

# Run pipeline
/path/to/VICAST/vicast-analyze/run_vicast_analyze_full.sh \
    $R1 $R2 NC_001477 8
```

**Submit:**
```bash
sbatch batch_job.slurm
```

### Passage Study Analysis

**Track variant frequencies across passages:**

```bash
# Run pipeline for each passage
for passage in P0 P1 P2 P3 P4 P5; do
    ./run_vicast_analyze_full.sh \
        ${passage}_R1.fastq.gz \
        ${passage}_R2.fastq.gz \
        NC_001477 8
done

# Compare filtered mutations
# Custom R/Python script to track frequency changes
```

**Example R analysis:**

```r
library(tidyverse)

# Load all passage variants
passages <- c("P0", "P1", "P2", "P3", "P4", "P5")

data_list <- lapply(passages, function(p) {
    read_tsv(paste0(p, "_results/", p, "_filtered_mutations.tsv")) %>%
        mutate(Passage = p)
})

all_data <- bind_rows(data_list)

# Track specific mutation
mutation_track <- all_data %>%
    filter(POS == 1523, REF == "A", ALT == "G") %>%
    select(Passage, Allele_Frequency)

# Plot emergence
ggplot(mutation_track, aes(x = Passage, y = Allele_Frequency, group = 1)) +
    geom_line() +
    geom_point() +
    labs(title = "Mutation A1523G Frequency Across Passages",
         y = "Allele Frequency")
```

### Custom Filtering

**Filter by specific criteria:**

```bash
# Only non-synonymous mutations
grep -v "synonymous_variant" sample_filtered_mutations.tsv > nonsynonymous.tsv

# Only HIGH impact variants
awk -F'\t' '$11 == "HIGH"' sample.snpEFF.ann.tsv > high_impact.tsv

# Specific gene
awk -F'\t' '$12 == "Protein_E"' sample.snpEFF.ann.tsv > protein_e.tsv
```

---

## Troubleshooting

### Issue: Low Mapping Rate

**Cause:** Wrong reference or contaminated sample

**Solution:**
```bash
# Check diagnostic BLAST results
cat diagnostic_sample/blast_results/viral_hits.txt

# Try different reference if strain mismatch
# Or use contamination screening guide
```

### Issue: No Variants Called

**Cause:** Perfect match to reference OR very low coverage

**Solution:**
```bash
# Check depth
cat sample_results/sample_depth.txt

# Check raw VCF
cat cleaned_seqs/variants/sample_vars.vcf

# Lower quality thresholds
--quality 100 --depth 50
```

### Issue: Too Many Variants (>1000)

**Cause:** Wrong reference, sequencing errors, or true diversity

**Solution:**
```bash
# Check if using correct reference accession
# Inspect variant types (SNPs vs indels)
grep -v "^#" sample.snpEFF.ann.tsv | cut -f8 | grep -o "TYPE=[^;]*"

# Increase quality threshold
--quality 5000 --freq 0.05
```

### Issue: SnpEff Annotation Fails

**Cause:** Chromosome name mismatch

**Solution:**
```bash
# Check VCF chromosome name
grep -v "^#" sample_vars.vcf | cut -f1 | uniq

# Check SnpEff database chromosome name
java -jar $SNPEFF_JAR dump NC_001477 | head

# Must match exactly!
```

---

## Best Practices

### For Publication-Quality Analysis

1. **Always run contamination screening**
   - Use `viral_diagnostic.sh` for every sample
   - Report in methods: "Contamination was assessed via de novo assembly and BLAST"

2. **Document parameter choices**
   - State quality/depth/frequency thresholds in methods
   - Justify choices based on coverage depth

3. **Validate key findings**
   - Manually inspect BAM files for important variants
   - Use Sanger sequencing for critical mutations
   - Check variant frequencies across biological replicates

4. **Report appropriately**
   - "Variant frequencies ≥1% with Q≥1000, depth≥200x"
   - "Consensus generated from variants ≥50% frequency"
   - "Contamination screening showed >95% reads mapped to expected virus"

---

## Quick Reference

### Essential Commands

```bash
# Full pipeline (all steps)
bash run_vicast_analyze_full.sh R1.fq.gz R2.fq.gz NC_001477 8

# Three-chunk workflow
bash run_vicast_analyze_qc_only.sh R1.fq.gz R2.fq.gz NC_001477 8
bash run_vicast_analyze_annotate_only.sh R1.fq.gz R2.fq.gz NC_001477

# Parse variants
python parse_snpeff_tsv.py input.tsv output.tsv --quality 1000 --depth 200 --freq 0.01

# Generate consensus
python generate_realistic_haplotype_consensus.py \
    --vcf filtered.tsv --reference ref.fa --accession ACC \
    --freq 0.50 --output-prefix consensus
```

### Decision Matrix

| Coverage | Mapping | Action |
|----------|---------|--------|
| >500x | >90% | Excellent - use default parameters |
| 200-500x | >80% | Good - standard analysis |
| 50-200x | >70% | Marginal - increase freq cutoff to 0.05 |
| <50x | <70% | Poor - reject or re-sequence |

---

## Next Steps

- **[Contamination Screening Guide](CONTAMINATION_SCREENING_GUIDE.md)** - Essential QC
- **[Haplotype Consensus Guide](HAPLOTYPE_CONSENSUS_GUIDE.md)** - Advanced consensus generation
- Process passage experiments
- Compare consensus genomes phylogenetically

---

**Last Updated:** 2026-02-05
**VICAST Version:** 2.2.0
