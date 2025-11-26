# VICAST Pipeline: Complete Workflow Architecture

## Overview: Three-Chunk Workflow with Manual Decision Points

The VICAST pipeline is structured as **3 automated chunks separated by 2 manual curation/decision points**:

```
┌─────────────────────────────────────────────────────────────────┐
│ CHUNK 1: QC & Diagnostics (Automated)                          │
│ Script: run_vicast_analyze_qc_only.sh                          │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 🔍 MANUAL DECISION POINT 1: Quality Assessment                 │
│ Evaluate: Depth coverage, contamination (bacteria/fungi/virus) │
│ Decision: Proceed with annotation? Reject sample? Re-sequence? │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ CHUNK 2: Variant Annotation (Automated)                        │
│ Script: run_vicast_analyze_annotate_only.sh                    │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ 🔍 MANUAL DECISION POINT 2: Parameter Optimization             │
│ Evaluate: Variant quality, depth distribution, frequency       │
│ Decision: Use defaults or adjust --quality/--depth/--freq?     │
└─────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────┐
│ CHUNK 3: Mutation Analysis & Haplotype Reconstruction          │
│ Scripts:                                                        │
│   1. parse_snpeff_tsv.py (Parse & Filter Mutations)            │
│   2. generate_realistic_haplotype_consensus.py                 │
│      (Viral Haplotype Reconstruction & Quasispecies Analysis)  │
└─────────────────────────────────────────────────────────────────┘
```

---

## CHUNK 1: QC & Diagnostics (Automated)

### Command
```bash
./run_vicast_analyze_qc_only.sh <R1.fastq.gz> <R2.fastq.gz> <ACCESSION> [threads]
```

### What It Does (Steps 1-6)
1. **Prepare reference genome** - Download from NCBI if needed
2. **Calculate read statistics** - Total reads, quality metrics
3. **Clean reads** - fastp QC (trimming, filtering)
4. **Map reads** - bwa mem alignment, samtools processing
5. **Call variants** - lofreq variant calling
6. **Generate diagnostics** - Depth file + comprehensive diagnostic report

### Outputs
```
{SAMPLE}_results/
  └── {SAMPLE}_depth.txt          # Per-position coverage depth

diagnostic_{SAMPLE}/
  ├── {SAMPLE}_diagnostic_report.txt              # Text summary
  ├── diagnostic_{SAMPLE}_presentation_ready_report.html  # HTML report
  ├── assembly_stats.txt                          # MEGAHIT assembly metrics
  └── blast_results/                              # Contamination detection
      ├── viral_hits.txt
      ├── bacterial_hits.txt
      └── fungal_hits.txt

cleaned_seqs/variants/
  └── {SAMPLE}_vars.vcf           # Raw variant calls
```

### Key Tools
- **fastp** - Quality control
- **bwa** - Read alignment
- **samtools** - BAM processing, depth calculation
- **lofreq** - Variant calling
- **MEGAHIT** - De novo assembly
- **BLAST** - Contamination detection

---

## 🔍 MANUAL DECISION POINT 1: Quality Assessment

### What to Review

#### 1. Depth Coverage (`{SAMPLE}_depth.txt`)
**Questions:**
- Is coverage uniform across the genome?
- Are there large gaps (low/zero coverage regions)?
- What is the median/mean depth?

**Actions:**
```bash
# Quick depth statistics
awk 'NR>1 {sum+=$3; count++} END {print "Mean depth:", sum/count}' {SAMPLE}_results/{SAMPLE}_depth.txt

# Find low coverage regions
awk 'NR>1 && $3 < 50 {print}' {SAMPLE}_results/{SAMPLE}_depth.txt | wc -l
```

#### 2. Contamination Analysis (`diagnostic_{SAMPLE}/`)
**Check Files:**
- `blast_results/viral_hits.txt` - Other viruses present?
- `blast_results/bacterial_hits.txt` - Bacterial contamination?
- `blast_results/fungal_hits.txt` - Fungal contamination?

**Questions:**
- Is the expected virus the dominant hit?
- Are there co-infections (multiple viruses)?
- Is there significant non-viral contamination?

#### 3. Diagnostic Report
**Review:**
- Mapping rate to reference
- Assembly completeness
- Overall quality metrics

### Decision Matrix

| Scenario | Depth | Contamination | Decision |
|----------|-------|---------------|----------|
| **Excellent** | >500x uniform | Target virus only | ✅ Proceed with annotation |
| **Good** | >200x, some gaps | Minimal contamination | ✅ Proceed with annotation |
| **Marginal** | 50-200x, gaps | Some contamination | ⚠️ Proceed with caution, may need parameter adjustment |
| **Poor** | <50x, large gaps | Major contamination | ❌ Consider re-sequencing or rejecting sample |
| **Co-infection** | Good coverage | Multiple viruses | ⚠️ Decide: analyze separately or as mixed population |

### Commands to Proceed (if QC passes)
```bash
# If satisfied with QC, proceed to annotation:
./run_vicast_analyze_annotate_only.sh <R1.fastq.gz> <R2.fastq.gz> <ACCESSION>

# Or submit as SLURM job:
sbatch --wrap="./run_vicast_analyze_annotate_only.sh <R1.fastq.gz> <R2.fastq.gz> <ACCESSION>"
```

---

## CHUNK 2: Variant Annotation (Automated)

### Command
```bash
./run_vicast_analyze_annotate_only.sh <R1.fastq.gz> <R2.fastq.gz> <ACCESSION>
```

### What It Does (Steps 7-9)
**Uses `--resume-from-vcf` flag** - Discovers existing VCF files, skips Steps 1-6 entirely (~17 sec vs hours)

7. **Filter variants** - Apply lofreq quality filters
8. **Annotate with snpEff** - Add functional annotations (gene, effect, amino acid changes)
9. **Parse to TSV** - Convert VCF format to tabular format

### Outputs
```
cleaned_seqs/variants/
  ├── {SAMPLE}.snpEFF.ann.vcf     # Annotated VCF
  ├── {SAMPLE}.snpEFF.ann.tsv     # Annotated TSV (all variants)
  └── {SAMPLE}_200.tsv            # Pre-filtered (depth >= 200)
```

### Key Tool
- **snpEff** - Variant effect prediction (missense, nonsense, synonymous, etc.)

---

## 🔍 MANUAL DECISION POINT 2: Parameter Optimization

### What to Review

#### 1. Variant Distribution (`{SAMPLE}.snpEFF.ann.tsv`)
**Questions:**
- How many total variants were called?
- What is the distribution of variant frequencies (AF)?
- What is the distribution of quality scores (QUAL)?
- What is the depth distribution (DP)?

**Analysis Commands:**
```bash
# Count total variants
grep -v "^#" cleaned_seqs/variants/{SAMPLE}.snpEFF.ann.tsv | wc -l

# Check allele frequency distribution
grep -v "^#" cleaned_seqs/variants/{SAMPLE}.snpEFF.ann.tsv | \
  awk -F'\t' '{print $8}' | grep -oP 'AF=[0-9.]+' | \
  cut -d= -f2 | sort -n | uniq -c

# Check depth distribution
grep -v "^#" cleaned_seqs/variants/{SAMPLE}.snpEFF.ann.tsv | \
  awk -F'\t' '{print $8}' | grep -oP 'DP=[0-9]+' | \
  cut -d= -f2 | sort -n | uniq -c

# Check quality score distribution
grep -v "^#" cleaned_seqs/variants/{SAMPLE}.snpEFF.ann.tsv | \
  awk -F'\t' '{print $6}' | sort -n | uniq -c
```

#### 2. Sample Quality Profile
Based on your analysis from Decision Point 1, determine appropriate thresholds:

**High Quality Sample (>500x depth, clean):**
```bash
--quality 1000    # Strict quality filtering
--depth 200       # High confidence depth
--freq 0.01       # Detect rare variants (1%)
```

**Medium Quality Sample (200-500x depth, some contamination):**
```bash
--quality 500     # Moderate quality filtering
--depth 100       # Moderate confidence depth
--freq 0.02       # Skip very rare variants (2%)
```

**Lower Quality Sample (50-200x depth):**
```bash
--quality 100     # Relaxed quality filtering
--depth 50        # Lower confidence acceptable
--freq 0.05       # Only higher frequency variants (5%)
```

#### 3. Analysis Goals

**Goal: Detect Quasispecies (intra-host diversity)**
- Use low `--freq` threshold (0.01 = 1%)
- Keep high `--quality` and `--depth` to ensure rare variants are real

**Goal: Consensus Genome Only**
- Use high `--freq` threshold (0.50 = 50%) for consensus generation
- Standard `--quality` and `--depth` acceptable

**Goal: SNP Calling (major variants)**
- Use moderate `--freq` threshold (0.20 = 20%)
- High `--quality` and `--depth` for confidence

### Parameter Definitions

| Parameter | Default | Definition | Typical Range |
|-----------|---------|------------|---------------|
| `--quality` | 1000 | Phred-scaled quality score (confidence variant is real) | 100-5000 |
| `--depth` | 200 (parse)<br>200 (consensus) | Minimum reads covering position | 50-500 |
| `--freq` | 0.01 (parse)<br>0.50 (consensus) | Minimum allele frequency (proportion of reads) | 0.01-0.80 |

---

## CHUNK 3: Mutation Analysis & Haplotype Reconstruction

### Script 1: Parse & Filter Mutations

#### Command
```bash
python parse_snpeff_tsv.py \
  cleaned_seqs/variants/{SAMPLE}.snpEFF.ann.tsv \
  {SAMPLE}_results/{SAMPLE}_filtered_mutations.tsv \
  --quality 1000 \
  --depth 200 \
  --freq 0.01
```

#### What It Does
1. Reads snpEff-annotated TSV
2. Extracts AF (allele frequency), DP (depth), DP4 (strand bias) from INFO field
3. Filters variants: `QUAL >= quality` AND `DP >= depth` AND `AF >= freq`
4. Outputs clean TSV with all annotation columns preserved

#### Parameters
- `--quality 1000` - Minimum quality score (default: 1000)
- `--depth 200` - Minimum coverage depth (default: 200)
- `--freq 0.01` - Minimum allele frequency for **quasispecies detection** (default: 0.01 = 1%)

#### Output
```
{SAMPLE}_results/{SAMPLE}_filtered_mutations.tsv
```

**Columns include:**
- CHROM, POS, REF, ALT, QUAL
- Total_Depth, Allele_Frequency, strand_bias, DP4
- EFFECT (e.g., missense_variant, synonymous_variant)
- PUTATIVE_IMPACT (HIGH, MODERATE, LOW, MODIFIER)
- GENE_NAME, FEATURE_ID
- HGVSc (nucleotide change), HGVSp (amino acid change)

**Use Cases:**
- Identify all variants above threshold (including rare quasispecies)
- Analyze mutation spectrum
- Detect selection patterns
- Input for downstream visualization

---

### Script 2: Realistic Viral Haplotype Reconstruction & Quasispecies Analysis

#### Command
```bash
python generate_realistic_haplotype_consensus.py \
  --vcf {SAMPLE}_results/{SAMPLE}_filtered_mutations.tsv \
  --reference cleaned_seqs/{ACCESSION}.fasta \
  --accession {ACCESSION} \
  --quality 1000 \
  --depth 200 \
  --freq 0.50 \
  --output-prefix {SAMPLE}_results/{SAMPLE}_consensus
```

#### What It Does
1. **Re-filters variants** with **majority rule** (`freq >= 0.50` by default)
2. **Applies mutations** to reference genome → consensus sequence
3. **Translates genome** to proteins using virus-specific gene coordinates
4. **Generates individual variant proteins** - creates separate sequences for each mutation
5. **Produces comprehensive summary report**

#### Parameters
- `--quality 1000` - Minimum quality score (default: 1000)
- `--depth 200` - Minimum coverage depth (default: 200)
- `--freq 0.50` - Minimum allele frequency for **consensus inclusion** (default: 0.50 = 50%)
  - **Key difference from Script 1:** Higher frequency = only majority variants in consensus

#### Outputs
```
{SAMPLE}_results/
  ├── {SAMPLE}_consensus_consensus.fasta         # Full genome consensus
  ├── {SAMPLE}_consensus_proteins.fasta          # Individual variant proteins
  └── {SAMPLE}_consensus_summary_report.txt      # Comprehensive report
```

#### `*_consensus.fasta` Contents
```
>{SAMPLE}_filtered_consensus Consensus with N mutations (Q>=1000, D>=200, F>=0.50)
ATGCGATCGATCG... (consensus genome sequence with majority variants applied)
```

#### `*_proteins.fasta` Contents
**Individual protein sequences for each mutation:**
```
>{SAMPLE}_ProteinE ProteinE protein [L287I] (from 862T>A)
MKTLILGAVILGVATAAQITAGIALHQ... (protein with single mutation)

>{SAMPLE}_ProteinM ProteinM protein [V194I] (from 580G>A)
MRHGVALMATLGTLFLFL... (protein with single mutation)

>{SAMPLE}_ProteinNS5 ProteinNS5 protein
MDPWRKGERNGSMKLTY... (reference protein, no mutations)
```

**Key Feature: Quasispecies Support**
- Creates **separate protein sequences** for each individual mutation
- Allows analysis of each variant independently
- Critical for understanding viral diversity and evolution

#### `*_summary_report.txt` Contents
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

ProteinNS5:
  No mutations

--------------------------------------------------------------------------------
INDIVIDUAL PROTEIN VARIANTS
--------------------------------------------------------------------------------

ProteinE: Variant with mutations: L287I
  From nucleotide change: 862T>A

ProteinE: Variant with mutations: K310R
  From nucleotide change: 929A>G

ProteinM: Variant with mutations: V194I
  From nucleotide change: 580G>A

ProteinNS5: Reference sequence (no mutations)
```

---

## Complete Workflow Example: WNV Analysis

### Step 1: Run QC (Chunk 1)
```bash
./run_vicast_analyze_qc_only.sh \
  WNV_SMS14_R1.fastq.gz \
  WNV_SMS14_R2.fastq.gz \
  AY532665.1 \
  8
```

**Outputs:**
- `WNV_SMS14_results/WNV_SMS14_depth.txt`
- `diagnostic_WNV_SMS14/WNV_SMS14_diagnostic_report.txt`
- `cleaned_seqs/variants/WNV_SMS14_vars.vcf`

---

### Step 2: Manual Review (Decision Point 1)

#### Check Coverage
```bash
# Mean depth
awk 'NR>1 {sum+=$3; count++} END {print "Mean depth:", sum/count}' WNV_SMS14_results/WNV_SMS14_depth.txt

# Low coverage positions (< 50x)
awk 'NR>1 && $3 < 50 {print}' WNV_SMS14_results/WNV_SMS14_depth.txt | wc -l
```

**Result:** Mean depth = 856x, only 12 positions < 50x → **Excellent coverage** ✅

#### Check Contamination
```bash
# Review diagnostic report
cat diagnostic_WNV_SMS14/WNV_SMS14_diagnostic_report.txt

# Check BLAST results
head diagnostic_WNV_SMS14/blast_results/viral_hits.txt
```

**Result:** Top hit = West Nile virus (AY532665.1), no significant contamination → **Clean sample** ✅

**Decision:** ✅ Proceed with annotation

---

### Step 3: Run Annotation (Chunk 2)
```bash
./run_vicast_analyze_annotate_only.sh \
  WNV_SMS14_R1.fastq.gz \
  WNV_SMS14_R2.fastq.gz \
  AY532665.1
```

**Outputs:**
- `cleaned_seqs/variants/WNV_SMS14.snpEFF.ann.vcf`
- `cleaned_seqs/variants/WNV_SMS14.snpEFF.ann.tsv`

**Time:** ~17 seconds (vs hours for full pipeline!)

---

### Step 4: Evaluate Parameters (Decision Point 2)

#### Check Variant Distribution
```bash
# Total variants
grep -v "^#" cleaned_seqs/variants/WNV_SMS14.snpEFF.ann.tsv | wc -l
# Result: 3,165 variants

# Allele frequency distribution
grep -v "^#" cleaned_seqs/variants/WNV_SMS14.snpEFF.ann.tsv | \
  awk -F'\t' '{print $8}' | grep -oP 'AF=[0-9.]+' | cut -d= -f2 | \
  awk '{if ($1 < 0.05) low++; else if ($1 < 0.20) med++; else high++}
       END {print "Low (<5%):", low, "Medium (5-20%):", med, "High (>20%):", high}'
# Result: Low: 2,890, Medium: 180, High: 95
```

**Analysis:**
- Most variants are low frequency (<5%) - strong quasispecies signal
- High depth (856x) supports detection of rare variants
- Sample quality is excellent

**Decision:** Use **strict parameters** to capture rare variants confidently
```bash
--quality 1000  # High confidence only
--depth 200     # Well above average depth
--freq 0.01     # Capture rare quasispecies (1%)
```

---

### Step 5: Parse & Filter Mutations (Chunk 3, Script 1)
```bash
conda activate vicast_analyze

python parse_snpeff_tsv.py \
  cleaned_seqs/variants/WNV_SMS14.snpEFF.ann.tsv \
  WNV_SMS14_results/WNV_SMS14_filtered_mutations.tsv \
  --quality 1000 \
  --depth 200 \
  --freq 0.01
```

**Output:**
```
Processing SnpEff TSV file: cleaned_seqs/variants/WNV_SMS14.snpEFF.ann.tsv
Quality cutoff: 1000
Depth cutoff: 200
Allele frequency cutoff: 0.01
Found 3165 total variants
After filtering (QUAL >= 1000, depth >= 200, freq >= 0.01): 287 variants
✅ Successfully parsed and filtered 287 mutations to WNV_SMS14_results/WNV_SMS14_filtered_mutations.tsv
```

**Result:** 287 high-confidence variants (including rare quasispecies ≥1%)

---

### Step 6: Generate Haplotype Consensus (Chunk 3, Script 2)
```bash
python generate_realistic_haplotype_consensus.py \
  --vcf WNV_SMS14_results/WNV_SMS14_filtered_mutations.tsv \
  --reference cleaned_seqs/AY532665.1.fasta \
  --accession AY532665.1 \
  --quality 1000 \
  --depth 200 \
  --freq 0.50 \
  --output-prefix WNV_SMS14_results/WNV_SMS14_consensus
```

**Output:**
```
================================================================================
INDIVIDUAL VARIANT PROTEIN GENERATOR
================================================================================
Reference genome length: 11029 bp
Found configuration for West Nile virus
Genes defined: 10
After filtering: 15 mutations pass criteria

🧬 Generating individual variant proteins...
Generated 18 unique protein variants

📝 Writing consensus genome...
Wrote consensus genome to: WNV_SMS14_results/WNV_SMS14_consensus_consensus.fasta

🧬 Writing protein sequences...
Wrote 18 protein sequences to: WNV_SMS14_results/WNV_SMS14_consensus_proteins.fasta

📄 Summary report saved to: WNV_SMS14_results/WNV_SMS14_consensus_summary_report.txt

✅ Individual variant protein generation complete!
Generated 18 unique protein variants
```

**Key Result:**
- 287 variants ≥1% frequency (quasispecies)
- 15 variants ≥50% frequency (consensus)
- 18 unique protein sequences (each mutation creates separate protein)

---

## Summary: Three-Chunk Architecture

| Chunk | Scripts | Type | Manual Review | Purpose |
|-------|---------|------|---------------|---------|
| **1** | `run_vicast_analyze_qc_only.sh` | Automated | → **Decision Point 1** | QC, mapping, variant calling, diagnostics |
| ↓ | | | **Depth & contamination** | |
| **2** | `run_vicast_analyze_annotate_only.sh` | Automated | → **Decision Point 2** | Variant annotation with snpEff |
| ↓ | | | **Parameter optimization** | |
| **3** | `parse_snpeff_tsv.py`<br>`generate_realistic_haplotype_consensus.py` | Manual (with parameters) | Final outputs | Mutation filtering, consensus, quasispecies analysis |

---

## Key Insights

### Why Two Different Frequency Thresholds?

**Script 1 (`parse_snpeff_tsv.py`):** `--freq 0.01` (1%)
- **Goal:** Capture ALL biologically relevant variants
- **Output:** Complete mutation spectrum including rare quasispecies
- **Use:** Diversity analysis, selection studies, evolutionary analysis

**Script 2 (`generate_realistic_haplotype_consensus.py`):** `--freq 0.50` (50%)
- **Goal:** Generate consensus representing dominant sequence
- **Output:** Majority-rule consensus + individual variant proteins
- **Use:** Reference sequence, protein structure prediction, phylogenetics

### Why Individual Variant Proteins?
Traditional consensus methods create a **single chimeric protein** that may not actually exist in any viral particle. The new approach creates **separate proteins for each mutation**, which:
- Represents actual quasispecies present in the sample
- Allows structural analysis of each variant independently
- Supports realistic evolutionary modeling
- Enables accurate fitness predictions

---

**Last Updated:** 2025-11-25
**Pipeline Version:** vicast-analyze (main branch)
**Key Feature:** `--resume-from-vcf` enables efficient workflow separation
