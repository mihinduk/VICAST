# VICAST Frequency-Stratified Consensus Generation - User Guide

**Module:** `generate_realistic_haplotype_consensus.py`

**Version:** 2.2.0

**Last Updated:** 2026-02-04

---

## Table of Contents

1. [Overview](#overview)
2. [What This Tool Does (and Doesn't Do)](#what-this-tool-does-and-doesnt-do)
3. [Scientific Background](#scientific-background)
4. [Parameters and Recommendations](#parameters-and-recommendations)
5. [Use Cases and Examples](#use-cases-and-examples)
6. [Interpreting Results](#interpreting-results)
7. [Limitations and Caveats](#limitations-and-caveats)
8. [Advanced Usage](#advanced-usage)

---

## Overview

This tool generates **frequency-stratified consensus sequences** from viral variant data. It groups variants by their frequency in the population and generates representative consensus sequences for different frequency tiers.

### What Problem Does This Solve?

**Standard single consensus genomes** have a limitation:
- They mask population diversity by reporting only the majority allele at each position
- A sample with mutations at 95%, 60%, and 15% all get treated differently
- You lose information about population structure

**This tool generates multiple consensi** representing:
- **Dominant strain** (high-frequency mutations, ≥95%)
- **Major variants** (emerging mutations, 80-94%)
- **Minority variants** (moderate frequency, 50-79%)
- **Low-frequency variants** (<50%, exploratory)

### Important Distinction

This tool performs **frequency-based grouping with statistical inference**, NOT true haplotype reconstruction.

**What it does:**
- Groups variants by frequency
- For high-frequency variants (≥95%), makes statistically justified assumption of linkage
- Generates consensus sequences representing different frequency tiers

**What it does NOT do:**
- Check BAM files for read-level co-occurrence
- Use paired-end reads to prove mutations are on the same genome
- Employ haplotype phasing algorithms (like PredictHaplo, QuasiRecomb, etc.)

---

## What This Tool Does (and Doesn't Do)

### ✅ What It DOES Do

| Capability | Description | Confidence Level |
|------------|-------------|------------------|
| **Dominant consensus** | Generates consensus with high-frequency mutations (≥95%) | HIGH - Statistically justified |
| **Frequency stratification** | Groups variants into tiers (dominant/major/minor/low) | HIGH - Direct from data |
| **Population diversity assessment** | Identifies presence of multiple variant frequencies | HIGH - Observational |
| **Exploratory haplotypes** | Suggests possible variant combinations | LOW - Assumption-based |
| **Passage tracking** | Monitors emergence and fixation of variants | HIGH - Frequency changes |

### ❌ What It Does NOT Do

| Limitation | Why Not | Alternative |
|------------|---------|-------------|
| **Prove linkage of minority variants** | Requires read-level evidence | Use PredictHaplo, QuasiRecomb, ViQuaS |
| **Determine exact haplotype frequencies** | Needs phasing information | Long-read sequencing + phasing tools |
| **Identify reassortment** | Requires cross-segment linkage | Specialized reassortment detection tools |
| **Distinguish PCR duplicates from biological replicates** | Needs UMI barcodes | Use UMI-based deduplication |

---

## Scientific Background

### The Statistical Basis for Linkage Inference

#### High-Frequency Variants (≥95%): High Confidence Linkage

**The Pigeonhole Principle** provides mathematical certainty:

```
Example:
- Sample has 1000 viral genomes
- Mutation A: 95% frequency (950 genomes)
- Mutation B: 95% frequency (950 genomes)

Question: How many genomes MUST have both mutations?

Answer: At minimum, 900 genomes (90%)

Proof:
- If you try to maximize separation:
  - 50 genomes with only A
  - 50 genomes with only B
  - Leaves 900 genomes that must have both A and B
```

**Conclusion:** At ≥95% frequency, linkage is **statistically certain** (≥90% must have both).

This is the **same principle that justifies standard consensus genome generation**.

#### Medium Frequency (80-94%): Moderate Confidence

```
Example:
- Mutation A: 80% (800/1000)
- Mutation B: 80% (800/1000)

Minimum overlap: 600 genomes (60%)
Maximum separation: 200 with A only, 200 with B only

Conclusion: Probable linkage, but 40% uncertainty
```

**Use case:** Tracking adaptive mutations sweeping through population during passage.

**Caveat:** State assumption in publications.

#### Low Frequency (50-79%): Cannot Assume Linkage

```
Example:
- Mutation A: 60% (600/1000)
- Mutation B: 60% (600/1000)

Minimum overlap: 200 genomes (20%)
Maximum overlap: 600 genomes (60%)

Conclusion: Complete uncertainty about linkage
```

**Use case:** Report individual variant frequencies only.

**Caveat:** Do NOT claim these form haplotypes without read-level evidence.

#### Very Low Frequency (<50%): No Linkage Inference

At frequencies below 50%, there is no statistical basis for linkage inference. These are minority variants; their combinations with other variants are unknown.

---

## Parameters and Recommendations

### Quality Filtering Parameters

#### `--quality` (default: 1000)

**What it does:** Minimum PHRED-scaled quality score for variant calls

**Biological context:**
- Quality 1000 = 99.9% confidence
- Quality 500 = 99.7% confidence
- Quality 100 = 99% confidence

**Recommendations by virus type:**

| Virus Type | Recommended Quality | Rationale |
|------------|-------------------|-----------|
| **ssRNA** (SARS-CoV-2, dengue, influenza) | 1000-5000 | High mutation rate; need to distinguish real variants from errors |
| **dsRNA** (rotavirus, reovirus) | 1000-2000 | Moderate mutation rate |
| **ssDNA** (parvovirus, circovirus) | 500-1000 | Lower mutation rate; less stringency needed |
| **dsDNA** (herpes, pox) | 500-1000 | Very stable genomes |

**When to adjust:**
- **Increase (5000+):** Deep sequencing, need very high confidence
- **Decrease (100-500):** Low coverage samples, well-characterized virus

---

#### `--depth` (default: 200)

**What it does:** Minimum read depth at variant position

**Biological context:**
- At depth 200, a 1% variant = 2 reads (may be real or error)
- At depth 200, a 5% variant = 10 reads (more confident)
- Higher depth = better frequency estimation

**Recommendations:**

| Coverage Level | Recommended Depth | Min Detectable Frequency |
|----------------|-------------------|-------------------------|
| Low (100-500x) | 50-100 | ~5-10% |
| Medium (500-2000x) | 200 | ~1% |
| High (2000-10000x) | 500 | ~0.5% |
| Ultra-deep (>10000x) | 1000 | ~0.1% |

**When to adjust:**
- **Increase:** Deep sequencing experiments, need confident rare variants
- **Decrease:** Low-coverage samples, risk missing real variants

---

#### `--freq` (default: 0.01 = 1%)

**What it does:** Minimum allele frequency to include a variant

**Biological context:**
- Illumina error rate: ~0.1-0.5% per base
- True minority variants in ssRNA viruses: often 1-40%
- Setting too low: include sequencing errors as "variants"
- Setting too high: miss real minority variants

**Recommendations by virus type:**

| Virus Type | Recommended Frequency | Rationale |
|------------|---------------------|-----------|
| **ssRNA** (high mutation) | 0.01-0.02 (1-2%) | Real quasispecies variants common |
| **dsRNA** | 0.02-0.05 (2-5%) | Moderate diversity |
| **ssDNA** | 0.05 (5%) | Fewer real variants |
| **dsDNA** (low mutation) | 0.05-0.10 (5-10%) | Very stable; low-freq likely errors |

**When to adjust:**
- **Decrease (0.005 = 0.5%):** High-quality sequencing, need to detect rare variants
- **Increase (0.05 = 5%):** Avoid false positives in stable genomes

---

### Frequency Threshold Parameters (NEW)

#### `--dominant-threshold` (default: 0.95 = 95%)

**What it does:** Minimum frequency for "dominant haplotype" with HIGH confidence linkage

**Statistical basis:** At 95%, pigeonhole principle guarantees ≥90% of genomes have all mutations

**Biological interpretation:**
- These mutations have **fixed or nearly fixed** in the population
- Statistically justified to assume co-occurrence
- Equivalent to majority consensus genome generation
- **Safe for publication claims**

**When to adjust:**
- **Increase to 0.99 (99%):** Need very high confidence (e.g., clinical samples)
- **Decrease to 0.90 (90%):** Acceptable with caveat; still 80% minimum overlap

**Publication language:**
> "Dominant consensus generated from high-frequency variants (≥95%) with statistical inference of linkage (pigeonhole principle: ≥90% minimum co-occurrence)."

---

#### `--major-threshold` (default: 0.80 = 80%)

**What it does:** Minimum frequency for "major variant" with MODERATE confidence linkage

**Statistical basis:** At 80% for two mutations, minimum 60% of genomes must have both

**Biological interpretation:**
- Mutations are **sweeping through population** but not yet fixed
- Common in passage studies (adaptive mutations emerging)
- Probable linkage, but 40% uncertainty
- **Use with caution in publications**

**When to adjust:**
- **Increase to 0.85-0.90:** More conservative, less uncertainty
- **Decrease to 0.70:** More exploratory, but state limitations

**Publication language:**
> "Major variants identified at 80-95% frequency. Linkage is assumed based on frequency (minimum 60% co-occurrence) but not proven. These represent emerging adaptive mutations during passage."

---

#### `--minor-threshold` (default: 0.50 = 50%)

**What it does:** Minimum frequency for "minority variant" reporting

**Statistical basis:** At 50%, NO inference about linkage with other variants possible

**Biological interpretation:**
- These are **minority variants** in the population
- Represent population diversity
- **Cannot assume these are linked to other mutations**
- Useful for tracking emergence over passages

**When to adjust:**
- **Increase to 0.60-0.70:** Focus on more abundant minorities
- **Decrease to 0.40:** Capture more of the variant spectrum

**Publication language:**
> "Minority variants (50-95%) were detected and reported individually. Linkage between minority variants is not inferred without read-level phasing evidence."

---

#### `--min-haplotype-freq` (default: 0.01 = 1%)

**What it does:** Minimum estimated frequency to report a putative haplotype

**Why it exists:**
- Prevents reporting of very rare combinations
- Reduces noise in output
- Focuses on biologically significant variants

**When to adjust:**
- **Increase to 0.05 (5%):** Focus on major populations only
- **Decrease to 0.005 (0.5%):** Report rare putative haplotypes (exploratory)

---

### Biological Differences Across Virus Types

Understanding your virus type is critical for choosing appropriate parameters. Different genome types have fundamentally different mutation rates and population structures.

#### Complete Comparison Table

| Property | ssRNA (+/−) | dsRNA | ssDNA | dsDNA |
|----------|-------------|-------|-------|-------|
| **Examples** | SARS-CoV-2, dengue, influenza, HIV, poliovirus | Rotavirus, reovirus, bluetongue virus | Parvovirus, circovirus, anellovirus, AAV | Herpesvirus, poxvirus, adenovirus |
| **Mutation rate** | High (10^-4 to 10^-6) | Intermediate (10^-5 to 10^-7) | Variable (10^-6 to 10^-9)* | Low (10^-7 to 10^-9) |
| **Polymerase/Replication** | RdRp (no proofreading) | RdRp (some correction) | Host DNA pol OR rolling circle* | DNA pol (3' exonuclease proofreading) |
| **Quasispecies diversity** | Large, diverse populations | Moderate diversity | Low to moderate | More homogeneous |
| **Minority variants** | Very common (5-40%) | Moderate (5-20%) | Uncommon (1-10%) | Rare (<5%) |
| **Population structure** | Complex quasispecies | Moderate diversity | Relatively simple | Usually clonal |
| **Sequencing error vs real variants** | Harder to distinguish | Moderate difficulty | Moderate difficulty | Easier to distinguish |
| **Genome stability** | Low (highly mutable) | Moderate | Moderate to high | High (very stable) |
| **Passage adaptation** | Rapid, frequent | Moderate rate | Slower | Rare/slow |

\*ssDNA mutation rate varies greatly depending on replication strategy:
- Using host DNA polymerase (with proofreading): Low (10^-8 to 10^-9)
- Using rolling circle replication: Higher (10^-6 to 10^-7)

#### Parameter Recommendations by Virus Type

| Parameter | ssRNA | dsRNA | ssDNA | dsDNA |
|-----------|-------|-------|-------|-------|
| **`--quality`** | 1000-5000 | 1000-2000 | 500-1000 | 500-1000 |
| **`--depth`** | 200+ | 200+ | 200+ | 200+ |
| **`--freq` (detection)** | 0.01-0.02 (1-2%) | 0.02-0.03 (2-3%) | 0.03-0.05 (3-5%) | 0.05-0.10 (5-10%) |
| **`--dominant-threshold`** | 0.95 | 0.95 | 0.95 | 0.95 |
| **`--major-threshold`** | 0.80 | 0.80-0.85 | 0.85-0.90 | 0.90 |
| **`--minor-threshold`** | 0.50 | 0.50-0.60 | 0.60-0.70 | 0.70 |

---

### ⚠️ Important Warnings for Specific Virus Types

#### Segmented Genomes (dsRNA and some ssRNA)

**Viruses affected:** Rotavirus (11 segments), Reovirus (10 segments), Influenza (8 segments), Bunyaviruses (3 segments)

**CRITICAL WARNING:**

> **Linkage assumptions do NOT apply across segments!**
>
> A mutation in segment 1 at 95% frequency and a mutation in segment 5 at 95% frequency do NOT necessarily co-occur on the same viral particle. Reassortment can shuffle segments between viruses.

**Recommendations:**
1. Run analysis per-segment
2. Interpret each segment independently
3. Do NOT assume mutations across segments are linked
4. For reassortment analysis, use specialized tools (not VICAST)

**Example workflow:**
```bash
# Analyze each segment separately
for segment in PB2 PB1 PA HA NP NA M NS; do
    python generate_realistic_haplotype_consensus.py \
        --vcf ${segment}_variants.vcf \
        --reference ${segment}_ref.fasta \
        --accession ${segment}_accession \
        --output-prefix ${segment}_consensus
done
```

---

#### ssDNA Viruses (Replication Strategy Matters)

**IMPORTANT:** ssDNA virus mutation rates vary dramatically based on replication mechanism.

**Two Categories:**

**1. Host DNA Polymerase Replication** (Low mutation rate)
- Examples: AAV, some parvoviruses
- Host DNA polymerase has 3' exonuclease proofreading
- Mutation rate: ~10^-8 to 10^-9 (similar to dsDNA)
- **Use stringent settings:** `--freq 0.05`, `--quality 500`

**2. Rolling Circle Replication** (Higher mutation rate)
- Examples: Circoviruses, some parvoviruses
- Less fidelity, more errors
- Mutation rate: ~10^-6 to 10^-7 (closer to dsRNA)
- **Use intermediate settings:** `--freq 0.03`, `--quality 1000`

**If unsure which mechanism your virus uses:**
- Start with stringent settings (`--freq 0.05`)
- Examine variant frequency distribution
- If many variants at 2-5%, consider using intermediate settings

---

### Virus-Specific Recommendations

#### ssRNA Viruses (SARS-CoV-2, Dengue, Influenza, HIV)

**Characteristics:**
- High mutation rate (10^-4 to 10^-6 per base per cycle)
- Large quasispecies populations
- Minority variants common (5-40%)
- Passage adaptation common

**Recommended settings:**
```bash
python generate_realistic_haplotype_consensus.py \
    --vcf variants.vcf \
    --reference genome.fasta \
    --accession NC_001477 \
    --quality 1000 \
    --depth 200 \
    --freq 0.01 \
    --dominant-threshold 0.95 \
    --major-threshold 0.80 \
    --minor-threshold 0.50 \
    --output-prefix sample_consensus
```

**Interpretation:**
- **Dominant haplotype (≥95%):** Fixed mutations, high confidence
- **Major variants (80-94%):** Emerging adaptive mutations (common in passages)
- **Minor variants (50-79%):** Population diversity, track over time
- **Low-freq (<50%):** Real quasispecies variants, but linkage unknown

---

#### dsRNA Viruses (Rotavirus, Reovirus)

**Characteristics:**
- Intermediate mutation rate (10^-5 to 10^-7)
- Moderate quasispecies diversity
- Segmented genomes (requires multi-segment analysis)

**Recommended settings:**
```bash
python generate_realistic_haplotype_consensus.py \
    --vcf variants.vcf \
    --reference genome.fasta \
    --accession NC_011500 \
    --quality 1000 \
    --depth 200 \
    --freq 0.02 \
    --dominant-threshold 0.95 \
    --major-threshold 0.80 \
    --minor-threshold 0.50 \
    --output-prefix sample_consensus
```

**Key difference:** Slightly higher `--freq 0.02` to reduce false positives

---

#### ssDNA Viruses (Parvovirus, Circovirus, AAV)

**Characteristics:**
- Variable mutation rate (10^-6 to 10^-9)
- Lower quasispecies diversity
- Minority variants less common

**Recommended settings:**
```bash
python generate_realistic_haplotype_consensus.py \
    --vcf variants.vcf \
    --reference genome.fasta \
    --accession NC_001401 \
    --quality 500 \
    --depth 200 \
    --freq 0.05 \
    --dominant-threshold 0.95 \
    --major-threshold 0.85 \
    --minor-threshold 0.60 \
    --output-prefix sample_consensus
```

**Key differences:**
- Lower quality threshold (500) - fewer real variants expected
- Higher freq threshold (0.05) - reduce false positives
- Higher thresholds overall - more homogeneous populations

---

#### dsDNA Viruses (Herpesvirus, Poxvirus, Adenovirus)

**Characteristics:**
- Low mutation rate (10^-7 to 10^-9)
- Very stable genomes
- Minimal quasispecies diversity
- Low-frequency "variants" often sequencing errors

**Recommended settings:**
```bash
python generate_realistic_haplotype_consensus.py \
    --vcf variants.vcf \
    --reference genome.fasta \
    --accession NC_001806 \
    --quality 500 \
    --depth 200 \
    --freq 0.05 \
    --dominant-threshold 0.95 \
    --major-threshold 0.85 \
    --minor-threshold 0.60 \
    --output-prefix sample_consensus
```

**Key differences:**
- Higher freq threshold (0.05) - critical to avoid errors as "variants"
- Most variants should be high frequency (fixed mutations)
- Expect fewer minority variants

---

### Quick Reference: Common Viruses

Use this table for quick parameter selection for frequently studied viruses:

| Virus | Type | Genome | `--freq` | `--quality` | Special Considerations |
|-------|------|--------|----------|-------------|----------------------|
| **SARS-CoV-2** | ssRNA+ | 30kb | 0.01 | 1000-2000 | High diversity, common quasispecies |
| **Influenza A** | ssRNA− | 8 seg | 0.01 | 1000 | **SEGMENTED** - analyze per segment |
| **Dengue virus** | ssRNA+ | 11kb | 0.01 | 1000-2000 | High quasispecies diversity |
| **Zika virus** | ssRNA+ | 11kb | 0.01 | 1000 | Similar to dengue |
| **HIV-1** | ssRNA+ (retro) | 10kb | 0.005-0.01 | 2000-5000 | Extremely high diversity |
| **Poliovirus** | ssRNA+ | 7.5kb | 0.01 | 1000 | High quasispecies diversity |
| **Hepatitis C** | ssRNA+ | 9.6kb | 0.01 | 1000-2000 | High quasispecies diversity |
| **Rotavirus** | dsRNA | 11 seg | 0.02 | 1000 | **SEGMENTED** - 11 segments, analyze separately |
| **Reovirus** | dsRNA | 10 seg | 0.02 | 1000 | **SEGMENTED** - 10 segments, analyze separately |
| **Parvovirus** | ssDNA | 5kb | 0.05 | 500 | Host DNA pol, low diversity |
| **AAV** | ssDNA | 4.7kb | 0.05 | 500 | Host DNA pol, very stable |
| **Circovirus** | ssDNA | 2kb | 0.03 | 1000 | Rolling circle, moderate rate |
| **Anellovirus** | ssDNA | 3.8kb | 0.03-0.05 | 500-1000 | Variable diversity |
| **Herpes simplex** | dsDNA | 152kb | 0.05 | 500 | Very stable, low diversity |
| **Varicella-zoster** | dsDNA | 125kb | 0.05 | 500 | Very stable |
| **Cytomegalovirus** | dsDNA | 235kb | 0.05 | 500 | Very stable |
| **Vaccinia virus** | dsDNA | 190kb | 0.05 | 500 | Poxvirus, very stable |
| **Adenovirus** | dsDNA | 36kb | 0.05 | 500 | Very stable |

**Notes:**
- **SEGMENTED viruses:** Run analysis separately for each segment
- **Retroviruses (HIV):** Go through DNA intermediate, use very stringent quality
- **Circovirus:** Rolling circle replication = higher mutation rate than other ssDNA

---

### Future Enhancement: Virus Type Presets

**Coming soon:** Automatic parameter selection based on virus type

```bash
# Instead of manually setting all parameters:
python generate_realistic_haplotype_consensus.py \
    --vcf variants.vcf \
    --reference genome.fasta \
    --virus-type ssRNA \
    --output-prefix sample

# Automatically sets:
# --quality 1000
# --freq 0.01
# --dominant-threshold 0.95
# --major-threshold 0.80
# --minor-threshold 0.50
```

**Planned presets:**
- `--virus-type ssRNA` (default for RNA viruses)
- `--virus-type dsRNA` (for rotavirus, reovirus)
- `--virus-type ssDNA` (for parvovirus, AAV)
- `--virus-type dsDNA` (for herpes, pox, adeno)
- `--virus-type ssRNA-segmented` (for influenza)
- `--virus-type retrovirus` (for HIV - extra stringent)

Users can still override individual parameters:
```bash
--virus-type ssRNA --quality 5000  # Use ssRNA preset but higher quality
```

---

## Use Cases and Examples

### Use Case 1: Passage Study (SARS-CoV-2)

**Scenario:** Tracking adaptive mutations across 10 passages

**Goal:** Identify emerging mutations and monitor their frequency changes

**Settings:**
```bash
--quality 1000 --depth 200 --freq 0.01
--dominant-threshold 0.95
--major-threshold 0.80
--minor-threshold 0.50
```

**Expected output:**
- **Passage 1:** Wild-type (100%), maybe a few low-freq variants (1-5%)
- **Passage 5:** Major variant emerges (40-80%)
- **Passage 10:** Major variant → Dominant (>95%), original wild-type gone

**Interpretation:**
> "An adaptive mutation in the spike protein emerged at passage 3 (15%), became a major variant by passage 5 (75%), and fixed by passage 10 (99%). The frequency trajectory suggests positive selection."

---

### Use Case 2: Clinical Sample (Dengue Virus)

**Scenario:** Patient sample with quasispecies diversity

**Goal:** Characterize population structure, identify drug resistance mutations

**Settings:**
```bash
--quality 2000 --depth 500 --freq 0.01
--dominant-threshold 0.95
--major-threshold 0.80
--minor-threshold 0.50
```

**Expected output:**
- Dominant consensus (95%+): Major circulating strain
- Several minority variants (5-40%): Quasispecies diversity
- Possibly drug resistance at 10-30%

**Interpretation:**
> "The dominant viral strain (95% frequency) contained no drug resistance mutations. However, a minority variant with NS5 I223V (associated with drug resistance) was detected at 12% frequency, suggesting potential for treatment failure."

**Caveat:**
> "Linkage of the resistance mutation with other variants cannot be determined without haplotype phasing. The 12% frequency indicates presence in at least 12% of the viral population."

---

### Use Case 3: Vaccine Strain Purity (Influenza)

**Scenario:** Quality control for vaccine production

**Goal:** Ensure minimal variant contamination

**Settings:**
```bash
--quality 5000 --depth 1000 --freq 0.005
--dominant-threshold 0.99
--major-threshold 0.95
--minor-threshold 0.80
```

**Expected output:**
- Dominant consensus (>99%): Vaccine strain
- Ideally, no variants >1%

**Interpretation:**
> "The vaccine preparation showed >99.5% sequence identity to the reference strain. Two minority variants were detected at 0.8% and 0.6% frequency, both synonymous substitutions with no phenotypic impact. Quality control passed."

---

## Interpreting Results

### Output Files

The tool generates three main output files:

#### 1. `[prefix]_realistic_haplotypes.fasta`

**Content:** Consensus sequences for each frequency tier

```fasta
>Sample_Dominant
ATGGCTAGC... (95%+ frequency mutations applied)

>Sample_Major_variant_1
ATGGCAAGC... (80-95% mutations applied)

>Sample_Wildtype
ATGGCTAGC... (reference, estimated remaining frequency)
```

**How to interpret:**
- **Dominant:** High confidence this represents the major strain
- **Major variant:** Probable emerging strain, state assumption
- **Wild-type:** Estimated frequency of un-mutated reference

---

#### 2. `[prefix]_realistic_proteins.fasta`

**Content:** Protein sequences translated from each consensus

```fasta
>Sample_Dominant_NS5
MFKFPGFLW... [L287I, V345M] (Dominant, est. 89.50%)

>Sample_Major_variant_1_NS5
MFKFPGLFW... [V345M] (Major_variant_1, est. 8.30%)
```

**How to interpret:**
- Each protein shows which amino acid changes it contains
- Frequency estimate shows abundance in population
- Compare across passages to track protein evolution

---

#### 3. `[prefix]_realistic_haplotype_report.txt`

**Content:** Detailed summary with mutation-by-mutation breakdown

**Key sections:**
- Total mutations detected
- Haplotype frequency estimates
- Mutations per haplotype
- Protein variant summary

**Use this for:** Detailed analysis and supplementary materials

---

### Statistical Confidence Interpretation

| Frequency Category | Linkage Confidence | Publication Language |
|-------------------|-------------------|----------------------|
| **≥95%** (Dominant) | **HIGH** | "High-frequency variants (≥95%) were assumed linked based on statistical certainty (pigeonhole principle guarantees ≥90% co-occurrence)." |
| **80-94%** (Major) | **MODERATE** | "Major variants (80-94%) represent probable linked mutations based on frequency, though linkage is not proven. Minimum 60% co-occurrence is guaranteed." |
| **50-79%** (Minor) | **LOW** | "Minority variants (50-79%) were reported individually. Linkage with other mutations cannot be inferred." |
| **<50%** (Low-freq) | **NONE** | "Low-frequency variants (<50%) detected. These may represent true quasispecies diversity, but no haplotype inference was performed." |

---

## Limitations and Caveats

### What You Can Claim

✅ **Acceptable claims:**
- "Dominant consensus generated from high-frequency variants (≥95%) using frequency-based inference"
- "Mutations at ≥95% frequency are statistically likely to co-occur (≥90% minimum)"
- "Population structure assessed through frequency stratification"
- "Tracking emergence and fixation of adaptive mutations across passages"

### What You Cannot Claim

❌ **Unacceptable claims without qualification:**
- "Haplotype reconstruction" (implies read-level phasing)
- "Identified competing viral lineages" (can't prove linkage of minority variants)
- "Determined exact haplotype frequencies" (assumptions involved)
- "Quasispecies reconstruction" (technical term meaning proven co-occurrence)

### When You Need True Haplotype Reconstruction

Use specialized phasing tools if you need to:
- **Prove linkage** of minority variants (<95%)
- **Determine exact haplotype** frequencies with confidence intervals
- **Identify reassortment** events (requires cross-segment phasing)
- **Track specific mutation combinations** through population

**Recommended tools:**
- **PredictHaplo:** Maximum likelihood haplotype reconstruction
- **QuasiRecomb:** Recombination-aware haplotype assembly
- **ViQuaS:** Viral quasispecies spectrum reconstruction
- **ShoRAH:** Short read assembly into haplotypes
- **Long-read sequencing:** PacBio/Nanopore for direct phasing

---

## Advanced Usage

### Combining with Read-Level Phasing (Future Feature)

**Proposed enhancement:** Check BAM files for variants within read/insert distance

For variants within 150bp (read length) or 500bp (insert size):
```bash
python check_read_level_cooccurrence.py \
    --bam aligned.bam \
    --vcf variants.vcf \
    --output cooccurrence.tsv
```

This would enable:
- ✅ **Partial haplotype reconstruction** for nearby variants
- ✅ Claims like "mutations at positions 100 and 150 co-occur on 45% of reads"
- ✅ Validation of frequency-based assumptions

---

### Integration with Phylogenetics

Consensus sequences can be used for phylogenetic analysis:

```bash
# Generate consensi for multiple samples/passages
for passage in P1 P2 P3 P4 P5; do
    python generate_realistic_haplotype_consensus.py \
        --vcf ${passage}.vcf \
        --reference ref.fasta \
        --accession NC_001477 \
        --dominant-threshold 0.95 \
        --output-prefix ${passage}
done

# Combine dominant consensi
cat P*_Dominant_realistic_haplotypes.fasta > all_passages.fasta

# Build phylogeny
mafft --auto all_passages.fasta > aligned.fasta
iqtree -s aligned.fasta -m TEST -bb 1000
```

---

## Best Practices Summary

### For Publication-Quality Analysis

1. **Use conservative thresholds:**
   - `--dominant-threshold 0.95` (default)
   - `--quality 1000+` for ssRNA, 500+ for DNA viruses
   - `--depth 200+` minimum

2. **State your methods clearly:**
   - Describe frequency-based inference
   - Cite pigeonhole principle for high-frequency linkage
   - Acknowledge limitations for minority variants

3. **Validate with independent methods when possible:**
   - Sanger sequencing for specific mutations
   - Replicate passage experiments
   - Compare with other bioinformatic tools

4. **Report transparently:**
   - Include all thresholds used
   - Report variant frequencies alongside consensi
   - State assumptions explicitly

---

## Troubleshooting

### Problem: Too many low-frequency variants

**Symptoms:** Output cluttered with many variants <5%

**Solutions:**
- Increase `--freq` to 0.05 (5%)
- Increase `--quality` to 2000+
- Check sequencing quality metrics

---

### Problem: No dominant haplotype

**Symptoms:** All variants are <95%

**Possible causes:**
1. **Mixed infection** (multiple strains)
2. **Highly diverse quasispecies** (true biology)
3. **Poor quality data** (sequencing issues)

**Solutions:**
- Lower `--dominant-threshold` to 0.85-0.90 (report as exploratory)
- Check for contamination
- Inspect read quality and coverage

---

### Problem: Unrealistic frequency estimates

**Symptoms:** Frequencies don't sum to 100%, or wild-type estimated at 0%

**Cause:** Overlapping frequency groups

**Solution:** This is expected behavior - tool normalizes frequencies in the end. Check the final report for normalized values.

---

## Citation

If you use this tool in your research, please cite:

```
Mihindukulasuriya KA, Handley SA. VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit.
GitHub: https://github.com/mihinduk/VICAST (2025)
```

And describe the method:
> "Frequency-stratified consensus sequences were generated using VICAST's haplotype consensus module. High-frequency variants (≥95%) were assumed linked based on the pigeonhole principle, which guarantees minimum 90% co-occurrence. Minority variants were reported individually without linkage inference."

---

## Support and Questions

- **GitHub Issues:** https://github.com/mihinduk/VICAST/issues
- **Documentation:** https://github.com/mihinduk/VICAST/wiki
- **Maintainer:** Kathie A. Mihindukulasuriya

---

**Last Updated:** 2026-02-04
**Version:** 2.2.0
