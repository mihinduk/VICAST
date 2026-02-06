# Read-Level Co-occurrence Analysis Guide

## Overview

The `check_read_cooccurrence.py` module enables **true haplotype reconstruction** by checking if variants actually appear together on the same sequencing reads. This provides **direct evidence** of variant linkage, moving beyond statistical inference.

## Scientific Background

### The Problem: Frequency-Based Inference vs. Proven Linkage

**Current approach** (`generate_realistic_haplotype_consensus.py`):
- Uses **frequency-based inference** to group variants
- Assumes high-frequency variants (≥95%) co-occur based on the pigeonhole principle
- **Limitation**: Cannot PROVE that mutations actually appear on the same viral genomes

**This solution** (`check_read_cooccurrence.py`):
- Examines actual sequencing reads to check variant co-occurrence
- Provides **direct evidence** of linkage for nearby variants
- Enables claims like: "mutations co-occur on 450/1000 reads (45% proven)"

### When Does This Work?

Read-level phasing works for variants within sequencing distance:

| Variant Distance | Phasing Method | Typical Success |
|-----------------|----------------|-----------------|
| < 150bp | Single read | ✓ Excellent |
| 150-500bp | Read pairs (insert size) | ✓ Good |
| 500-1000bp | Long reads or pairs | △ Limited |
| > 1000bp | Not possible with short reads | ✗ No data |

For **typical Illumina paired-end sequencing**:
- Read length: ~150bp
- Insert size: ~300-500bp
- **Optimal range**: Variants < 500bp apart

## Installation

### Dependencies

```bash
# Install via conda (recommended)
conda install -c bioconda pysam
conda install pandas

# Or via pip
pip install pysam pandas
```

### Verify Installation

```bash
python3 -c "import pysam; print(f'pysam version: {pysam.__version__}')"
python3 -c "import pandas; print(f'pandas version: {pandas.__version__}')"
```

## Usage

### Basic Usage

```bash
# Analyze all variant pairs in VCF file
python check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf
```

This will:
1. Parse variants from VCF
2. Find all variant pairs within 500bp (default)
3. Check which reads contain both variants
4. Output co-occurrence statistics to `cooccurrence.tsv`

### Common Use Cases

#### 1. High-Quality Variants Only

```bash
# Only analyze high-confidence variants
python check_read_cooccurrence.py \
    --bam aligned.bam \
    --vcf variants.vcf \
    --min-qual 1000 \
    --min-depth 200 \
    --min-freq 0.05
```

#### 2. Custom Distance Threshold

```bash
# Analyze only very close variants (within read length)
python check_read_cooccurrence.py \
    --bam aligned.bam \
    --vcf variants.vcf \
    --max-distance 150
```

#### 3. Using VICAST Filtered Output

```bash
# VICAST outputs filtered TSV files - these work directly
python check_read_cooccurrence.py \
    --bam sample_aligned_sorted.bam \
    --vcf sample_annotated_filtered.tsv \
    --output sample_cooccurrence.tsv
```

#### 4. Single Reads Only (No Pairs)

```bash
# Don't use read pair information (stricter, but more conservative)
python check_read_cooccurrence.py \
    --bam aligned.bam \
    --vcf variants.vcf \
    --no-pairs
```

## Interpreting Results

### Output Format

The output TSV file contains:

| Column | Description | Example |
|--------|-------------|---------|
| variant1_pos | Position of first variant (1-based) | 1000 |
| variant2_pos | Position of second variant (1-based) | 1234 |
| distance | Distance between variants (bp) | 234 |
| evidence_level | Type of evidence (read-level/insert-level) | read-level |
| total_spanning_reads | Reads covering both positions | 1000 |
| informative_reads | Reads with determinable alleles | 950 |
| both_alt | Reads with both alternate alleles | 450 |
| variant1_only | Reads with only first variant | 450 |
| variant2_only | Reads with only second variant | 25 |
| both_ref | Reads with both reference alleles | 25 |
| undetermined | Reads with unclear alleles | 50 |
| cooccurrence_rate | Proportion with both variants | 0.4737 |
| linkage_proven | Whether co-occurrence was observed | True |

### Interpreting Co-occurrence Rates

**High co-occurrence rate (>80%)**:
- Variants likely on the same haplotype
- Strong evidence for linkage
- Example: Two mutations that arose together and are maintained

**Medium co-occurrence rate (40-80%)**:
- Variants may co-occur in some genomes but not all
- Possible evidence of multiple haplotypes
- Example: One mutation at 90%, another at 50%, they co-occur in the 50% population

**Low co-occurrence rate (<40%)**:
- Variants are mostly independent
- May be on different haplotypes
- Example: Two different sub-populations with different mutations

**No co-occurrence (0%)**:
- Variants are mutually exclusive
- Strong evidence for separate haplotypes
- Example: Two alternative mutations at the same locus

### Evidence Levels

**read-level** (distance < 150bp):
- Both variants on the same physical read
- Highest confidence
- No possibility of artifact from read pairing

**insert-level** (distance 150-500bp):
- Variants on read pair (same DNA fragment)
- High confidence
- Requires proper read pairing

**fragment-level** (distance > 500bp):
- Limited or no evidence possible
- Depends on long reads or special libraries

---

## Using Co-occurrence Data for Manual Variant Review

### Overview

Co-occurrence analysis provides **evidence-based decision support** for Manual Checkpoint 2 in the VICAST workflow - the variant filtering review step. Instead of relying solely on frequency thresholds and quality scores, you can use read-level evidence to distinguish real variants from artifacts.

### Workflow Integration

#### Standard VICAST Workflow (Without BAM)

```bash
# Step 7-9: Variant calling and filtering
./run_vicast_analyze_annotate_only.sh sample_R1.fq.gz sample_R2.fq.gz NC_001477.1

# MANUAL CHECKPOINT 2: Review variants
# - Check frequencies, quality scores
# - Remove suspected artifacts
# - Decisions based on thresholds alone
```

#### Enhanced Workflow (With BAM Co-occurrence)

```bash
# Step 7-9: Variant calling and filtering
./run_vicast_analyze_annotate_only.sh sample_R1.fq.gz sample_R2.fq.gz NC_001477.1

# NEW: Generate co-occurrence evidence
python ../scripts/check_read_cooccurrence.py \
    --bam cleaned_seqs/mapping/sample.lofreq.realign.bam \
    --vcf cleaned_seqs/variants/sample_vars.filt.vcf \
    --output sample_cooccurrence.tsv

# MANUAL CHECKPOINT 2: Enhanced variant review
# - Review variants AND their linkage patterns
# - Use read-level evidence for decisions
# - Distinguish real haplotypes from artifacts
```

---

### Decision Rules for Variant Filtering

#### Rule 1: Validate High-Frequency Variant Clusters

**Pattern to look for:**
- Multiple high-frequency variants (≥80%) in nearby positions
- Co-occurrence rates ≥95% between them

**Example:**
```
Position 22673  Freq=95%  T>C  │ Co-occurrence Matrix:
Position 22674  Freq=92%  C>T  │ 22673+22674: 99.74% (69,436/69,616)
Position 22679  Freq=93%  T>C  │ 22673+22679: 99.96% (68,670/68,697)
Position 22686  Freq=90%  C>T  │ 22673+22686: 99.97% (66,150/66,169)
                               │ 22674+22679: 99.95% (68,755/68,791)
                               │ 22679+22686: 99.99% (68,512/68,521)
```

**Interpretation:** These variants form a coherent haplotype representing ~90-95% of the viral population.

**Decision:** ✅ **KEEP ALL** - Strong evidence of real linked variants on same haplotype

---

#### Rule 2: Flag Suspicious Mid-Frequency Singletons

**Pattern to look for:**
- Variant at 30-70% frequency
- NO high co-occurrence (≥95%) with any nearby variants
- Doesn't fit into expected haplotype structure

**Example:**
```
Position 5000   Freq=45%  G>A  │ Co-occurrence with nearby variants:
Position 4950   Freq=85%  C>T  │ 5000+4950: 12% (5/42)
Position 5100   Freq=88%  A>G  │ 5000+5100: 8% (3/37)
                               │
Interpretation: 45% variant doesn't belong to major haplotype
```

**Possible causes:**
- Sequencing artifact in difficult region
- Alignment error (homopolymer, repeat)
- PCR error amplified during library prep
- Real recombinant (rare for most viruses)

**Decision:** 🚩 **FLAG FOR INVESTIGATION**

**Checks to perform:**
```bash
# 1. Check coverage depth at position
samtools depth -r NC_001477.1:4990-5010 sample.bam

# 2. Examine alignment quality
samtools view sample.bam NC_001477.1:4990-5010 | less

# 3. Check for nearby homopolymers or repeats
# If region is AAAAAAAA or similar → likely artifact
```

**Final decision:**
- If region has alignment issues: ❌ **REMOVE**
- If alignment looks clean: ⚠️ **MARK AS UNCERTAIN**, report conservatively

---

#### Rule 3: Resolve Ambiguous Low-Frequency Variants

**Pattern A - Linked to major haplotype:**
```
Position 10000  Freq=8%   C>T  │ Co-occurrence:
Position 10050  Freq=85%  A>G  │ 10000+10050: 92% (45/49 reads)
                               │
Interpretation: 8% variant appears WITH the 85% major variant
→ Emerging mutation on the dominant haplotype
```

**Decision:** ✅ **KEEP** - Real emerging variant, biologically plausible (e.g., passage adaptation)

---

**Pattern B - Isolated singleton:**
```
Position 15000  Freq=6%   T>C  │ Co-occurrence:
Position 14950  Freq=87%  G>A  │ 15000+14950: 2% (1/48)
Position 15100  Freq=89%  C>T  │ 15000+15100: 3% (1/35)
                               │
Interpretation: 6% variant is independent of major haplotype
→ Likely sequencing error or rare artifact
```

**Decision:** ❌ **REMOVE** - Insufficient evidence for biological reality

**Threshold guideline:**
- Low-frequency variants (<10%) with <10% co-occurrence with major haplotype → Remove
- Low-frequency variants with >80% co-occurrence with major haplotype → Keep

---

#### Rule 4: Identify Distinct Haplotypes (Mutually Exclusive Variants)

**Pattern to look for:**
- Two variants at similar mid-range frequencies (30-70%)
- ZERO or very low co-occurrence (<5%) between them
- Each variant shows high co-occurrence with DIFFERENT sets of variants

**Example:**
```
Position 10162  Freq=45%  G>A  │ Haplotype A marker
Position 10163  Freq=55%  C>T  │ Haplotype B marker
                               │
Co-occurrence: 10162+10163: 0% (0/79 reads)
→ Mutually exclusive = different haplotypes

Downstream cluster (10179, 10200, 10212, 10230):
10162 + cluster: 77-80% co-occurrence → Haplotype A
10163 + cluster: 19-23% co-occurrence → Haplotype B (inverse)
```

**Interpretation:**
- **Haplotype A**: Carries 10162 variant + downstream cluster (45% frequency)
- **Haplotype B**: Carries 10163 variant alone (55% frequency)

**Decision:** ✅ **KEEP BOTH** - Clear evidence of two distinct viral subpopulations

**Publication language:**
> "Read-level co-occurrence analysis revealed two distinct haplotypes: Haplotype A (45%) carrying position 10162 variant with downstream mutations at 10179-10230 (77-80% linkage), and Haplotype B (55%) defined by the mutually exclusive 10163 variant (0% co-occurrence with 10162)."

---

### Variant Review Checklist Template

Use this checklist during Manual Checkpoint 2:

```
┌─────────────────────────────────────────────────────────────────────────┐
│ VARIANT REVIEW CHECKLIST - BAM Co-occurrence Evidence                  │
├──────────┬──────┬──────┬───────────────────────┬──────────┬────────────┤
│ Position │ Freq │ Type │ Co-occurrence Pattern │ Decision │ Notes      │
├──────────┼──────┼──────┼───────────────────────┼──────────┼────────────┤
│ 22673    │ 95%  │ SNP  │ 99% with 22674-22686  │ KEEP     │ Major hap  │
│ 22674    │ 92%  │ SNP  │ 99% with cluster      │ KEEP     │ Major hap  │
│ 5000     │ 45%  │ SNP  │ No linkage found      │ FLAG     │ Check aln  │
│ 10162    │ 45%  │ SNP  │ 77% with cluster      │ KEEP     │ Hap A      │
│ 10163    │ 55%  │ SNP  │ Inverse to 10162      │ KEEP     │ Hap B      │
│ 15000    │ 6%   │ SNP  │ Isolated singleton    │ REMOVE   │ Artifact   │
└──────────┴──────┴──────┴───────────────────────┴──────────┴────────────┘

Decision Key:
  KEEP     - Strong evidence, retain in final VCF
  REMOVE   - Likely artifact, exclude from analysis
  FLAG     - Requires manual inspection
  UNCERTAIN - Report but note low confidence
```

**Download template:** See `variant_review_template.tsv` in VICAST repository

---

### Step-by-Step Review Process

#### Step 1: Load Both Files

```bash
# Open co-occurrence results
less sample_cooccurrence.tsv

# Open filtered VCF for reference
less cleaned_seqs/variants/sample_vars.filt.vcf
```

#### Step 2: Identify High-Frequency Clusters

```bash
# Extract high-frequency variants (≥80%)
awk -F'\t' 'NR==1 || $NF >= 0.80' cleaned_seqs/variants/sample_vars.filt.vcf > high_freq_vars.vcf

# Find which pairs have >95% co-occurrence
awk -F'\t' 'NR==1 || $12 >= 0.95' sample_cooccurrence.tsv > high_cooccur.tsv
```

**Expected result:** Clusters of nearby variants with high co-occurrence

#### Step 3: Flag Suspicious Singletons

Look for variants that:
- Have moderate frequency (30-70%)
- Don't appear in `high_cooccur.tsv`
- Aren't part of any cluster

Mark these for manual inspection.

#### Step 4: Classify Low-Frequency Variants

```bash
# Extract low-frequency variants (5-20%)
awk -F'\t' 'NR==1 || ($NF >= 0.05 && $NF < 0.20)' sample_vars.filt.vcf > low_freq_vars.vcf
```

For each low-frequency variant:
- Check if it co-occurs (>80%) with any high-frequency variant
  - YES → Keep (emerging mutation)
  - NO → Remove (likely artifact)

#### Step 5: Create Filtered VCF

```bash
# Create final filtered VCF with decisions applied
python ../scripts/apply_cooccurrence_filters.py \
    --vcf cleaned_seqs/variants/sample_vars.filt.vcf \
    --cooccurrence sample_cooccurrence.tsv \
    --review-checklist variant_review.tsv \
    --output sample_vars.final.vcf
```

---

### Common Patterns and Interpretations

#### Pattern 1: Perfect Linkage Cluster (Publication-Quality)

```
Positions: 1000, 1025, 1050, 1080
Frequencies: 95%, 94%, 96%, 93%
Co-occurrence: All pairs >99.5%

Interpretation: Single dominant haplotype at ~95% frequency
Confidence: ★★★★★ Excellent
Publication: "Variants at positions 1000-1080 showed >99.5% co-occurrence
              (n=10,000-12,000 reads per pair), confirming linkage on a
              single dominant haplotype."
```

#### Pattern 2: Two Distinct Haplotypes

```
Haplotype A markers: 2000, 2050  (45% each, 99% co-occurrence)
Haplotype B markers: 2005, 2055  (55% each, 98% co-occurrence)
A vs B co-occurrence: <5%

Interpretation: Two competing viral populations
Confidence: ★★★★☆ High
Publication: "BAM analysis revealed two distinct haplotypes with mutually
              exclusive markers (A: 45%, B: 55%; <5% co-occurrence between
              lineages, n=500-800 reads)."
```

#### Pattern 3: Emerging Variant on Major Haplotype

```
Major haplotype: Positions 5000-5200 at 90% frequency
New variant: Position 5100 at 12% frequency
Co-occurrence: 5100 with major haplotype: 85%

Interpretation: Emerging mutation arising on dominant background
Confidence: ★★★☆☆ Moderate
Publication: "Position 5100 variant (12%) showed 85% co-occurrence with
              the major haplotype (90%), suggesting emergence through
              mutation rather than co-infection."
```

#### Pattern 4: Suspicious Artifact (Remove)

```
Position 8000: 35% frequency
Co-occurrence with all nearby variants: <10%
Region: Homopolymer (AAAAAAA)
Quality: Mean Q=25

Interpretation: Likely sequencing artifact in difficult region
Confidence: ★☆☆☆☆ Low (variant should be removed)
Publication: Do not report
```

---

### Quality Control Metrics

Track these metrics during variant review:

```
Review Session Metrics:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
Total variants in filtered VCF: 244
High-frequency (≥80%): 45
  ├─ In high-cooccur clusters: 42 (93%) ✓ Good
  └─ Singletons: 3 (7%) → FLAG

Mid-frequency (20-80%): 28
  ├─ Linked to haplotypes: 20 (71%) ✓ Good
  └─ Orphaned: 8 (29%) → INVESTIGATE

Low-frequency (<20%): 171
  ├─ Linked to major haplotype: 12 (7%) → KEEP
  └─ Singletons: 159 (93%) → REMOVE

Variants with co-occurrence data: 89/244 (36%)
Variants too distant (>500bp): 155 (64%) → Use frequency alone

Final retained variants: 74
  ├─ With read-level evidence: 62 (84%)
  └─ Frequency-based only: 12 (16%)
```

---

### Integration with VICAST-Analyze Workflow

#### Recommended Insertion Point

Add BAM co-occurrence as **Step 8B** in the VICAST-Analyze pipeline:

**Current workflow:**
```
Step 7: Variant calling (lofreq)
Step 8: Automated filtering
Step 9: SnpEff annotation
```

**Enhanced workflow:**
```
Step 7: Variant calling (lofreq)
Step 8A: Automated filtering
Step 8B: BAM co-occurrence analysis (NEW - optional)
        MANUAL CHECKPOINT 2: Review with co-occurrence evidence
Step 9: SnpEff annotation
```

#### Command-Line Integration

Add to `run_vicast_analyze_annotate_only.sh`:

```bash
# After Step 8A: Filtering
echo "Step 8B: Generating co-occurrence evidence (optional)..."
if [ -f "${PIPELINE_BASE}/scripts/check_read_cooccurrence.py" ]; then
    python "${PIPELINE_BASE}/scripts/check_read_cooccurrence.py" \
        --bam "${OUTPUT_DIR}/mapping/${SAMPLE}.lofreq.realign.bam" \
        --vcf "${OUTPUT_DIR}/variants/${SAMPLE}_vars.filt.vcf" \
        --output "${OUTPUT_DIR}/variants/${SAMPLE}_cooccurrence.tsv"

    echo "Co-occurrence results: ${OUTPUT_DIR}/variants/${SAMPLE}_cooccurrence.tsv"
    echo "MANUAL REVIEW: Check variants with co-occurrence evidence before proceeding"
    echo "Press Enter to continue to SnpEff annotation, or Ctrl+C to review first..."
    read
fi
```

---

## Integration with Haplotype Consensus

### Workflow

1. **Run variant calling** (VICAST pipeline)
   ```bash
   bash run_vicast_analyze_full.sh sample
   ```

2. **Check co-occurrence** (this tool)
   ```bash
   python check_read_cooccurrence.py \
       --bam sample_aligned_sorted.bam \
       --vcf sample_annotated_filtered.tsv \
       --output sample_cooccurrence.tsv
   ```

3. **Generate frequency-based haplotypes**
   ```bash
   python generate_realistic_haplotype_consensus.py \
       --vcf sample_annotated_filtered.tsv \
       --reference reference.fasta \
       --accession MN908947.3 \
       --output-prefix sample
   ```

4. **Compare and validate**
   - Check if frequency-based predictions match read-level evidence
   - Update haplotype descriptions with proven linkage information
   - Identify any discrepancies that might indicate errors

### Example Validation

**Frequency-based prediction:**
```
Dominant haplotype: 85% frequency
  Mutations: S:D614G, S:N501Y (both at 85% frequency)
  Inference: Likely co-occur based on matching frequencies
```

**Read-level validation:**
```
variant_pair: 23403A>G + 23063A>T (S:D614G + S:N501Y)
distance: 340bp
both_alt: 780/900 reads (86.7%)
linkage_proven: TRUE
```

**Conclusion:** ✓ Frequency prediction VALIDATED by read-level evidence

## Advanced Usage

### Programmatic Use

```python
from check_read_cooccurrence import VariantCooccurrenceAnalyzer, parse_vcf

# Load variants
variants = parse_vcf("variants.vcf", min_qual=1000, min_freq=0.05)

# Initialize analyzer
analyzer = VariantCooccurrenceAnalyzer(
    "aligned.bam",
    min_base_quality=20,
    min_mapping_quality=20
)

# Analyze specific variant pair
result = analyzer.analyze_variant_pair(
    chrom="MN908947.3",
    variant1={'pos': 1000, 'ref': 'A', 'alt': 'T', 'id': '1001A>T'},
    variant2={'pos': 1200, 'ref': 'C', 'alt': 'G', 'id': '1201C>G'},
    use_pairs=True
)

print(f"Co-occurrence rate: {result['cooccurrence_rate']:.2%}")
print(f"Linkage proven: {result['linkage_proven']}")

# Analyze all pairs
results = analyzer.analyze_all_variant_pairs(
    variants,
    max_distance=500,
    use_pairs=True
)
```

### Custom Filtering

```python
# Filter results for high-confidence linkage
proven_links = [r for r in results if r['linkage_proven']]
high_cooccurrence = [r for r in proven_links if r['cooccurrence_rate'] > 0.8]
strong_evidence = [r for r in high_cooccurrence if r['informative_reads'] > 100]

for result in strong_evidence:
    print(f"{result['variant_pair']}: {result['cooccurrence_rate']:.2%}")
```

## Troubleshooting

### BAM Index Not Found

**Error:** `BAM index not found at aligned.bam.bai`

**Solution:**
```bash
samtools index aligned.bam
```

The script will attempt to create the index automatically, but you may need to do it manually if there are permission issues.

### No Spanning Reads Found

**Error:** `No spanning reads found`

**Possible causes:**
1. Variants are too far apart (> insert size)
2. Low coverage at variant positions
3. BAM file doesn't match VCF reference

**Solutions:**
- Check variant positions: `samtools depth aligned.bam | grep -E "1000|1200"`
- Verify coverage is adequate: `samtools coverage aligned.bam`
- Confirm BAM and VCF reference match: `samtools view -H aligned.bam`

### Low Co-occurrence Rate

**Observation:** Co-occurrence rate much lower than variant frequencies

**Possible explanations:**
1. **Different haplotypes**: Variants are on separate viral populations
2. **Sequencing bias**: Technical artifact affecting variant calling
3. **Recombination**: Evidence of viral recombination event
4. **Mapping errors**: Reads incorrectly mapped

**Investigation:**
```bash
# Check individual variant frequencies
samtools mpileup -r test_ref:1000-1001 aligned.bam
samtools mpileup -r test_ref:1200-1201 aligned.bam

# Visualize in IGV (Integrative Genomics Viewer)
# Look for patterns suggesting separate haplotypes
```

### Memory Issues

**Error:** `MemoryError` or `Killed`

**Solution:**
- Process smaller regions at a time
- Increase available memory
- Use `--max-distance` to limit analysis

```bash
# Analyze only very close variants
python check_read_cooccurrence.py \
    --bam aligned.bam \
    --vcf variants.vcf \
    --max-distance 200
```

## Performance Considerations

### Speed

- **Fast**: Analyzing 10-20 variants with 100x coverage: ~1-5 seconds
- **Moderate**: Analyzing 50-100 variants with 500x coverage: ~30-60 seconds
- **Slow**: Analyzing 500+ variants with 1000x+ coverage: several minutes

### Optimization Tips

1. **Pre-filter variants**: Only analyze high-quality variants
2. **Reduce max distance**: Use `--max-distance 200` for read-level only
3. **Use indexed BAM**: Ensure BAM is sorted and indexed
4. **Filter by frequency**: Use `--min-freq 0.05` to skip rare variants

## Biological Interpretation

### Strong Linkage (>90% co-occurrence)

**Interpretation**: Mutations arose together and are maintained as a unit

**Examples**:
- Compensatory mutations (one mutation buffers cost of another)
- Mutations in same functional domain
- Founder effect (single ancestral virus with both mutations)

### Partial Linkage (50-90% co-occurrence)

**Interpretation**: Multiple viral sub-populations with different mutation patterns

**Examples**:
- Mixed infection with different strains
- Within-host evolution producing variants
- Population bottleneck followed by expansion

### No Linkage (<50% co-occurrence)

**Interpretation**: Independent mutations on different haplotypes

**Examples**:
- Mutations at different time points
- Different selective pressures
- Recombination between variants

### Mutually Exclusive (0% co-occurrence)

**Interpretation**: Mutations cannot coexist

**Examples**:
- Alternative mutations at same position
- Lethal combination
- Strong negative epistasis

## Citation

If you use this tool in your research, please cite:

```
VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit
Handley Lab, Washington University in St. Louis
https://github.com/handley-lab/VICAST
```

## Further Reading

### Haplotype Phasing

- **Read-based phasing**: Direct observation of variants on reads
- **Statistical phasing**: Using population genetics to infer haplotypes
- **Long-read sequencing**: PacBio/Nanopore for direct long-range phasing

### Viral Quasispecies

- **Within-host diversity**: Multiple viral variants in single host
- **Haplotype reconstruction**: Determining which mutations co-occur
- **Linkage disequilibrium**: Statistical association between variants

### Related Tools

- **WhatsHap**: Read-based phasing for diploid organisms
- **HapCUT2**: Haplotype assembly from sequencing data
- **PacBio/Nanopore**: Long-read technologies for direct phasing

## Contact

For questions, issues, or feature requests:
- GitHub Issues: https://github.com/handley-lab/VICAST/issues
- Lab website: [Add website]
- Email: [Add contact]

## License

[Specify license - typically MIT or GPL for academic tools]

## Changelog

### Version 1.0.0 (2026-02-05)
- Initial release
- Basic co-occurrence analysis
- VCF/TSV input support
- Read and read-pair analysis
- Comprehensive documentation

### Future Enhancements
- [ ] Indel handling improvements
- [ ] Multi-sample analysis
- [ ] Visualization output (HTML report)
- [ ] Integration with generate_realistic_haplotype_consensus.py
- [ ] Support for long-read data (PacBio/Nanopore)
- [ ] Statistical significance testing
- [ ] Haplotype network visualization
