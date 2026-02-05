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
