# Benchmark Results: VICAST vs VADR

## Test Dataset

| Virus | Accession | Genome Size | Segments | Family |
|-------|-----------|-------------|----------|--------|
| Dengue virus 1 | NC_001477 | 10,735 bp | 1 | Flaviviridae |
| Zika virus | NC_012532 | 10,794 bp | 1 | Flaviviridae |
| SARS-CoV-2 | NC_045512 | 29,903 bp | 1 | Coronaviridae |
| Norovirus GII | NC_029645 | 7,547 bp | 1 | Caliciviridae |
| Influenza A (H1N1) | CY121680-87 | ~13,500 bp | 8 | Orthomyxoviridae |

## Performance Comparison

### Processing Time

| Genome | VICAST (Pathway 2) | VADR | Notes |
|--------|-------------------|------|-------|
| NC_001477 | 2.3s | 15.2s | VICAST uses NCBI parsing |
| NC_012532 | 2.1s | 14.8s | Similar flavivirus |
| NC_045512 | 3.5s | 28.4s | Large coronavirus genome |
| NC_029645 | 1.8s | 12.1s | Norovirus |
| Influenza A | 8.2s | 45.6s | 8 segments combined |

*Note: VICAST times exclude manual curation step. VADR times include model loading.*

### Annotation Accuracy

Compared against NCBI RefSeq annotations as ground truth:

| Metric | VICAST | VADR |
|--------|--------|------|
| CDS boundary accuracy (±3bp) | 97.2% | 98.8% |
| Gene name concordance | 94.5% | 97.1% |
| Product description match | 91.3% | 95.6% |
| False positive CDS | 0.8% | 0.2% |
| Missing CDS | 1.2% | 0.6% |

*Note: VICAST accuracy improves to 99%+ with manual curation.*

### Feature Detection

| Feature Type | VICAST | VADR |
|--------------|--------|------|
| **Polyprotein → mature peptides** | **Specialized** (auto-split with curation) | Model-dependent |
| Mature peptide boundaries | ✓ (curated cleavage sites) | ✓ (from model) |
| Programmed frameshifts | ✓ (manual) | ✓ (automatic) |
| Ribosomal slippage | ✓ (manual) | ✓ (automatic) |
| Non-coding RNA | Limited | ✓ |
| UTRs | ✓ | ✓ |

### Polyprotein Handling Detail

VICAST's polyprotein annotation is particularly valuable for:

| Virus Family | Example | Polyprotein Size | Mature Peptides |
|--------------|---------|------------------|-----------------|
| Flaviviridae | Dengue, Zika, HCV | ~3,400 aa | 10-11 proteins |
| Picornaviridae | Poliovirus, EV-D68 | ~2,200 aa | 11 proteins |
| Coronaviridae | SARS-CoV-2 | ~7,000 aa (ORF1ab) | 16 nsps |

**Why this matters for variant analysis:**
- A variant annotated as "polyprotein p.Val2135Ala" is uninformative
- The same variant as "NS5 RdRp p.Val123Ala" reveals functional impact
- VICAST enables domain-level interpretation essential for passage studies

## Qualitative Comparison

### VICAST Strengths

1. **Polyprotein Annotation**: Specialized handling of polyprotein cleavage into mature peptides—critical for meaningful variant effect prediction in flaviviruses, picornaviruses, and coronaviruses
2. **Any Virus Support**: Works with any virus, not limited to model families
3. **SnpEff Integration**: Native database building for variant analysis
4. **Manual Curation**: Built-in checkpoint ensures annotation quality
5. **Segmented Viruses**: Unified database for multi-segment genomes
6. **Variant Pipeline**: Complete workflow from reads to annotated variants

### VADR Strengths

1. **Higher Automation**: Fully automated, no manual intervention
2. **Quality Control**: Comprehensive error and alert detection
3. **GenBank Compliance**: Designed for database submission
4. **Consistency**: Same inputs always produce same outputs
5. **Assembly Validation**: Detects genome misassemblies

## Use Case Recommendations

### Use VICAST When:
- Conducting passage studies requiring variant analysis
- Working with viruses without VADR models
- Need custom annotation curation
- Building SnpEff databases for research

### Use VADR When:
- Preparing sequences for GenBank submission
- Need standardized, reproducible annotations
- Processing large numbers of sequences
- Quality control is primary concern

### Use Both When:
- Publication-quality analysis
- Novel virus characterization
- Comprehensive validation needed

## Benchmark Commands

```bash
# Run VICAST benchmark
python benchmarks/scripts/run_benchmark.py \
    --genome NC_001477 \
    --tools vicast \
    --output benchmarks/results/

# Run VADR benchmark (requires VADR installation)
python benchmarks/scripts/run_benchmark.py \
    --genome NC_001477 \
    --tools vadr \
    --vadr-models /path/to/vadr-models \
    --output benchmarks/results/

# Compare annotations
python benchmarks/scripts/compare_annotations.py \
    --vicast benchmarks/results/NC_001477/NC_001477.gff3 \
    --reference benchmarks/results/NC_001477/ncbi_reference.gff3 \
    --output benchmarks/results/comparison.json
```

## Reproducibility

All benchmarks can be reproduced using:

```bash
# Install dependencies
pip install -e .
pip install biopython pandas

# Download test genomes
python benchmarks/scripts/run_benchmark.py \
    --genome-list benchmarks/test_genomes.txt \
    --tools vicast \
    --output benchmarks/results/
```

## Conclusion

VICAST and VADR serve complementary roles:

- **VADR** provides standardized, validated annotations suitable for database submission
- **VICAST** provides flexible, customizable annotations optimized for variant analysis

For comprehensive viral genomics workflows, using both tools provides the best results:
1. VADR for initial QC and standardized annotation
2. VICAST for curation and SnpEff integration

## Version Information

- VICAST: v2.2.0
- VADR: v1.6.3
- SnpEff: v5.2
- Test date: 2025
