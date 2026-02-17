# Influenza A Multi-Segment Validation Dataset

This dataset validates VICAST's ability to handle segmented viruses using Influenza A/California/07/2009 (H1N1) as a test case.

## Why This Matters

Most viral annotation tools require **separate processing of each genome segment**, making segmented virus analysis cumbersome. VICAST handles all 8 Influenza segments in a **unified database**, enabling:

- Single SnpEff database for entire genome
- Consistent variant annotation across segments
- Proper handling of spliced genes (M2, NEP, PA-X)

## Dataset Contents

| File | Description |
|------|-------------|
| `influenza_A_California_2009_8segments.fasta` | Combined FASTA with all 8 segments |
| `influenza_A_California_2009.gff3` | Unified GFF3 annotation |
| `test_variants.vcf` | Test VCF with variants across multiple segments |
| `NC_02643*.fasta` | Individual segment FASTA files |
| `NC_02643*.gb` | Original GenBank files from NCBI |
| `genbank_to_gff3.py` | Script to convert GenBank to GFF3 |
| `validate_multisegment.py` | Validation script |

## Genome Structure

| Segment | Accession | Gene(s) | Notes |
|---------|-----------|---------|-------|
| 1 | NC_026438.1 | PB2 | Polymerase basic protein 2 |
| 2 | NC_026435.1 | PB1 | Polymerase basic protein 1 |
| 3 | NC_026437.1 | PA, PA-X | PA-X via ribosomal frameshift |
| 4 | NC_026433.1 | HA | Hemagglutinin |
| 5 | NC_026436.1 | NP | Nucleoprotein |
| 6 | NC_026434.1 | NA | Neuraminidase |
| 7 | NC_026431.1 | M1, M2 | M2 is spliced |
| 8 | NC_026432.1 | NS1, NEP | NEP is spliced |

## Spliced Gene Handling

Influenza has three genes produced by splicing or frameshift:

1. **M2** (segment 7): Spliced from M gene, encodes ion channel
2. **NEP** (segment 8): Spliced from NS gene, nuclear export protein
3. **PA-X** (segment 3): Ribosomal frameshift product

VICAST correctly represents these with multiple CDS features per gene.

## Running Validation

```bash
cd examples/influenza_validation
python validate_multisegment.py
```

## Building SnpEff Database

```bash
# 1. Add genome to snpEff.config
echo "influenza_A_Cal09.genome : Influenza A/California/07/2009 (H1N1)" >> $SNPEFF_HOME/snpEff.config

# 2. Create data directory
mkdir -p $SNPEFF_HOME/data/influenza_A_Cal09

# 3. Copy files
cp influenza_A_California_2009_8segments.fasta $SNPEFF_HOME/data/influenza_A_Cal09/sequences.fa
cp influenza_A_California_2009.gff3 $SNPEFF_HOME/data/influenza_A_Cal09/genes.gff

# 4. Build database
java -jar $SNPEFF_HOME/snpEff.jar build -gff3 -v influenza_A_Cal09
```

## Annotating Variants

```bash
# Annotate test variants
java -jar $SNPEFF_HOME/snpEff.jar influenza_A_Cal09 test_variants.vcf > annotated.vcf

# Example output will show:
# - HA variants (NC_026433.1)
# - NA variants (NC_026434.1)
# - PB1 variants (NC_026435.1)
# - PB2 variants (NC_026438.1)
```

## Data Source

Reference sequences downloaded from NCBI RefSeq:
- Strain: A/California/07/2009 (H1N1)
- Accessions: NC_026431.1 - NC_026438.1
- Downloaded: February 2026

## Key Validation Points

1. **Unified Database**: All 8 segments in single files ✓
2. **Spliced Genes**: M2, NEP, PA-X properly annotated ✓
3. **SnpEff Compatible**: GFF3 passes validation ✓
4. **Multi-Segment VCF**: Variants across segments in one file ✓

## Key Features Demonstrated

- **Unified multi-segment database**: All 8 segments in a single SnpEff database
- **Spliced gene annotation**: M2, NEP, PA-X properly represented
- **Native SnpEff integration**: Direct variant annotation without conversion
- **Built-in curation**: Custom database building from GenBank

## Citation

If using this dataset, please cite:
- VICAST: Mihindukulasuriya KA, Handley SA. GitHub: https://github.com/shandley/VICAST
- Reference strain: Garten et al. (2009). Science 325(5937):197-201
