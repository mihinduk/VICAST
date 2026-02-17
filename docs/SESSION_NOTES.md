# VICAST Development Session Notes

## Session: February 2026 - Publication Preparation

### Overview

This session focused on creating real-world validation datasets and improving the tooling for setting up SnpEff databases.

---

## Completed Work

### 1. VCF REF Base Validation

**Problem**: Test VCF files had incorrect REF bases that didn't match the reference genome. SnpEff would annotate them but with warnings, and the annotations could be incorrect.

**Solution**: Added `validate_vcf_ref_bases()` to `src/vicast/validation.py`

```python
from vicast.validation import validate_vcf_ref_bases

is_valid, errors, warnings = validate_vcf_ref_bases(
    vcf_path="variants.vcf",
    fasta_path="reference.fasta",
    max_errors=100  # Stop after 100 errors
)

if not is_valid:
    for error in errors:
        print(f"ERROR: {error}")
```

**Key Implementation Details**:
- Uses `_parse_fasta_sequences()` helper to load full sequences
- Compares VCF REF field against actual bases at that position
- Returns tuple of (is_valid, errors, warnings)
- Handles multi-base REF alleles correctly
- Case-insensitive comparison

---

### 2. SnpEff Setup Tools

#### tools/setup_snpeff.sh

Automated SnpEff installation script:

```bash
# Install SnpEff only
./tools/setup_snpeff.sh --install-only

# Install and build a custom genome
./tools/setup_snpeff.sh --build-genome my_virus genome.fasta genes.gff3
```

#### vicast-annotate/genbank_to_gff3.py

Converts GenBank files to SnpEff-compatible GFF3:

```bash
# Single file
python genbank_to_gff3.py NC_001477.gb -o dengue.gff3

# Multiple files with combined FASTA
python genbank_to_gff3.py segment*.gb -o flu.gff3 --fasta-out flu.fasta --validate
```

#### vicast-annotate/setup_genome.py

One-command SnpEff genome setup:

```bash
# From GenBank
python setup_genome.py dengue_1 --genbank NC_001477.gb --validate

# From GFF3/FASTA
python setup_genome.py my_virus --gff genes.gff3 --fasta genome.fasta
```

---

### 3. Influenza A Validation Dataset

**Location**: `examples/influenza_validation/`

**Purpose**: Validates multi-segment virus handling (8 segments in unified database)

**Files**:
- `NC_02643*.gb` - Individual segment GenBank files
- `NC_02643*.fasta` - Individual segment FASTA files
- `influenza_A_California_2009_8segments.fasta` - Combined FASTA
- `influenza_A_California_2009.gff3` - Combined GFF3
- `test_variants.vcf` - 5 test variants across segments
- `setup_and_validate.sh` - End-to-end validation script

**Test Variants**:
| Segment | Position | Gene | Type |
|---------|----------|------|------|
| NC_026433.1 | 100 | HA | Synonymous |
| NC_026433.1 | 300 | HA | Synonymous |
| NC_026434.1 | 200 | NA | Missense |
| NC_026435.1 | 500 | PB1 | Synonymous |
| NC_026438.1 | 1000 | PB2 | Missense |

---

### 4. SARS-CoV-2 Validation Dataset

**Location**: `examples/sars_cov2_validation/`

**Purpose**: Validates polyprotein annotation (ORF1ab → 16 mature peptides)

**Files**:
- `NC_045512.gb` - GenBank from NCBI
- `NC_045512.fasta` - Reference genome
- `NC_045512.gff3` - GFF3 with mature peptides
- `test_variants.vcf` - 6 clinically significant variants
- `genbank_to_gff3_polyprotein.py` - Custom converter for mature peptides
- `setup_and_validate.sh` - End-to-end validation script

**Test Variants**:
| Position | Mutation | Gene | AA Change | Significance |
|----------|----------|------|-----------|--------------|
| 23403 | A>G | Spike | D614G | Early fitness mutation |
| 23063 | A>T | Spike | N501Y | Receptor binding (Alpha/Omicron) |
| 23012 | G>A | Spike | E484K | Immune escape (Beta/Gamma) |
| 28881 | G>A | N | R203K | Common in variants |
| 10500 | G>A | ORF1ab | G3412D | nsp5 (3CLpro) region |
| 15000 | T>C | ORF1ab | V4912A | nsp12 (RdRp) region |

**Polyprotein Handling**:

The `genbank_to_gff3_polyprotein.py` script:
1. Extracts CDS features (skips ORF1a, keeps ORF1ab with frameshift)
2. Extracts mat_peptide features (nsp1-16)
3. Deduplicates by coordinate (pp1a and pp1ab share same peptides)
4. Maps product names to nsp numbers using regex and lookup table

Key code for nsp name extraction:
```python
# Priority: product > note (note can have misleading "former nsp1" text)
nsp_match = re.match(r'^(nsp\d+)\b', product, re.IGNORECASE)
if nsp_match:
    nsp_name = nsp_match.group(1).lower()

# Fallback to known product mappings
product_to_nsp = {
    'leader protein': 'nsp1',
    '3c-like proteinase': 'nsp5',
    'rna-dependent rna polymerase': 'nsp12',
    'helicase': 'nsp13',
    # ...
}
```

---

## Technical Notes

### SnpEff Database Building

SnpEff requires specific file names:
- `sequences.fa` - Reference genome
- `genes.gff` - GFF3 annotations

Directory structure:
```
$SNPEFF_HOME/data/genome_name/
├── sequences.fa
└── genes.gff
```

Config entry in `snpEff.config`:
```
genome_name.genome : Human-readable description
```

Build command:
```bash
java -jar snpEff.jar build -gff3 -v genome_name
```

### GFF3 Format for SnpEff

Required hierarchy:
```
gene → mRNA → CDS
```

Each feature needs:
- Unique ID attribute
- Parent attribute (except gene)
- gene attribute on CDS/mRNA

Example:
```
chr1  NCBI  gene  100  500  .  +  .  ID=gene_S;Name=S
chr1  NCBI  mRNA  100  500  .  +  .  ID=mRNA_S;Parent=gene_S;gene=S
chr1  NCBI  CDS   100  500  .  +  0  ID=cds_S;Parent=mRNA_S;gene=S;product=spike
```

### Ribosomal Frameshift Handling

SARS-CoV-2 ORF1ab uses ribosomal frameshift at position 13468. In GenBank:
```
CDS  join(266..13468,13468..21555)
```

The overlapping position (13468) is where the frameshift occurs. SnpEff handles this as two CDS parts with the same parent mRNA.

---

## Test Commands

```bash
# All tests (112 total as of this session)
PYTHONPATH=src pytest tests/ -v

# Just validation tests
PYTHONPATH=src pytest tests/test_validation.py -v

# Just conservation tests
PYTHONPATH=src pytest tests/test_conservation.py -v

# Validate a VCF
PYTHONPATH=src python -c "
from vicast.validation import validate_vcf_ref_bases
valid, errs, warns = validate_vcf_ref_bases('test.vcf', 'ref.fasta')
print('PASS' if valid else 'FAIL')
for e in errs: print(f'  {e}')
"

# Run Influenza validation
./examples/influenza_validation/setup_and_validate.sh

# Run SARS-CoV-2 validation
./examples/sars_cov2_validation/setup_and_validate.sh
```

---

## Known Issues / Warnings

1. **SnpEff WARNING_TRANSCRIPT_MULTIPLE_STOP_CODONS**: Expected for ORF1ab due to frameshift representation. Not a real problem.

2. **macOS grep -P**: Perl regex not available. Scripts use `sed` instead.

3. **Python environment**: Tests require `PYTHONPATH=src` since package isn't pip-installed.

---

## Next Steps

1. **Performance benchmarks**: Quantitative runtime and accuracy profiling
2. **Manuscript**: Start Bioinformatics Application Note draft
3. **Additional validation**: HCV or Picornavirus for another polyprotein example
