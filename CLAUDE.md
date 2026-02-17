# VICAST Claude Code Session Guide

This document provides context for Claude Code sessions working on VICAST.

## Project Overview

**VICAST** (Viral Cultured-virus Annotation and SnpEff Toolkit) is a viral genome annotation and variant analysis pipeline. It's being prepared for publication as a Bioinformatics Application Note.

### Repository Structure

```
vicast/
├── src/vicast/              # Core Python library
│   ├── __init__.py
│   ├── config.py            # Configuration management
│   ├── validation.py        # GFF3/VCF validation (including REF base checking)
│   └── conservation/        # MSA conservation scoring module
│       ├── models.py        # Data classes
│       ├── msa_parser.py    # Parse A3M, Stockholm, Clustal, FASTA
│       ├── scores.py        # Shannon entropy, percent identity
│       ├── mapper.py        # HGVSp -> MSA position mapping
│       └── annotator.py     # Add conservation to TSV
├── vicast-annotate/         # Annotation pipeline scripts
│   ├── genbank_to_gff3.py   # GenBank to GFF3 converter
│   └── setup_genome.py      # SnpEff genome setup tool
├── vicast-analyze/          # Analysis pipeline scripts
│   └── add_conservation_scores.py  # Conservation annotation CLI
├── tools/
│   └── setup_snpeff.sh      # SnpEff installation script
├── tests/                   # pytest test suites
│   ├── test_validation.py   # 16 tests for validation module
│   ├── test_conservation.py # 42 tests for conservation module
│   └── ...
├── examples/
│   ├── influenza_validation/  # 8-segment Influenza A dataset
│   ├── sars_cov2_validation/  # Polyprotein validation dataset
│   └── data/                  # Dengue test data
├── ROADMAP.md               # Development roadmap
└── pyproject.toml           # Python packaging
```

## Current Branch

**`publication-prep`** - All publication preparation work is on this branch.

## Recent Session Work (February 2026)

### Completed This Session

1. **VCF REF Base Validation** (`src/vicast/validation.py`)
   - `validate_vcf_ref_bases()` function checks VCF REF alleles match reference genome
   - Prevents silent annotation errors from SnpEff
   - 16 tests in `tests/test_validation.py`

2. **CLI Tools** (`vicast-annotate/`, `tools/`)
   - `tools/setup_snpeff.sh` - Automated SnpEff installation
   - `vicast-annotate/genbank_to_gff3.py` - GenBank to GFF3 converter
   - `vicast-annotate/setup_genome.py` - One-command SnpEff genome setup

3. **SARS-CoV-2 Validation Dataset** (`examples/sars_cov2_validation/`)
   - Demonstrates polyprotein annotation (ORF1ab → 16 nsps)
   - Test VCF with D614G, N501Y, E484K, R203K mutations
   - Custom `genbank_to_gff3_polyprotein.py` for mature peptide handling

4. **Influenza Validation Dataset** (`examples/influenza_validation/`)
   - Demonstrates multi-segment handling (8 segments)
   - VCF REF bases corrected and validated

### Previous Session Work

1. **MSA Conservation Module** (`src/vicast/conservation/`)
   - Parses MSA files (A3M, Stockholm, Clustal, FASTA)
   - Calculates conservation scores (Shannon entropy, percent identity)
   - Maps HGVSp protein positions to MSA columns
   - 42 tests in `tests/test_conservation.py`

## Remaining Tasks (from ROADMAP.md)

### Phase 1 (Near-term)
- [ ] Additional validation dataset (HCV or Picornavirus polyprotein)
- [ ] Integration tests with simulated reads

### Phase 2 (Medium-term)
- [ ] Structural context integration (BFVD)
- [ ] Performance benchmarks
- [ ] Batch processing mode

### Phase 3 (Publication)
- [ ] Bioinformatics Application Note manuscript
- [ ] Documentation site

## Key Commands

```bash
# Run all tests
PYTHONPATH=src pytest tests/ -v

# Run specific test file
PYTHONPATH=src pytest tests/test_validation.py -v

# Run Influenza validation
./examples/influenza_validation/setup_and_validate.sh

# Run SARS-CoV-2 validation
./examples/sars_cov2_validation/setup_and_validate.sh

# Set up a new SnpEff genome
PYTHONPATH=src python vicast-annotate/setup_genome.py GENOME_NAME --genbank input.gb

# Validate VCF REF bases
PYTHONPATH=src python -c "
from vicast.validation import validate_vcf_ref_bases
is_valid, errors, warnings = validate_vcf_ref_bases('test.vcf', 'ref.fasta')
print('Valid' if is_valid else 'Invalid')
"
```

## Git Status

- Branch: `publication-prep`
- Remote: `origin/publication-prep` (GitHub)
- PR: #1 (open, for publication preparation)

### Recent Commits
```
50cfa45 Add SARS-CoV-2 polyprotein validation dataset
26d238c Add CLI tools for SnpEff setup and genome database building
7a5bdd1 Update ROADMAP.md with Influenza validation dataset completion
2925afb Add VCF REF base validation to prevent genome mismatch errors
```

## Dependencies

- Python 3.8+
- BioPython (for GenBank/FASTA parsing)
- pandas (for TSV handling)
- SnpEff (external, installed to `tools/snpEff/`)

## Test Data Locations

| Dataset | Path | Purpose |
|---------|------|---------|
| Dengue | `examples/data/` | Basic workflow testing |
| Influenza A | `examples/influenza_validation/` | Multi-segment handling |
| SARS-CoV-2 | `examples/sars_cov2_validation/` | Polyprotein annotation |
| Test MSA | `examples/data/test_virus_msa.a3m` | Conservation module testing |

## Important Notes

1. **SnpEff Installation**: Located at `tools/snpEff/` (gitignored). Run `tools/setup_snpeff.sh --install-only` to install.

2. **VCF REF Validation**: Always validate VCF files before annotation:
   ```python
   from vicast.validation import validate_vcf_ref_bases
   is_valid, errors, warnings = validate_vcf_ref_bases(vcf_path, fasta_path)
   ```

3. **Running Tests**: Must use `PYTHONPATH=src` since package isn't installed:
   ```bash
   PYTHONPATH=src pytest tests/ -v
   ```

4. **Polyprotein GFF3**: SARS-CoV-2 uses custom converter (`genbank_to_gff3_polyprotein.py`) to properly handle mature peptides.

## Session Resumption Checklist

When resuming work:

1. Check git status: `git status`
2. Check current branch: `git branch`
3. Review ROADMAP.md for next tasks
4. Run tests to verify state: `PYTHONPATH=src pytest tests/ -v`
5. Check for uncommitted work

## Contact

- Repository: https://github.com/shandley/VICAST
- Authors: Mihindukulasuriya KA, Handley SA
