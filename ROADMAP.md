# VICAST Development Roadmap

## Overview

This roadmap outlines completed work, items under review, and planned future enhancements for VICAST (Viral Cultured-virus Annotation and SnpEff Toolkit).

---

## Completed (publication-prep branch)

### Infrastructure & Portability
- [x] **Remove hardcoded paths** - Replaced all HTCF-specific paths with portable alternatives
- [x] **Unified configuration system** - Created `vicast-annotate/config.py` for centralized settings
- [x] **Docker/Singularity support** - Added `Dockerfile` and `singularity.def` for containerized deployment
- [x] **Modern Python packaging** - Created `pyproject.toml` with proper dependencies and entry points
- [x] **GitHub Actions CI/CD** - Added `.github/workflows/test.yml` for automated testing

### Testing & Validation
- [x] **pytest test suite for vicast-annotate** - 31 tests covering annotation pipeline
- [x] **pytest test suite for vicast-analyze** - 23 tests covering analysis pipeline
- [x] **Example dataset** - Added NC_001477 (Dengue virus) with reproducible workflow
- [x] **Validation module** - Created `examples/validate_workflow.py` for end-to-end verification

### Code Quality
- [x] **Logging module** - Replaced print statements with `vicast-annotate/logging_config.py`
- [x] **CLI improvements** - Replaced interactive prompts with command-line flags (`--auto`, `--skip-curation`)
- [x] **Script portability** - Updated shell scripts to use `#!/usr/bin/env bash` and portable paths

### Documentation
- [x] **Feature documentation** - Comprehensive comparison highlighting VICAST strengths:
  - Polyprotein annotation (mature peptide splitting)
  - Contamination screening pipeline
  - Segmented virus unified database handling
  - Quasispecies/haplotype reconstruction
  - Complex indel parsing
  - Annotation gap detection & repair

---

## In Review (PR pending)

The `publication-prep` branch contains all completed work above and is ready for review. Key commits:

1. Add portable packaging and CI infrastructure
2. Make Python and shell scripts portable
3. Add comprehensive pytest test suites
4. Add example dataset and validation module
5. Add logging module for standardized output
6. Add feature documentation

---

## Planned Enhancements

### Phase 1: Core Improvements (Next)
- [x] **MSA Conservation Module** - Leverage ColabFold MSAs for conservation scoring
  - Parse existing MSAs from ColabFold/AlphaFold databases (A3M, Stockholm, Clustal, FASTA)
  - Calculate per-position conservation scores (Shannon entropy, percent identity)
  - Annotate variants with conservation impact (7 new TSV columns)
  - CLI script: `vicast-analyze/add_conservation_scores.py`
  - 42 unit tests in `tests/test_conservation.py`

- [x] **Real-world validation dataset** - Influenza A multi-segment validation
  - Complete 8-segment Influenza A/California/07/2009 (H1N1) dataset
  - GenBank to GFF3 converter for SnpEff database building
  - VCF REF base validation to prevent genome mismatch errors
  - End-to-end setup script: `examples/influenza_validation/setup_and_validate.sh`
  - 16 validation tests in `tests/test_validation.py`

- [x] **SARS-CoV-2 polyprotein validation** - Demonstrates polyprotein handling
  - Complete reference genome with ORF1ab polyprotein annotation
  - 16 mature peptides (nsp1-16) in GFF3
  - Test VCF with clinically significant mutations (D614G, N501Y, E484K)
  - End-to-end setup script: `examples/sars_cov2_validation/setup_and_validate.sh`

- [ ] **Additional validation datasets** - Expand test coverage:
  - Polyprotein virus (HCV or Picornavirus)

- [ ] **Integration tests** - End-to-end tests with real sequencing data
  - Simulated reads for controlled testing
  - Known variant detection validation

### Phase 2: Enhanced Analysis
- [ ] **Structural context integration** - Link variants to protein structures
  - Integration with BFVD (Big Fantastic Virus Database)
  - Domain-aware variant annotation
  - Active site proximity scoring

- [ ] **Performance benchmarks** - Quantitative comparisons
  - Runtime profiling
  - Memory usage profiling
  - Scalability testing with large datasets

- [ ] **Batch processing mode** - Improved handling of multiple samples
  - Parallel processing support
  - Progress tracking
  - Summary reports

### Phase 3: Publication
- [ ] **Bioinformatics Application Note** - Manuscript preparation
  - Methods section with workflow diagrams
  - Benchmark results and comparisons
  - Use case examples

- [ ] **Documentation site** - Comprehensive user guide
  - Installation instructions
  - Tutorial with examples
  - API documentation

---

## Contributing

We welcome contributions! To contribute:

1. Fork the repository
2. Create a feature branch from `main`
3. Make your changes with tests
4. Submit a pull request

Please ensure all tests pass before submitting:
```bash
pytest tests/ -v
```

---

## Version History

- **v2.2.0** (current) - Publication preparation release
  - Portable infrastructure
  - Comprehensive testing
  - Docker/Singularity support

- **v2.1.x** - Previous versions (HTCF-specific)

---

## Contact

For questions or collaboration:
- Open an issue on GitHub
- Contact the maintainers directly

## License

See LICENSE file for details.
