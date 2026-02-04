# VICAST Benchmarks

This directory contains benchmarking tools and documentation for comparing VICAST against other viral annotation tools, particularly VADR.

## Overview

VICAST addresses a different use case than tools like VADR. While both annotate viral genomes, they serve complementary purposes in the viral genomics workflow.

## Contents

```
benchmarks/
├── README.md                           # This file
├── vadr_comparison/
│   ├── COMPARISON.md                   # Detailed feature comparison
│   ├── use_cases.md                    # When to use each tool
│   └── benchmark_results.md            # Benchmark results summary
├── scripts/
│   ├── run_benchmark.py                # Main benchmark runner
│   ├── compare_annotations.py          # Compare annotation outputs
│   └── generate_metrics.py             # Calculate comparison metrics
└── results/
    └── (benchmark output files)
```

## Quick Start

```bash
# Run full benchmark suite
python benchmarks/scripts/run_benchmark.py --genome NC_001477 --tools vicast,vadr

# Compare specific annotation outputs
python benchmarks/scripts/compare_annotations.py \
    --vicast output/vicast_annotations.gff3 \
    --vadr output/vadr_annotations.gff3
```

## Key Findings

See [vadr_comparison/COMPARISON.md](vadr_comparison/COMPARISON.md) for detailed analysis.

### Summary

| Feature | VICAST | VADR |
|---------|--------|------|
| **Primary Use Case** | Passage study annotation + SnpEff | GenBank submission validation |
| **Virus Coverage** | Any virus | Specific families (models required) |
| **Annotation Source** | NCBI/BLASTx/manual | HMM models from RefSeqs |
| **Manual Curation** | Built-in checkpoint | Post-processing required |
| **SnpEff Integration** | Native | Requires conversion |
| **Variant Calling** | Integrated pipeline | Not included |
| **Segmented Viruses** | Native support | Model-dependent |

### Complementary Workflow

For comprehensive viral analysis, VICAST and VADR can be used together:

1. **VADR** → Validate genome assembly quality and basic annotation
2. **VICAST** → Curate annotations and integrate with SnpEff for variant analysis

## Citation

If you use these benchmarks, please cite:

```
Mihindukulasuriya KA, Handley SA. VICAST: Viral Cultured-virus Annotation
and SnpEff Toolkit. GitHub: https://github.com/mihinduk/VICAST

Schäffer AA, et al. (2020). VADR: validation and annotation of virus
sequence submissions to GenBank. BMC Bioinformatics. 21:211.
```
