# VICAST Environment Setup Guide

## Why Two Environments?

VICAST uses **two separate conda environments** due to BLAST version requirements:

- **Pathways 1, 2, 4** work with BLAST 2.13.0
- **Pathway 3 (VADR)** requires BLAST >=2.15.0

VADR has strict dependency requirements that conflict with some tools in the base environment.

## Environment Overview

| Environment | Use For | Key Tools |
|-------------|---------|-----------|
| `vicast` | Pathways 1, 2, 4 | SnpEff, BLAST 2.13, Python tools |
| `vicast-vadr` | Pathway 3 only | VADR, BLAST 2.15+, SnpEff |

## Installation

### Base Environment (vicast)

For most users - supports Pathways 1, 2, and 4:

```bash
# Create environment
conda env create -f environment.yml

# Activate
conda activate vicast

# Verify
python3 -c "import Bio; print('Biopython OK')"
snpeff -version
blast+ -version
```

**Use this for:**
- Pathway 1: Checking existing SnpEff databases
- Pathway 2: Well-annotated genomes (standard pipeline)
- Pathway 4: BLASTx-based annotation
- Segmented virus handling

### VADR Environment (vicast-vadr)

Only needed if using Pathway 3 (VADR-enhanced annotation):

```bash
# Create environment
conda env create -f environment_vadr.yml

# Activate
conda activate vicast-vadr

# Verify VADR
v-annotate.pl -h
python3 -c "import Bio; print('Biopython OK')"
```

**Use this for:**
- Pathway 3: VADR-enhanced annotation
- When using `--use-vadr` flag

## Usage Workflow

### For Pathways 1, 2, 4 (Base Environment)

```bash
# Activate base environment
conda activate vicast

# Pathway 1: Check if genome exists
python3 vicast-annotate/step0_check_snpeff.py NC_001477

# Pathway 2: Standard annotation
python3 vicast-annotate/vicast_annotate.py NC_001477

# Pathway 4: BLASTx annotation
python3 vicast-annotate/step1_blastx_annotate.py genome.fasta --blast-db nr

# Segmented viruses
python3 vicast-annotate/vicast_annotate_segmented.py influenza_h1n1 \
  --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
  --names PB2,PB1,PA,HA,NP,NA,M,NS
```

### For Pathway 3 (VADR Environment)

```bash
# Activate VADR environment
conda activate vicast-vadr

# Run with VADR
python3 vicast-annotate/step1_parse_viral_genome.py NC_009942 --use-vadr

# Continue with step2 in same environment
python3 vicast-annotate/step2_add_to_snpeff.py NC_009942 NC_009942_vadr_curated.tsv
```

## Troubleshooting

### "VADR not found" Error

**Problem:** Running Pathway 3 in base `vicast` environment

**Solution:** Switch to `vicast-vadr` environment:
```bash
conda deactivate
conda activate vicast-vadr
```

### "Wrong BLAST version" Warning

**Problem:** Using wrong environment for pathway

**Solution:**
- Pathways 1, 2, 4: Use `vicast` environment
- Pathway 3: Use `vicast-vadr` environment

### Environment Creation Fails

**Common causes:**
1. **Memory limit**: Conda solver needs substantial memory
   - Solution: Run on login node or request more memory
   
2. **Channel conflicts**: Package versions incompatible
   - Solution: Update conda channels: `conda config --add channels bioconda`
   
3. **Disk space**: Environments can be large (~5-7 GB each)
   - Solution: Check `df -h` and clean old environments

### Which Environment Am I In?

```bash
# Check active environment
conda env list

# Shows active environment with *:
#   base                     /path/to/conda
#   vicast                *  /path/to/conda/envs/vicast
#   vicast-vadr             /path/to/conda/envs/vicast-vadr
```

## Storage Requirements

| Component | Size |
|-----------|------|
| Base `vicast` environment | ~5 GB |
| VADR `vicast-vadr` environment | ~7 GB |
| **Total for both environments** | **~12 GB** |

## Alternative: Single Environment with Updated BLAST

If disk space is limited, you can create a single environment with BLAST 2.15+:

```bash
# Modify environment.yml
# Change: blast=2.13.0
# To: blast>=2.15.0

# Add VADR
# Add line: - vadr>=1.6.4

# Create unified environment
conda env create -f environment.yml
```

**Pros:** One environment, simpler
**Cons:** Larger install, longer solve time, potential conflicts

## Best Practices

1. **Start with base environment** - Most users only need Pathways 1 & 2
2. **Add VADR environment only if needed** - When you need Pathway 3
3. **Keep environments updated** - `conda update --all`
4. **Clean unused environments** - `conda env remove -n old_env`
5. **Export working environments** - `conda env export > working_env.yml`

## Quick Reference

```bash
# Create both environments
conda env create -f environment.yml
conda env create -f environment_vadr.yml

# Use base environment (most common)
conda activate vicast
python3 vicast-annotate/vicast_annotate.py NC_001477

# Use VADR environment (Pathway 3 only)
conda activate vicast-vadr
python3 vicast-annotate/step1_parse_viral_genome.py NC_009942 --use-vadr

# Switch environments
conda deactivate
conda activate vicast  # or vicast-vadr
```

## For HPC/SLURM Users

Create environments on login node, then use in jobs:

```bash
#!/bin/bash
#SBATCH --job-name=vicast
#SBATCH --mem=8G

# Activate appropriate environment
source /path/to/conda/bin/activate
conda activate vicast  # or vicast-vadr

# Run VICAST
python3 vicast-annotate/vicast_annotate.py NC_001477
```

---

For more details, see: [DEPENDENCIES.md](DEPENDENCIES.md)
