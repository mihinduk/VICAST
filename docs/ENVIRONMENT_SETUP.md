# VICAST Environment Setup Guide

## Overview

VICAST uses a **single conda environment** that includes all necessary tools for viral genome annotation and SnpEff database creation.

## Installation

### Create VICAST Environment

```bash
# Create environment from YAML file
conda env create -f environment.yml

# Activate environment
conda activate vicast

# Verify installation
python3 -c "import Bio; print('Biopython OK')"
snpeff -version
blastp -version
```

**This environment supports:**
- Pathway 1: Checking existing SnpEff databases
- Pathway 2: Well-annotated genomes (standard pipeline)
- Pathway 3: BLASTx-based annotation
- Pathway 4: Segmented virus handling

## What's Included

The `vicast` environment includes:

| Tool | Version | Purpose |
|------|---------|---------|
| Python | 3.8+ | Core scripting |
| Biopython | Latest | GenBank/FASTA parsing |
| Pandas | Latest | Data manipulation |
| SnpEff | 5.2+ | Variant annotation |
| BLAST+ | 2.15+ | Homology annotation |
| pip | Latest | Python package management |

## Usage Workflow

### Standard Usage

```bash
# Activate environment
conda activate vicast

# Pathway 1: Check if genome exists in SnpEff
python3 vicast-annotate/step0_check_snpeff.py NC_001477

# Pathway 2: Standard annotation (well-annotated genomes)
python3 vicast-annotate/step1_parse_viral_genome.py NC_001477
python3 vicast-annotate/step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv

# Pathway 3: BLASTx annotation (poorly annotated genomes)
python3 vicast-annotate/step1_blastx_annotate.py genome.fasta --blast-db nr
python3 vicast-annotate/step2_add_to_snpeff.py genome genome_blastx.tsv

# Pathway 4: Segmented viruses (multi-chromosome)
python3 vicast-annotate/vicast_annotate_segmented.py influenza_h1n1 \
  --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
  --names PB2,PB1,PA,HA,NP,NA,M,NS
```

## Troubleshooting

### Environment Creation Fails

**Common causes:**

1. **Memory limit**: Conda solver needs substantial memory
   - Solution: Run on login node or request more memory

2. **Channel conflicts**: Package versions incompatible
   - Solution: Update conda channels:
     ```bash
     conda config --add channels defaults
     conda config --add channels bioconda
     conda config --add channels conda-forge
     ```

3. **Disk space**: Environment can be large (~5-7 GB)
   - Solution: Check available space: `df -h`
   - Clean old environments: `conda env remove -n old_env`

### Package Import Errors

**Problem:** `ImportError: No module named 'Bio'`

**Solution:** Ensure environment is activated:
```bash
conda activate vicast
python3 -c "import Bio; print('Success!')"
```

### BLAST Database Issues

**Problem:** BLASTx cannot find database

**Solution:** Set BLASTDB environment variable:
```bash
export BLASTDB=/path/to/blast/databases
# Or specify full path in command:
python3 step1_blastx_annotate.py genome.fasta --blast-db /full/path/to/nr
```

### Which Environment Am I In?

```bash
# Check active environment
conda env list

# Shows active environment with *:
#   base                     /path/to/conda
#   vicast                *  /path/to/conda/envs/vicast
```

## Storage Requirements

| Component | Size |
|-----------|------|
| VICAST environment | ~5-7 GB |
| SnpEff data directory | Varies by genomes added |
| BLAST databases (if using Pathway 3) | Can be large (nr ~200 GB) |

## For HPC/SLURM Users

### Environment Setup on HPC

Create environment on login node (one-time):

```bash
# On HTCF, activate conda first
source /ref/sahlab/software/anaconda3/bin/activate

# Create VICAST environment
conda env create -f environment.yml
```

### Using in SLURM Jobs

```bash
#!/bin/bash
#SBATCH --job-name=vicast_annotate
#SBATCH --mem=8G
#SBATCH --time=2:00:00

# Activate conda (HTCF-specific)
source /ref/sahlab/software/anaconda3/bin/activate

# Activate VICAST environment
conda activate vicast

# Run VICAST pipeline
python3 vicast-annotate/step1_parse_viral_genome.py NC_001477
python3 vicast-annotate/step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv
```

## Best Practices

1. **Keep environment updated**
   ```bash
   conda activate vicast
   conda update --all
   ```

2. **Export working environment** (for reproducibility)
   ```bash
   conda env export > vicast_working_$(date +%Y%m%d).yml
   ```

3. **Clean package cache** (save disk space)
   ```bash
   conda clean --all
   ```

4. **Use environment.yml for installation** (not manual package installs)
   - Ensures all dependencies are compatible
   - Makes environment reproducible

5. **Test after installation**
   ```bash
   conda activate vicast
   validate_vicast_setup  # If available
   ```

## Updating VICAST

When updating VICAST from GitHub:

```bash
# Pull latest changes
cd VICAST
git pull origin main

# Update environment if environment.yml changed
conda env update -f environment.yml --prune
```

## Environment Variables

VICAST uses these environment variables (set via `setup/snpeff_env.sh`):

```bash
export SNPEFF_JAR=/path/to/snpEff.jar
export SNPEFF_DATA=/path/to/snpEff/data
export SCRATCH_DIR=/path/to/scratch  # Optional
```

Load them with:
```bash
source setup/snpeff_env.sh
```

Or add to your `~/.bashrc` for automatic loading:
```bash
# Add to ~/.bashrc
source /path/to/VICAST/setup/snpeff_env.sh
```

## Quick Reference

```bash
# One-time setup
conda env create -f environment.yml

# Each session
conda activate vicast
source setup/snpeff_env.sh

# Run VICAST
python3 vicast-annotate/step0_check_snpeff.py GENOME_ID
python3 vicast-annotate/step1_parse_viral_genome.py GENOME_ID
python3 vicast-annotate/step2_add_to_snpeff.py GENOME_ID GENOME_ID_no_polyprotein.tsv

# When done
conda deactivate
```

---

For dependency details, see: [DEPENDENCIES.md](DEPENDENCIES.md)
