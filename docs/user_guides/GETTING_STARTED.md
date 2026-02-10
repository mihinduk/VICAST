# VICAST Getting Started Guide

**Version:** 2.2.0
**Last Updated:** 2026-02-05

---

## Table of Contents

1. [What is VICAST?](#what-is-vicast)
2. [Installation Options](#installation-options)
3. [Docker Installation](#docker-installation)
4. [Local Installation (Conda)](#local-installation-conda)
5. [HPC Installation](#hpc-installation)
6. [First Analysis Example](#first-analysis-example)
7. [Verification Steps](#verification-steps)
8. [Next Steps](#next-steps)

---

## What is VICAST?

**VICAST** (Viral Cultured-virus Annotation and SnpEff Toolkit) is a comprehensive suite of pipelines for cultured virus genomic analysis, specializing in:

- **Genome annotation curation** - Add poorly annotated viral genomes to SnpEff databases
- **Variant calling** - Detect low-frequency variants in tissue culture passage studies
- **Contamination screening** - De novo assembly and BLAST-based quality control
- **Quasispecies analysis** - Track viral population diversity and evolution

### VICAST Components

| Component | Purpose | Use When |
|-----------|---------|----------|
| **VICAST-annotate** | Genome annotation → SnpEff database | Adding new virus to pipeline |
| **VICAST-analyze** | FASTQ → VCF variant calling | Analyzing passage experiments |
| **viral_diagnostic** | Contamination screening | Quality control before publication |

---

## Installation Options

Choose the installation method that fits your environment:

| Method | Best For | Install Time | Complexity |
|--------|----------|--------------|------------|
| **Docker** | Quick start, reproducibility | 5-10 min | Low |
| **Conda (Local)** | Development, flexibility | 15-30 min | Medium |
| **HPC/SLURM** | Large-scale analysis | 30-60 min | High |

---

## Docker Installation

### Quick Start with Docker

**Prerequisites:**
- Docker installed ([Get Docker](https://docs.docker.com/get-docker/))
- Internet connection

**Installation:**

```bash
# Option 1: Pull pre-built image (recommended)
docker pull ghcr.io/mihinduk/vicast:latest

# Option 2: Build from source
git clone https://github.com/mihinduk/VICAST.git
cd VICAST
docker build -t vicast:latest .
```

**Verify Installation:**

```bash
# Run interactive shell
docker run -it -v $(pwd):/data vicast:latest bash

# Inside container, test VICAST
python -c "from vicast.config import get_config; print('VICAST ready!')"
```

### Running VICAST with Docker

**Mount your data directory:**

```bash
# Linux/Mac
docker run -it -v $(pwd):/data vicast:latest bash

# Windows PowerShell
docker run -it -v ${PWD}:/data vicast:latest bash

# Windows CMD
docker run -it -v %cd%:/data vicast:latest bash
```

**Example: Annotate a genome**

```bash
docker run -v $(pwd):/data vicast:latest \
    python /opt/vicast/vicast-annotate/step1_parse_viral_genome.py NC_001477
```

**Example: Run variant calling**

```bash
docker run -v $(pwd):/data vicast:latest \
    bash /opt/vicast/vicast-analyze/run_vicast_analyze_full.sh \
    sample_R1.fastq.gz sample_R2.fastq.gz NC_001477 4
```

### Docker Environment Variables

The Docker image sets these automatically:

```bash
VICAST_HOME=/opt/vicast
SNPEFF_HOME=/opt/snpEff
SNPEFF_JAR=/opt/snpEff/snpEff.jar
SNPEFF_DATA=/opt/snpEff/data
```

---

## Local Installation (Conda)

### Prerequisites

- **Operating System:** Linux, macOS, or WSL2 on Windows
- **Conda or Mamba:** [Install Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- **Git:** For cloning repository
- **Disk Space:**
  - Annotation only: ~4 GB
  - Full pipeline: ~8-10 GB

### Step 1: Clone Repository

```bash
git clone https://github.com/mihinduk/VICAST.git
cd VICAST
```

### Step 2: Create Conda Environment

Choose the environment for your needs:

**Option A: Annotation Only** (faster, lighter)
```bash
conda env create -f environment_vicast_annotate.yml
conda activate vicast_annotate
```

**Option B: Full Analysis Pipeline** (annotation + variant calling)
```bash
conda env create -f environment_vicast_analyze.yml
conda activate vicast_analyze
```

**Option C: Everything** (complete toolkit)
```bash
conda env create -f environment.yml
conda activate vicast
```

**Environment Comparison:**

| Feature | annotate | analyze | full |
|---------|----------|---------|------|
| Install time | 5-10 min | 10-15 min | 20-30 min |
| Disk space | ~4 GB | ~6 GB | ~8-10 GB |
| Genome annotation | ✅ | ❌ | ✅ |
| Variant calling | ❌ | ✅ | ✅ |

### Step 3: Install VICAST Python Package

```bash
pip install -e .
```

### Step 4: Install SnpEff

**Download SnpEff:**

```bash
cd ~/software  # or your preferred location
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

**Set Environment Variables:**

```bash
# Add to ~/.bashrc or ~/.zshrc
export SNPEFF_HOME="${HOME}/software/snpEff"
export SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
export SNPEFF_DATA="${SNPEFF_HOME}/data"

# Reload shell configuration
source ~/.bashrc  # or source ~/.zshrc
```

**Verify Java version (SnpEff requires Java 21+):**

```bash
java -version
# Should show: openjdk version "21" or higher

# If Java 21 not installed, install via conda:
conda install -c conda-forge openjdk=21
```

### Step 5: Configure VICAST

```bash
# Create configuration directory
mkdir -p ~/.vicast

# Copy configuration template
cp vicast_config.template.sh ~/.vicast/config.sh

# Edit configuration
nano ~/.vicast/config.sh
```

**Example configuration (`~/.vicast/config.sh`):**

```bash
# SnpEff Configuration
export SNPEFF_HOME="${HOME}/software/snpEff"
export SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
export SNPEFF_DATA="${SNPEFF_HOME}/data"

# Java (use conda-provided if available)
export JAVA_HOME="${CONDA_PREFIX}"

# NCBI Configuration
export NCBI_EMAIL="your_email@institution.edu"

# Optional: Performance settings
export VICAST_THREADS=4
export VICAST_SCRATCH="/tmp/vicast"
```

**Load configuration automatically:**

```bash
# Add to ~/.bashrc or ~/.zshrc
if [ -f ~/.vicast/config.sh ]; then
    source ~/.vicast/config.sh
fi
```

### Step 6: Validate Installation

```bash
python -c "from vicast.config import get_config; get_config().print_status()"
```

**Expected output:**

```
VICAST Configuration Status
===========================
✓ SNPEFF_JAR: /home/user/software/snpEff/snpEff.jar
✓ SNPEFF_DATA: /home/user/software/snpEff/data
✓ Python packages: installed
✓ Configuration: valid
```

---

## HPC Installation

### For SLURM/PBS Systems (HTCF, etc.)

**⚠️ IMPORTANT:** Do NOT install on login nodes! Use compute nodes to avoid filling shared `/tmp`.

### Step 1: Clone Repository

```bash
# On login node (safe)
cd /scratch/$USER  # or your working directory
git clone https://github.com/mihinduk/VICAST.git
cd VICAST
```

### Step 2: Create Environment (Submit as Job)

**Option A: Use provided SLURM script**

```bash
# For annotation environment
sbatch create_vicast_env.slurm environment_vicast_annotate.yml

# Or for analysis environment
sbatch create_vicast_env.slurm environment_vicast_analyze.yml

# Monitor job
squeue -u $USER
```

**Option B: Interactive session**

```bash
# Request interactive node
srun --mem=32G --time=2:00:00 --tmp=20G --pty bash

# Activate conda
source /ref/sahlab/software/anaconda3/bin/activate  # adjust for your system

# Create environment
conda env create -f environment_vicast_analyze.yml

# Exit when done
exit
```

### Step 3: Configure for Shared Environment

**For HTCF users (shared lab environment):**

```bash
# Environment is already created at:
# /ref/sahlab/software/envs/vicast_analyze

# No additional setup needed!
```

**For other HPC systems:**

```bash
# Create shared environment location
export CONDA_ENVS_DIRS="/shared/lab/envs:$HOME/.conda/envs"

# Create environment
conda env create -f environment_vicast_analyze.yml
```

### Step 4: Set Up Module System (Optional)

Create a module file for easy loading:

```bash
# Create modulefile: /shared/modulefiles/vicast/2.2.0
#%Module1.0

proc ModulesHelp { } {
    puts stderr "VICAST - Viral Cultured-virus Annotation and SnpEff Toolkit"
}

module-whatis "VICAST viral genomics toolkit"

set root /path/to/VICAST

prepend-path PATH ${root}/vicast-annotate
prepend-path PATH ${root}/vicast-analyze
setenv VICAST_HOME ${root}
setenv SNPEFF_JAR /shared/snpEff/snpEff.jar
setenv SNPEFF_DATA /shared/snpEff/data
```

**Load module:**

```bash
module load vicast/2.2.0
```

### Step 5: Configure Pipeline Paths

```bash
cd vicast-analyze

# Copy template
cp pipeline_config.template.sh pipeline_config.sh

# Edit for your HPC system
nano pipeline_config.sh
```

**Example HPC configuration:**

```bash
# Conda/Mamba command
MAMBA_CMD="conda run -n vicast_analyze"

# SnpEff paths
SNPEFF_DIR="/ref/sahlab/software/snpEff"
SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"

# Java
JAVA_PATH="java"

# Pipeline location
PIPELINE_BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
```

---

## First Analysis Example

Let's run a complete analysis of Dengue virus (DENV) sample data.

### Example Dataset

We'll use a small test dataset (if you don't have your own):

```bash
# Create test directory
mkdir -p ~/vicast_test
cd ~/vicast_test

# Option 1: Use your own data
# R1=your_sample_R1.fastq.gz
# R2=your_sample_R2.fastq.gz
# ACCESSION=NC_001477  # Dengue virus

# Option 2: Download example data (if available)
# wget https://example.com/test_R1.fastq.gz
# wget https://example.com/test_R2.fastq.gz
```

### Step 1: Annotate Genome (One-Time Setup)

**Activate environment:**

```bash
conda activate vicast_annotate  # or vicast_analyze, or vicast
```

**Check if genome is already in SnpEff:**

```bash
cd /path/to/VICAST/vicast-annotate

python step0_check_snpeff.py NC_001477
```

**If not found, annotate it:**

```bash
# Pathway 2: Well-annotated NCBI genome
python step1_parse_viral_genome.py NC_001477

# Output: NC_001477_no_polyprotein.tsv
# Review and edit TSV if needed

# Add to SnpEff
python step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv
```

**Expected output:**

```
✓ Downloaded NC_001477.fasta
✓ Parsed annotations (10 CDS features)
✓ Created NC_001477_no_polyprotein.tsv
✓ Generated GFF3 file
✓ Built SnpEff database
✓ NC_001477 ready for use!
```

### Step 2: Run Variant Calling (QC First)

**Activate analysis environment:**

```bash
conda activate vicast_analyze  # or vicast
```

**Run QC workflow (Steps 1-6):**

```bash
cd /path/to/VICAST/vicast-analyze

./run_vicast_analyze_qc_only.sh \
    /path/to/test_R1.fastq.gz \
    /path/to/test_R2.fastq.gz \
    NC_001477 \
    4  # threads
```

**Check outputs:**

```bash
# Depth coverage
cat sample_results/sample_depth.txt | head

# Diagnostic report
cat diagnostic_sample/sample_diagnostic_report.txt

# Contamination check
cat diagnostic_sample/sample_viral_blast.tsv
```

### Step 3: Review QC and Continue

**If QC looks good, run annotation (Steps 7-9):**

```bash
./run_vicast_analyze_annotate_only.sh \
    /path/to/test_R1.fastq.gz \
    /path/to/test_R2.fastq.gz \
    NC_001477
```

**Or run complete pipeline (all steps):**

```bash
./run_vicast_analyze_full.sh \
    /path/to/test_R1.fastq.gz \
    /path/to/test_R2.fastq.gz \
    NC_001477 \
    4
```

### Step 4: Examine Results

**Check annotated variants:**

```bash
# Annotated VCF
less cleaned_seqs/variants/sample.snpEFF.ann.vcf

# Annotated TSV (easier to read)
less cleaned_seqs/variants/sample.snpEFF.ann.tsv

# Count variants
grep -v "^#" cleaned_seqs/variants/sample.snpEFF.ann.tsv | wc -l
```

**Filter and parse variants:**

```bash
python parse_snpeff_tsv.py \
    cleaned_seqs/variants/sample.snpEFF.ann.tsv \
    sample_results/sample_filtered_mutations.tsv \
    --quality 1000 \
    --depth 200 \
    --freq 0.01
```

**Generate consensus genome:**

```bash
python generate_realistic_haplotype_consensus.py \
    --vcf sample_results/sample_filtered_mutations.tsv \
    --reference cleaned_seqs/NC_001477.fasta \
    --accession NC_001477 \
    --freq 0.50 \
    --output-prefix sample_results/sample_consensus
```

---

## Verification Steps

### Test 1: Configuration Validation

```bash
python -c "from vicast.config import get_config; get_config().print_status()"
```

**Expected:** All checks pass with ✓

### Test 2: SnpEff Database

```bash
java -jar $SNPEFF_JAR databases | grep NC_001477
```

**Expected:** Shows NC_001477 in database list

### Test 3: Conda Environment

```bash
# Check key tools
which bwa
which samtools
which lofreq
which fastp

# Check Python packages
python -c "import Bio; print('Biopython:', Bio.__version__)"
python -c "import pandas; print('Pandas:', pandas.__version__)"
```

**Expected:** All commands found, packages imported

### Test 4: VICAST Commands

```bash
# Test annotation script
python vicast-annotate/step0_check_snpeff.py --help

# Test analysis script
bash vicast-analyze/run_vicast_analyze_full.sh
```

**Expected:** Help messages displayed without errors

### Test 5: End-to-End Pipeline

```bash
# Run test script (if provided)
bash test_vicast_env.sh
```

**Expected:** All tests pass

---

## Next Steps

### Learn the Workflows

1. **[Genome Annotation Guide](VICAST_ANNOTATE_GUIDE.md)** - Add new viruses to SnpEff
2. **[Variant Calling Guide](VICAST_ANALYZE_GUIDE.md)** - Complete analysis pipeline
3. **[Contamination Screening](CONTAMINATION_SCREENING_GUIDE.md)** - Quality control

### Common Workflows

**Adding a New Virus:**
```bash
cd vicast-annotate
python vicast_annotate.py YOUR_ACCESSION
```

**Analyzing Passage Experiment:**
```bash
cd vicast-analyze
for sample in P0 P5 P10; do
    ./run_vicast_analyze_full.sh ${sample}_R1.fq.gz ${sample}_R2.fq.gz NC_001477 4
done
```

**Contamination Screening:**
```bash
./viral_diagnostic.sh sample_R1.fq.gz sample_R2.fq.gz NC_001477 sample_name 4
```

---

## Troubleshooting

### Issue: "conda: command not found"

**Solution:**
```bash
# Initialize conda
source ~/miniconda3/bin/activate  # adjust path
conda init bash  # or zsh
# Restart terminal
```

### Issue: "Java version too old"

**Solution:**
```bash
# Install Java 21 via conda
conda install -c conda-forge openjdk=21
```

### Issue: "SNPEFF_JAR not found"

**Solution:**
```bash
# Set environment variables
export SNPEFF_JAR=/path/to/snpEff.jar
export SNPEFF_DATA=/path/to/snpEff/data

# Make permanent
echo 'export SNPEFF_JAR=/path/to/snpEff.jar' >> ~/.bashrc
source ~/.bashrc
```

### Issue: "Disk quota exceeded" (HPC)

**Solution:**
```bash
# Clean conda cache
conda clean --all

# Use scratch space for environments
export CONDA_ENVS_DIRS="/scratch/$USER/envs:$HOME/.conda/envs"
```

### Issue: "Pipeline fails at Step X"

**Solution:**
1. Check error message in output
2. Verify input files exist: `ls -lh your_file.fastq.gz`
3. Check SnpEff database: `java -jar $SNPEFF_JAR databases | grep YOUR_ACCESSION`
4. See [Troubleshooting Guide](TROUBLESHOOTING.md) for detailed solutions

---

## Getting Help

- **Documentation:** [User Guides](README.md)
- **Issues:** https://github.com/mihinduk/VICAST/issues
- **Email:** Check repository for contact information

---

## Quick Reference

### Essential Commands

```bash
# Activate environment
conda activate vicast_analyze

# Check configuration
python -c "from vicast.config import get_config; get_config().print_status()"

# Annotate genome
python vicast-annotate/step1_parse_viral_genome.py NC_001477

# Run variant calling
bash vicast-analyze/run_vicast_analyze_full.sh R1.fq.gz R2.fq.gz NC_001477 4

# Screen contamination
bash vicast-analyze/viral_diagnostic.sh R1.fq.gz R2.fq.gz NC_001477 sample 4
```

### Key Files and Directories

```
VICAST/
├── vicast-annotate/     # Genome annotation scripts
├── vicast-analyze/      # Variant calling pipeline
├── environment*.yml     # Conda environments
├── vicast_config.template.sh  # Configuration template
└── docs/user_guides/    # Documentation
```

---

**Ready to start?** Continue with:
- **[Genome Annotation Guide](VICAST_ANNOTATE_GUIDE.md)** to add your virus
- **[Variant Calling Guide](VICAST_ANALYZE_GUIDE.md)** to analyze samples

---

**Last Updated:** 2026-02-05
**VICAST Version:** 2.2.0
