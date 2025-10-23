# VICAST Environment Files

VICAST provides three conda environment files for different use cases.

## Quick Recommendation

**For genome annotation (Pathways 1-4):**
```bash
conda env create -f environment_vicast_annotate.yml
```

**For variant calling and passage analysis:**
```bash
conda env create -f environment_vicast_analyze.yml
```

**For everything (annotation + variant calling):**
```bash
conda env create -f environment.yml
```

---

## environment_vicast_annotate.yml ⚡ (MOST COMMON)

**Best for:** Genome annotation workflows

**Includes:**
- Python 3.9
- Biopython, Pandas, NumPy
- BLAST (Pathway 3 - homology annotation)
- EMBOSS (Pathway 2 QC - gap detection/ORF finding)
- SnpEff (variant annotation)

**Advantages:**
- ✅ **Lightweight** (~3-4 GB, fast install)
- ✅ **Purpose-built** for annotation
- ✅ **Less memory** (works with 16GB RAM)
- ✅ **Quick setup** (~5-10 minutes)

**Use this if:**
- Adding viral genomes to SnpEff
- Curating polyprotein annotations
- BLASTx-based annotation transfer
- Only need annotation (no variant calling)

---

## environment_vicast_analyze.yml 📊 (VARIANT CALLING)

**Best for:** Passage studies and variant calling

**Includes everything in annotate, plus:**
- bwa, samtools, bcftools (alignment)
- lofreq, freebayes (variant callers)
- fastp, fastqc, multiqc (QC)
- megahit, spades (assembly)
- bedtools, mosdepth (coverage)
- Additional analysis tools

**Advantages:**
- ✅ **Complete toolkit** for passage experiments
- ✅ **Low-frequency variants** (lofreq)
- ✅ **Coverage analysis**
- ✅ **Multiple variant callers**

**Use this if:**
- Tracking virus evolution across passages
- Detecting minority variants
- Analyzing deep sequencing data
- Need full genomic analysis pipeline

---

## environment.yml 🔬 (EVERYTHING)

**Best for:** Complete VICAST (annotation + variant calling)

**Includes:**
- All vicast_annotate tools
- All vicast_analyze tools
- Everything in one environment

**Advantages:**
- ✅ **One environment** for all workflows
- ✅ **No switching** between environments
- ✅ **Convenient** for diverse projects

**Use this if:**
- You do both annotation and variant calling
- You want simplicity over size
- Disk space isn't a concern
- You prefer one environment for everything

---

## Comparison

| Feature | annotate | analyze | full |
|---------|----------|---------|------|
| Creation time | 5-10 min | 10-15 min | 20-30 min |
| Disk space | ~3-4 GB | ~5-6 GB | ~8-10 GB |
| Memory needed | 16 GB | 16-24 GB | 32 GB |
| Packages | ~15 | ~35 | ~50 |
| Annotation (Pathways 1-4) | ✅ | ❌ | ✅ |
| Variant calling | ❌ | ✅ | ✅ |
| Passage analysis | ❌ | ✅ | ✅ |

---

## Installation Instructions

### On HPC/SLURM systems (like HTCF):

**⚠️ IMPORTANT: DO NOT install on login nodes!**

Conda uses `/tmp` in your home directory by default, which is small and shared with all users. Installing on login nodes can fill this space and cause problems for everyone.

**Option 1: Submit as a SLURM job (RECOMMENDED)**

```bash
# For annotation environment
sbatch create_vicast_env.slurm environment_vicast_annotate.yml

# Or for variant calling environment
sbatch create_vicast_env.slurm environment_vicast_analyze.yml
```

The SLURM script uses the specified environment file.

**Option 2: Use an interactive session**

```bash
# Request interactive session with enough memory and tmp space
srun --mem=32G --time=2:00:00 --tmp=20G --pty bash

# Once in the interactive session, create environment
source /ref/sahlab/software/miniforge3/bin/activate
conda env create -f environment_vicast_annotate.yml

# Exit when done
exit
```

**Why interactive sessions?**
- ✅ Doesn't fill shared `/tmp` on login nodes
- ✅ Has dedicated temporary space
- ✅ Won't affect other users
- ✅ Has enough memory for conda solver
- ✅ Proper HPC citizenship!

### On local machine or with sufficient memory:

**Annotation environment:**
```bash
conda env create -f environment_vicast_annotate.yml
conda activate vicast_annotate
```

**Variant calling environment:**
```bash
conda env create -f environment_vicast_analyze.yml
conda activate vicast_analyze
```

**Full environment (everything):**
```bash
conda env create -f environment.yml
conda activate vicast
```

---

## Troubleshooting

### "Disk quota exceeded"

This means your home directory is full. This is the **most common issue** on HPC systems.

**Quick fix (works 90% of the time):**

```bash
# Clean conda package cache (completely safe)
conda clean --all
```

This removes cached package files (~2-5 GB typically) that have already been installed. **Your existing environments are NOT affected.**

**Check space before/after:**
```bash
# Before cleaning
du -sh ~/miniforge3/pkgs/cache/
# Example output: 1.9G

# Clean
conda clean --all

# After cleaning
du -sh ~/miniforge3/pkgs/cache/
# Example output: 12M  (much smaller!)
```

**If still getting disk quota errors:**

1. **Check overall home directory usage:**
   ```bash
   du -sh ~
   quota -s  # If available on your system
   ```

2. **Use a different location** (like scratch space):
   ```bash
   conda env create -f environment_minimal.yml -p /scratch/your_username/vicast_env
   conda activate /scratch/your_username/vicast_env
   ```

3. **Contact your system administrator** for more quota

**Prevention:**
```bash
# Set conda to use scratch space for package cache
conda config --add pkgs_dirs /scratch/your_username/.conda/pkgs
```

### "/tmp is full" or "No space left on device" (during conda install)

**Problem:** Conda fills `/tmp` in your home directory on login nodes.

**Why this happens:**
- Conda uses `/tmp` for extracting packages during installation
- Login node `/tmp` is small and shared by all users
- Installing on login nodes affects everyone!

**Solution: Use interactive session or SLURM job**

```bash
# Option 1: Interactive session (get your own /tmp space)
srun --mem=32G --time=2:00:00 --tmp=20G --pty bash
conda env create -f environment_minimal.yml
exit

# Option 2: Submit as job (recommended)
sbatch create_vicast_env.slurm
```

**Or configure conda to use a different tmp directory:**
```bash
# Create tmp directory in scratch space
mkdir -p /scratch/your_username/tmp

# Tell conda to use it
export TMPDIR=/scratch/your_username/tmp
conda env create -f environment_minimal.yml

# Make permanent by adding to ~/.bashrc
echo "export TMPDIR=/scratch/your_username/tmp" >> ~/.bashrc
```

**Note:** This is good HPC citizenship! Don't fill shared resources on login nodes.

### "Solving environment takes forever"

Try these in order:

1. **Use environment_minimal.yml** (faster to solve)

2. **Increase memory allocation** in SLURM script:
   ```bash
   #SBATCH --mem=32G  # or even 64G
   ```

3. **Install mamba** (much faster solver):
   ```bash
   conda install -n base mamba
   ```

### "Package conflicts"

If you get conflicts, try:

1. **Update conda:**
   ```bash
   conda update -n base conda
   ```

2. **Start fresh:**
   ```bash
   conda env remove -n vicast
   conda clean --all
   conda env create -f environment_minimal.yml
   ```

---

## Required: SnpEff Environment Variables

**⚠️ IMPORTANT: Set these before using VICAST!**

VICAST requires SnpEff environment variables to be configured. Without these, step2 will fail.

### On HTCF (typical setup):

```bash
# Add to your ~/.bashrc for permanent setup
echo 'export SNPEFF_JAR=/ref/sahlab/software/snpEff/snpEff.jar' >> ~/.bashrc
echo 'export SNPEFF_DATA=/ref/sahlab/software/snpEff/data' >> ~/.bashrc

# Reload your bashrc
source ~/.bashrc

# Verify
echo $SNPEFF_JAR
java -jar $SNPEFF_JAR -version
```

### On other systems:

```bash
# Replace paths with your SnpEff installation location
export SNPEFF_JAR=/path/to/snpEff.jar
export SNPEFF_DATA=/path/to/snpEff/data

# Make permanent by adding to ~/.bashrc
echo 'export SNPEFF_JAR=/path/to/snpEff.jar' >> ~/.bashrc
echo 'export SNPEFF_DATA=/path/to/snpEff/data' >> ~/.bashrc
```

**Without these variables:**
- ✅ Pathway 1 works (just checks existing databases)
- ✅ Step1 of Pathway 2 works (parsing GenBank)
- ❌ Step2 of Pathway 2 **FAILS** (can't build database)
- ❌ Pathways 3-4 **FAIL** at step2

## Using Existing Environments

**You don't need to create a new environment if you already have one with:**
- Python 3.8+
- Biopython
- Pandas
- (Optional) BLAST for Pathway 3

Just activate your existing environment and use VICAST:
```bash
conda activate your_existing_env

# Make sure SnpEff variables are set (see above!)
echo $SNPEFF_JAR  # Should show path to snpEff.jar

# Then use VICAST
python3 vicast-annotate/step1_parse_viral_genome.py NC_001477
```

---

## Which Environment Should I Use?

```
What's your primary workflow?
│
├─ Genome annotation (Pathways 1-4)
│  └→ environment_vicast_annotate.yml ⚡
│
├─ Variant calling & passage analysis
│  └→ environment_vicast_analyze.yml 📊
│
└─ Both annotation + variant calling
   └→ environment.yml 🔬
```

**Quick decision tree:**
- Adding viruses to SnpEff? → **vicast_annotate**
- Analyzing passage experiments? → **vicast_analyze**
- Doing everything? → **full (environment.yml)**
- Want fastest install? → **vicast_annotate**

---

## After Installation

Verify your installation:
```bash
conda activate vicast
./test_vicast_env.sh
```

All tests should pass with green checkmarks (✓).

---

**Questions?** See [TESTING_VICAST_ENV.md](TESTING_VICAST_ENV.md) or open an issue on GitHub.
