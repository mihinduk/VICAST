# VICAST Environment Files

VICAST provides two conda environment files to suit different needs.

## Quick Recommendation

**For most users (VICAST-annotate only):**
```bash
conda env create -f environment_minimal.yml
```

**For full pipeline (annotation + variant calling):**
```bash
conda env create -f environment.yml
```

---

## environment_minimal.yml ⚡ (RECOMMENDED)

**Best for:** VICAST-annotate pipelines (Pathways 1-4)

**Includes:**
- Python 3.9
- Biopython, Pandas, NumPy
- BLAST (for Pathway 3 - BLASTx annotation)
- SnpEff (for variant annotation)

**Advantages:**
- ✅ **Faster to create** (~5-10 minutes vs 20-30 minutes)
- ✅ **Less memory needed** (works with 16GB)
- ✅ **Smaller disk footprint** (~3-4 GB vs ~7-8 GB)
- ✅ **Fewer conflicts** (fewer packages = easier to solve)
- ✅ **Everything you need** for annotation workflows

**Use this if:**
- You only need VICAST-annotate (genome annotation)
- You want faster environment creation
- You have limited disk space
- You're testing/learning VICAST

---

## environment.yml 🔬 (FULL)

**Best for:** Complete VICAST (annotation + variant calling)

**Includes everything in minimal, plus:**
- bwa, samtools, bcftools (alignment and variant calling)
- megahit (assembly)
- lofreq (low-frequency variant calling)
- fastp, seqkit (read processing)
- Additional analysis tools

**Advantages:**
- ✅ **Complete toolkit** for viral genomics
- ✅ **Ready for VICAST-analyze** (coming soon)
- ✅ **All-in-one environment**

**Use this if:**
- You need variant calling pipelines
- You want to analyze passage experiments
- You're doing full viral genomics workflows
- You have plenty of disk space and time

---

## Comparison

| Feature | minimal | full |
|---------|---------|------|
| Creation time | 5-10 min | 20-30 min |
| Disk space | ~3-4 GB | ~7-8 GB |
| Memory needed | 16 GB | 32 GB |
| Packages | ~30 | ~50 |
| VICAST-annotate | ✅ | ✅ |
| VICAST-analyze | ❌ | ✅ |

---

## Installation Instructions

### On HPC/SLURM systems (like HTCF):

**⚠️ IMPORTANT: DO NOT install on login nodes!**

Conda uses `/tmp` in your home directory by default, which is small and shared with all users. Installing on login nodes can fill this space and cause problems for everyone.

**Option 1: Submit as a SLURM job (RECOMMENDED)**

```bash
# This automatically runs in an interactive session with proper resources
sbatch create_vicast_env.slurm
```

The SLURM script automatically uses `environment_minimal.yml` and falls back to `environment.yml` if needed.

**Option 2: Use an interactive session**

```bash
# Request interactive session with enough memory and tmp space
srun --mem=32G --time=2:00:00 --tmp=20G --pty bash

# Once in the interactive session, create environment
source /ref/sahlab/software/anaconda3/bin/activate
conda env create -f environment_minimal.yml

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

**Minimal environment:**
```bash
conda env create -f environment_minimal.yml
conda activate vicast
```

**Full environment:**
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

## Using Existing Environments

**You don't need to create a new environment if you already have one with:**
- Python 3.8+
- Biopython
- Pandas
- (Optional) BLAST for Pathway 3

Just activate your existing environment and use VICAST:
```bash
conda activate your_existing_env
python3 vicast-annotate/step1_parse_viral_genome.py NC_001477
```

---

## Which Environment Should I Use?

```
Do you need variant calling pipelines?
│
├─ YES → Use environment.yml (full)
│
└─ NO → Use environment_minimal.yml (recommended)
    │
    └─ Just genome annotation? → minimal ✓
```

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
