# Testing VICAST Environment Creation

This guide walks through creating and testing the `vicast` conda environment from `environment.yml`.

## Important: HPC Best Practices

**⚠️ DO NOT create conda environments on login nodes!**

**Why?**
- Conda uses `/tmp` in your home directory by default
- Login node `/tmp` is small and **shared by all users**
- Installing on login nodes fills this space and causes problems for everyone
- Conda solver also needs substantial memory (16-32 GB)

**Always use:**
- ✅ SLURM batch job (recommended - see below)
- ✅ Interactive session with `srun`
- ❌ NEVER directly on login node

## Step 1: Submit Environment Creation Job

Since conda environment creation can be memory-intensive and uses significant /tmp space, we'll use SLURM:

```bash
# SSH to HTCF
ssh htcf

# Navigate to VICAST directory
cd /scratch/sahlab/kathie/vicast/VICAST

# Submit the job
sbatch create_vicast_env.slurm
```

**Expected output:**
```
Submitted batch job 12345678
```

## Step 2: Monitor Job Progress

```bash
# Check job status
squeue -u mihindu

# Watch the log file in real-time
tail -f create_vicast_env_*.out

# Or check periodically
cat create_vicast_env_*.out
```

**Alternative: Use interactive session**

If you prefer to watch the installation interactively:

```bash
# Request interactive session with resources
srun --mem=32G --time=2:00:00 --tmp=20G --pty bash

# Once in interactive session
source /ref/sahlab/software/anaconda3/bin/activate
conda env create -f environment_minimal.yml

# Exit when done
exit
```

**Job typically takes:** 10-30 minutes

**Success indicators to look for:**
- ✓ "SUCCESS! Environment created"
- ✓ Biopython, Pandas, NumPy versions shown
- ✓ blastx version shown
- ✓ "To use this environment:" instructions appear

## Step 3: Activate and Test Environment

Once the job completes successfully:

```bash
# Activate conda
source /ref/sahlab/software/anaconda3/bin/activate

# Activate vicast environment
conda activate vicast

# Run test suite
./test_vicast_env.sh
```

**Expected output from test script:**
```
==========================================
Testing VICAST Environment
==========================================

✓ Environment: vicast

Test 1: Python Package Imports
-------------------------------
✓ Biopython 1.79
✓ Pandas 1.5.3
✓ NumPy 1.23.5

Test 2: BLAST Installation
---------------------------
✓ blastx found
blastx: 2.15.0+

Test 3: VICAST Script Execution
--------------------------------
Testing step0_check_snpeff.py...
✓ step0_check_snpeff.py loads
Testing step1_parse_viral_genome.py...
✓ step1_parse_viral_genome.py loads
Testing step1_blastx_annotate.py...
✓ step1_blastx_annotate.py loads
Testing step2_add_to_snpeff.py...
✓ step2_add_to_snpeff.py loads

Test 4: Pathway Detection (Quick Functional Test)
--------------------------------------------------
[Shows pathway detection output for NC_045512]

==========================================
All Tests Passed! ✓
==========================================
```

## Step 4: Quick Functional Tests

Try a real VICAST workflow:

### Test Pathway 2 (Well-annotated genome)

```bash
# Create test directory
mkdir -p /scratch/sahlab/kathie/vicast/test/env_validation
cd /scratch/sahlab/kathie/vicast/test/env_validation

# Download and parse Dengue genome (should auto-download)
python3 ../../VICAST/vicast-annotate/step1_parse_viral_genome.py NC_001477

# Check output files
ls -lh NC_001477*
```

**Expected files:**
- `NC_001477.gb` (GenBank file)
- `NC_001477.fasta` (genome sequence)
- `NC_001477_no_polyprotein.gff3` (GFF3 annotation)
- `NC_001477_no_polyprotein.tsv` (editable TSV)

### Test Pathway 3 (BLASTx) - Optional

Only if you have a BLAST database available:

```bash
# Test with small genome and local database
python3 ../../VICAST/vicast-annotate/step1_blastx_annotate.py \
  NC_001477.fasta \
  --blast-db /path/to/blast/db \
  --evalue 1e-5
```

## Troubleshooting

### Job Gets Killed

**Symptom:** `create_vicast_env_*.out` ends with "Killed" or "Out of memory"

**Solution:** Increase memory in SLURM script:
```bash
# Edit create_vicast_env.slurm
#SBATCH --mem=32G  # Increase from 16G to 32G
```

### Conda Solve Takes Forever

**Symptom:** Job runs for >1 hour, stuck on "Solving environment"

**Solutions:**
1. Try mamba (faster solver):
   ```bash
   conda install -n base conda-libmamba-solver
   conda config --set solver libmamba
   ```

2. Or simplify environment.yml (remove version pins)

### Package Import Fails

**Symptom:** `✗ Biopython import failed`

**Solution:** Install missing package manually:
```bash
conda activate vicast
conda install -c bioconda biopython
```

### BLAST Not Found

**Symptom:** `✗ blastx not found`

**Solution:**
```bash
conda activate vicast
conda install -c bioconda blast
```

## Verifying for New Users

If you're testing this to ensure it works for new users, pay attention to:

1. **Installation time:** Should complete in <30 minutes
2. **Disk space:** Check environment size
   ```bash
   du -sh /ref/sahlab/software/anaconda3/envs/vicast
   ```
   Expected: ~5-7 GB

3. **Package conflicts:** No warnings about incompatible versions

4. **All scripts load:** All 4 VICAST scripts pass in Test 3

5. **Functional test works:** Pathway 2 test successfully downloads and parses genome

## Success Criteria

Environment is ready for new users if:
- ✅ Job completes without errors
- ✅ All packages import successfully
- ✅ BLAST 2.15+ is installed
- ✅ All 4 VICAST scripts load without errors
- ✅ Pathway 2 successfully processes a test genome
- ✅ Total disk usage is reasonable (<10 GB)

## Cleanup (if needed)

To remove the test environment:

```bash
conda deactivate
conda env remove -n vicast
```

---

**Note:** Once validated, document any environment.yml changes needed in a commit message.
