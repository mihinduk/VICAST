# VICAST Testing Guide

## Quick Test Instructions

Copy and paste these commands to test VICAST on HTCF.

---

## Setup

### 1. Navigate to VICAST directory
```bash
cd /scratch/sahlab/kathie/vicast/VICAST
```

### 2. Activate conda (choose one method)

**Method A: Create vicast environment (recommended)**
```bash
source /ref/sahlab/software/anaconda3/bin/activate
conda env create -f environment.yml
conda activate vicast
```

**Method B: Use existing environment with Python/Biopython**
```bash
# If vicast creation fails, try using an existing environment
source /ref/sahlab/software/anaconda3/bin/activate
conda activate viral_genomics  # or py311, or proj_virome
```

**Method C: Test without conda (basic functionality)**
```bash
# Use system Python if available
module load python3  # if available on your cluster
```

---

## Test 1: Pathway Detection (Easiest)

**Purpose:** Check if genome is in SnpEff or needs annotation

**Command:**
```bash
python3 vicast-annotate/step0_check_snpeff.py NC_045512
```

**Expected Output:**
```
STEP 0: Checking Genome Annotation Status
Genome ID: NC_045512
Checking SnpEff database...
  ✓ Found in SnpEff (built-in)
    
RECOMMENDATION
Pathway 1: Already in SnpEff
```

**What it tests:**
- Python script execution
- Pathway detection logic
- SnpEff environment check

**If it works:** ✓ Pathway 1 detection works!

---

## Test 2: Download and Parse Genome (Medium)

**Purpose:** Test standard annotation pipeline (Pathway 2)

**Step 2a: Download genome**
```bash
# Create test directory
mkdir -p test_dengue
cd test_dengue

# Download Dengue virus genome
python3 -c "
from Bio import Entrez, SeqIO
Entrez.email = 'vicast_test@example.com'
handle = Entrez.efetch(db='nucleotide', id='NC_001477', rettype='gb', retmode='text')
with open('NC_001477.gb', 'w') as f:
    f.write(handle.read())
handle.close()

handle = Entrez.efetch(db='nucleotide', id='NC_001477', rettype='fasta', retmode='text')
with open('NC_001477.fasta', 'w') as f:
    f.write(handle.read())
handle.close()
print('Downloaded NC_001477 (Dengue virus)')
"
```

**Step 2b: Parse GenBank file**
```bash
python3 ../vicast-annotate/step1_parse_viral_genome.py NC_001477
```

**Expected Output:**
```
STEP 1: Parse Viral Genome and Create Editable TSV
Parsing GenBank file and skipping polyproteins...
  Skipping polyprotein: polyprotein
  Features written: 10
  Polyproteins skipped: 1

Converting to tab-delimited format for editing...
Created tab-delimited file: NC_001477_no_polyprotein.tsv
  Total features: 10
```

**Expected Files:**
- `NC_001477.gb` - GenBank file
- `NC_001477.fasta` - Genome sequence
- `NC_001477_no_polyprotein.gff3` - GFF3 annotation
- `NC_001477_no_polyprotein.tsv` - Editable TSV

**Verify TSV:**
```bash
head NC_001477_no_polyprotein.tsv
```

**What it tests:**
- Biopython functionality
- NCBI download
- GenBank parsing
- Polyprotein skipping
- TSV generation

**If it works:** ✓ Pathway 2 (Step 1) works!

---

## Test 3: Full Pathway 2 with SnpEff (Harder - requires SnpEff setup)

**Purpose:** Complete annotation pipeline

**Prerequisites:**
- SnpEff installed and configured
- Environment variables set

**Step 3a: Setup SnpEff environment**
```bash
# Check if SnpEff is available
export SNPEFF_HOME=/path/to/snpeff  # Update this path
export SNPEFF_DATA=$SNPEFF_HOME/data
export SNPEFF_JAR=$SNPEFF_HOME/snpEff.jar

# Verify
ls -l $SNPEFF_JAR
```

**Step 3b: Run step2**
```bash
cd test_dengue
python3 ../vicast-annotate/step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv
```

**Expected Output:**
```
STEP 2: Add Edited Viral Genome to snpEff
Converting edited TSV to GFF3...
Validating GFF3 file...
Adding genome to snpEff database...
Building snpEff database for NC_001477...
[SUCCESS] Successfully built snpEff database for NC_001477
```

**What it tests:**
- GFF3 validation
- SnpEff integration
- Database building

**If it works:** ✓ Pathway 2 complete!

---

## Test 4: Segmented Virus (Advanced)

**Purpose:** Test influenza multi-segment handling

**Command:**
```bash
cd /scratch/sahlab/kathie/vicast/VICAST
mkdir -p test_influenza
cd test_influenza

python3 ../vicast-annotate/vicast_annotate_segmented.py test_influenza \
  --segments CY121680,CY121681 \
  --names PB2,PB1
```

**Note:** This tests with just 2 segments (faster). For full influenza, use all 8 segments.

**Expected Output:**
```
VICAST-annotate: Segmented Virus Pipeline
Genome ID: test_influenza
Downloading segments from NCBI...
  Downloading CY121680...
    ✓ Downloaded CY121680
  Downloading CY121681...
    ✓ Downloaded CY121681
Combining 2 FASTA files...
  ✓ Combined 2 sequences
Combining 2 GenBank files to GFF3...
  ✓ Total features: X
```

**What it tests:**
- Multi-segment download
- FASTA combining
- GFF3 combining
- SnpEff multi-chromosome database

---

## Test 5: Master Controller (Full Automation)

**Purpose:** Test automatic pathway selection

**Command:**
```bash
cd /scratch/sahlab/kathie/vicast/VICAST
python3 vicast-annotate/vicast_annotate.py NC_001477
```

**Expected Flow:**
1. Runs step0 (pathway detection)
2. Detects Pathway 2 (well-annotated)
3. Runs step1 automatically
4. Pauses for manual curation
5. Runs step2 after confirmation

**What it tests:**
- Full pipeline integration
- Automatic decision-making

---

## Troubleshooting

### "ModuleNotFoundError: No module named 'Bio'"
**Solution:** Biopython not installed
```bash
conda activate vicast
# or
pip install biopython
```

### "SNPEFF_HOME not set"
**Solution:** Set SnpEff environment variables
```bash
export SNPEFF_HOME=/path/to/snpeff
source setup/snpeff_env.sh
```

### "conda env create fails / gets killed"
**Reason:** Memory limit on login node

**Solution:** Submit as job
```bash
#!/bin/bash
#SBATCH --job-name=create_vicast_env
#SBATCH --mem=16G
#SBATCH --time=2:00:00

source /ref/sahlab/software/anaconda3/bin/activate
cd /scratch/sahlab/kathie/vicast/VICAST
conda env create -f environment.yml
```

### "Download failed"
**Solution:** Check internet connectivity or use manual download
```bash
# Manual download from NCBI
wget "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=NC_001477&rettype=gb&retmode=text" -O NC_001477.gb
```

---

## Quick Success Checklist

- [ ] Test 1: Pathway detection works (NC_045512)
- [ ] Test 2a: Download genome works (NC_001477)
- [ ] Test 2b: Parse GenBank works
- [ ] Test 2c: TSV file created and readable
- [ ] Test 3: SnpEff integration (if SnpEff available)
- [ ] Test 4: Segmented virus handling
- [ ] Test 5: Master controller runs

---

## Minimal Test (No Dependencies)

If conda fails, test core Python functionality:

```bash
cd /scratch/sahlab/kathie/vicast/VICAST

# Test Python can import scripts
python3 -c "
import sys
sys.path.insert(0, 'vicast-annotate')
from step0_check_snpeff import check_snpeff_database
print('✓ Python imports work')
"

# Test file structure
ls -la vicast-annotate/
echo "✓ Files present"

# Test executable permissions
ls -l vicast-annotate/*.py | grep "^-rwx"
echo "✓ Scripts executable"
```

---

## Support

For issues or questions:
- GitHub Issues: https://github.com/mihinduk/VICAST/issues
- GitHub Discussions: https://github.com/mihinduk/VICAST/discussions

---

## Testing Recommendations

**Start with easiest tests first:**

1. **Test 1** - Pathway detection (no setup needed)
2. **Test 2a/2b** - Download and parse (needs Biopython)
3. **Test 2c** - Verify TSV generation
4. **Test 3** - Full pipeline (needs SnpEff)
5. **Test 4** - Segmented viruses (most advanced)
6. **Test 5** - Master controller (full automation)

**Minimum viable test:** Tests 1-2 will validate core functionality
**Impressive demo:** Test 4 with full 8-segment influenza
