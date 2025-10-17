# VICAST Dependencies and Requirements

## Core Requirements (All Pathways)

### Conda Environment
Install the base environment:
```bash
conda env create -f environment.yml
conda activate vicast
```

**Included in environment.yml:**
- Python 3.9
- Biopython 1.79
- Pandas 1.5.3
- NumPy 1.23.5
- Matplotlib 3.6.2
- PyYAML 6.0
- Samtools 1.17
- BCFtools 1.17
- SeqKit 2.3.1
- SnpEff 5.1
- Entrez-direct (optional, for NCBI downloads)

## Pathway-Specific Requirements

### Pathway 1: Already in SnpEff
**No additional requirements** - uses existing SnpEff installation

**Software Needed:**
- SnpEff (included in conda env)
- Java (for SnpEff)

---

### Pathway 2: Well-Annotated (NCBI)
**Standard annotation pipeline**

**Software Needed:**
- Python 3.9+ (included)
- Biopython (included)
- Pandas (included)
- SnpEff (included)

**Optional:**
- Entrez-direct for automatic NCBI downloads

**No additional installations required beyond conda environment**

---

### Pathway 3: VADR-Enhanced Annotation
**Requires VADR installation**

**Software Needed:**
- All Pathway 2 requirements
- **VADR 1.6+** (must be installed separately)
- Perl 5.12+
- Various Perl modules (installed by VADR installer)

**Installation:**
```bash
# Install VADR using provided script
bash setup/install_vadr.sh /path/to/software

# Set environment variable
export VADR_DIR=/path/to/software/vadr
```

**VADR Models:**
The installation script downloads models for:
- Flaviviruses (flavi)
- Coronaviruses (corona)
- Caliciviruses (calici)
- Influenza A (flua)
- And others

**Size Warning:**
- VADR installation: ~2 GB
- Model databases: ~5 GB

---

### Pathway 4: BLASTx-Based Annotation
**Requires BLAST+ and protein databases**

**Software Needed:**
- All Pathway 2 requirements
- **BLAST+ 2.13.0** (included in conda env)

**BLAST Databases:**
Option 1 - NCBI nr (comprehensive but large):
```bash
# Download NCBI nr database (~100 GB)
mkdir -p /path/to/blast_db
cd /path/to/blast_db
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nr.*.tar.gz"
for file in nr.*.tar.gz; do tar -xzf $file; done
```

Option 2 - RefSeq Protein (smaller):
```bash
# Download RefSeq protein database (~40 GB)
mkdir -p /path/to/blast_db
cd /path/to/blast_db
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein.*.tar.gz"
for file in refseq_protein.*.tar.gz; do tar -xzf $file; done
```

Option 3 - Custom Viral Database (recommended):
```bash
# Create custom viral protein database
# Download viral proteins from NCBI
# Format as BLAST database

# Example:
makeblastdb -in viral_proteins.faa -dbtype prot -out viral_db
```

**Set BLAST database path:**
```bash
export BLASTDB=/path/to/blast_db
```

---

## SnpEff Configuration

### Standard Setup
```bash
# Set environment variables
export SOFTWARE_DIR=/path/to/software
source setup/snpeff_env.sh
```

### Advanced: SnpEff 5.2+ with Java 21
For latest SnpEff features:
```bash
# Install Java 21
bash setup/install_java21_conda.sh

# Update SnpEff to 5.2+
# (Follow SnpEff documentation)
```

---

## Environment Variables

### Required for all pathways:
```bash
export SNPEFF_JAR=/path/to/snpEff.jar
export SNPEFF_DATA=/path/to/snpEff/data
export SNPEFF_HOME=/path/to/snpEff
```

### Required for Pathway 3:
```bash
export VADR_DIR=/path/to/vadr
```

### Required for Pathway 4:
```bash
export BLASTDB=/path/to/blast_db
```

### Optional:
```bash
export SCRATCH_DIR=/path/to/scratch
```

---

## Validation

Test your installation:
```bash
# Activate environment
conda activate vicast

# Source SnpEff environment
source setup/snpeff_env.sh

# Run validation
python3 -c "from vicast_validation import validate_vicast_setup; validate_vicast_setup()"
```

This will check:
- ✓ Python packages installed
- ✓ SnpEff accessible
- ✓ VADR installed (if using Pathway 3)
- ✓ BLAST+ installed (if using Pathway 4)
- ✓ Environment variables set

---

## Storage Requirements

### Minimum (Pathways 1-2):
- Conda environment: ~3 GB
- SnpEff databases: ~5 GB
- Working space: ~1 GB per genome
- **Total: ~10 GB**

### With VADR (Pathway 3):
- Add: ~7 GB for VADR and models
- **Total: ~17 GB**

### With BLAST databases (Pathway 4):
- Option 1 (nr): Add ~100 GB
- Option 2 (refseq_protein): Add ~40 GB
- Option 3 (custom viral): Add ~1-5 GB
- **Total: 20-110 GB depending on database choice**

---

## Recommended Setup

For a complete VICAST installation supporting all pathways:

```bash
# 1. Create conda environment
conda env create -f environment.yml
conda activate vicast

# 2. Setup SnpEff
export SOFTWARE_DIR=/path/to/software
source setup/setup_snpeff_custom_paths.sh $SOFTWARE_DIR

# 3. Install VADR (for Pathway 3)
bash setup/install_vadr.sh $SOFTWARE_DIR

# 4. Setup BLAST database (for Pathway 4)
# Choose one option:
# a) Download nr (comprehensive but large)
# b) Download refseq_protein (balanced)
# c) Create custom viral database (recommended)

# 5. Validate installation
python3 -c "from vicast_validation import validate_vicast_setup; validate_vicast_setup()"
```

---

## Troubleshooting

### "VADR not found"
```bash
# Check VADR_DIR
echo $VADR_DIR

# Reinstall if needed
bash setup/install_vadr.sh $SOFTWARE_DIR
```

### "BLAST not found"
```bash
# Check conda environment
conda list | grep blast

# Reinstall if needed
conda install -c bioconda blast
```

### "SnpEff errors"
```bash
# Check Java version
java -version  # Need Java 8+ (Java 21 for SnpEff 5.2+)

# Check SnpEff paths
echo $SNPEFF_JAR
ls -lh $SNPEFF_JAR
```

### "Python package missing"
```bash
# Reinstall environment
conda env remove -n vicast
conda env create -f environment.yml
```

---

## Minimal Installation (Pathway 2 only)

If you only need standard annotation (Pathway 2):

```bash
# 1. Conda environment
conda env create -f environment.yml
conda activate vicast

# 2. SnpEff setup
source setup/snpeff_env.sh

# That's it! No VADR or BLAST databases needed
```

**Storage: ~10 GB**
**Most users start here and add pathways as needed**

---

## Quick Reference

| Pathway | Conda Env | SnpEff | VADR | BLAST DB | Total Size |
|---------|-----------|--------|------|----------|------------|
| 1 - Already in SnpEff | ✓ | ✓ | - | - | ~10 GB |
| 2 - Well-annotated | ✓ | ✓ | - | - | ~10 GB |
| 3 - VADR | ✓ | ✓ | ✓ | - | ~17 GB |
| 4 - BLASTx | ✓ | ✓ | - | ✓ | 20-110 GB |
| All pathways | ✓ | ✓ | ✓ | ✓ | 27-127 GB |

---

For installation help, see: [VICAST README](../README.md)
