# VICAST

**V**iral **C**ultured-virus **A**nnotation and **S**npEff **T**oolkit

![VICAST Logo](vicast-annotate/VICAST_logo.png)

[![CI](https://github.com/mihinduk/VICAST/actions/workflows/ci.yml/badge.svg)](https://github.com/mihinduk/VICAST/actions/workflows/ci.yml)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A comprehensive suite of semi-automated pipelines for cultured virus genomic analysis, specializing in annotation curation and variant calling for viral passage studies.

## Overview

VICAST provides systematic workflows for:
- Adding poorly annotated viral genomes to SnpEff databases with quality curation
- Comprehensive variant calling pipelines optimized for tissue culture passage analysis
- Integrated validation and testing functions
- Automated project initialization and environment setup

## Components

### VICAST-annotate
Pipeline for curating viral genome annotations and integrating them into SnpEff databases. Handles NCBI genome downloads, BLASTx annotation, and custom database creation.

**Key Features:**
- Automated NCBI genome retrieval
- BLASTx-based homology annotation
- Semi-automated curation with manual checkpoints
- SnpEff database integration
- Segmented virus support
- Built-in testing and validation functions

### VICAST-analyze
Variant calling pipeline for cultured virus passage studies, optimized for identifying low-frequency variants and tracking evolutionary changes across passages.

**Key Features:**
- Quality control with fastp
- Read mapping with bwa
- Low-frequency variant calling with lofreq
- SnpEff variant annotation
- Contamination detection with de novo assembly
- Multi-tier variant analysis (high/medium/low frequency)
- BAM read-level co-occurrence analysis for variant validation
- Consensus genome generation from dominant variants

## Installation

### Option 1: Docker (Recommended for Portability)

Docker provides a self-contained environment with all dependencies pre-installed.

#### Building the Image

```bash
# Clone repository
git clone https://github.com/mihinduk/VICAST.git
cd VICAST

# Build Docker image
docker build -t vicast:latest .
```

#### Running VICAST in Docker

**Prerequisites:**

VICAST requires two persistent directories on your host system:
1. **Data directory**: For input files (FASTQ) and output results (VCF, HTML reports)
2. **Database directory**: For SnpEff viral genome databases (reusable across analyses)

**Setup - Create Directories:**

```bash
# Option 1: Home directory (recommended for laptops/workstations)
mkdir -p ~/vicast_data ~/vicast_databases

# Option 2: Scratch space (recommended for HPC/servers with limited home quota)
mkdir -p /scratch/$USER/vicast_data /scratch/$USER/vicast_databases

# Option 3: Custom location (adjust paths as needed)
mkdir -p /path/to/your/vicast_data /path/to/your/vicast_databases
```

**Quick Start - Pre-built Database Workflow:**

VICAST provides 20+ pre-built viral genome databases available from GitHub.

```bash
# Start VICAST container with directory mounts
# Replace ~/vicast_data and ~/vicast_databases with your actual paths if using Option 2 or 3 above
docker run -it --rm \
  -v ~/vicast_data:/data \
  -v ~/vicast_databases:/opt/vicast/snpeff_data_custom \
  -w /data \
  --user $(id -u):$(id -g) \
  --hostname vicast \
  vicast:latest

# You are now inside the container
# Prompt shows: (base) I have no name!@vicast:/data$
# The "I have no name!" is expected when using --user flag - it's harmless

# Step 1: List available pre-built databases
install_prebuilt_database.sh --list
# Shows 20+ viruses: DENV-1/2, SARS-CoV-2, Influenza A, WNV, ZIKV, etc.

# Step 2: Install a database (example: Dengue virus 2)
install_prebuilt_database.sh --install NC_001474.2

# Step 3: Download test data from SRA (optional)
download_sra_data.sh SRR5992153

# Step 4: Run complete analysis pipeline
run_vicast_analyze_full.sh SRR5992153_1.fastq.gz SRR5992153_2.fastq.gz NC_001474.2 8

# Results appear in ~/vicast_data on your host:
# - QC reports (HTML)
# - Variant calls (VCF)
# - Annotated variants (TSV)
# - Coverage plots

# Exit container when done
exit
```

**Understanding the Docker Command:**

```bash
docker run -it --rm \
  -v ~/vicast_data:/data \
  -v ~/vicast_databases:/opt/vicast/snpeff_data_custom \
  -w /data \
  --user $(id -u):$(id -g) \
  --hostname vicast \
  vicast:latest
```

**What each parameter does:**
- `-it`: Interactive terminal (allows you to type commands)
- `--rm`: Remove container after exit (saves disk space)
- `-v ~/vicast_data:/data`: Mount your data directory into the container at `/data`
- `-v ~/vicast_databases:/opt/vicast/snpeff_data_custom`: Mount database directory (databases persist!)
- `-w /data`: Set working directory to `/data` inside container
- `--user $(id -u):$(id -g)`: Run as your user ID (ensures correct file permissions)
- `--hostname vicast`: Set container hostname to "vicast" (shows in prompt)
- `vicast:latest`: The Docker image to run

**Note:** The prompt will show `(base) I have no name!@vicast:/data$` - the "I have no name!" is expected when using `--user` flag and is harmless. It indicates you're running with your host user ID for correct file permissions.

**Where Your Files Are:**

After running the analysis, find your results on the host system:

```bash
# On your host machine (not in container):
ls ~/vicast_data/
# Shows: FASTQ files, VCF, HTML reports, annotated TSV

ls ~/vicast_databases/
# Shows: NC_001474.2/ (database directory)
# This database is now available for all future runs!
```

**Reusing Installed Databases:**

Once a database is installed, it persists in `~/vicast_databases/`. Next time you run the container with the same database mount, the database is already available - no need to reinstall!

```bash
# Restart container later with same mounts
docker run -it --rm \
  -v ~/vicast_data:/data \
  -v ~/vicast_databases:/opt/vicast/snpeff_data_custom \
  -w /data \
  --user $(id -u):$(id -g) \
  --hostname vicast \
  vicast:latest

# Database is already there!
ls /opt/vicast/snpeff_data_custom/NC_001474.2/

# Run analysis immediately
run_vicast_analyze_full.sh new_sample_R1.fastq.gz new_sample_R2.fastq.gz NC_001474.2 8
```

**Custom Genome Annotation:**

For genomes not in the pre-built collection, use VICAST-annotate pathways:

```bash
# Same container setup - REPLACE paths if you used Option 2 or 3 above
docker run -it --rm \
  -v ~/vicast_data:/data \
  -v ~/vicast_databases:/opt/vicast/snpeff_data_custom \
  -w /data \
  --user $(id -u):$(id -g) \
  --hostname vicast \
  vicast:latest

# Inside container - download and annotate a new genome from NCBI
step1_parse_viral_genome.py AF252854.1  # Example: Bovine adenovirus 2

# This creates: AF252854.1.tsv
# Edit the TSV file if you need to curate gene names/features
# (Use exit to leave container, edit on host, then restart container)

# Add the curated genome to SnpEff
step2_add_to_snpeff.py AF252854.1 AF252854.1.tsv

# Database now saved to ~/vicast_databases/AF252854.1/
# Run analysis with your new genome
run_vicast_analyze_full.sh R1.fastq.gz R2.fastq.gz AF252854.1 8
```

**Available Pre-built Databases (v2.3.0+):**
- DENV-1, DENV-2 (Dengue virus)
- SARS-CoV-2 (COVID-19)
- Influenza A (PR8, 8 segments)
- West Nile virus (multiple strains)
- Zika, Sindbis, Ross River, VEEV
- And 10+ more viruses

Run `install_prebuilt_database.sh --list` to see all available databases.

**Important Usage Notes:**
- ✅ **User mapping (`--user`)**: Required for correct file permissions on host
- ✅ **Database mount**: `-v ~/vicast_databases:/opt/vicast/snpeff_data_custom` persists databases across runs
- ✅ **Data mount**: `-v ~/vicast_data:/data` for input/output files
- ✅ **Space requirements**: ~20-50GB for analysis intermediate files
- ✅ **Database persistence**: Installed databases available in all future runs with same mount

### Option 2: Conda Installation

#### Prerequisites
- Linux/Unix environment (Linux, macOS, WSL2)
- Conda or Mamba package manager
- Git

#### Quick Start

1. **Clone the repository:**
```bash
git clone https://github.com/mihinduk/VICAST.git
cd VICAST
```

2. **Create conda environment:**

For annotation only (lighter, ~3-4 GB):
```bash
conda env create -f environment_vicast_annotate.yml
conda activate vicast_annotate
```

For full analysis pipeline (~5-6 GB):
```bash
conda env create -f environment_vicast_analyze.yml
conda activate vicast_analyze
```

3. **Install VICAST Python package:**
```bash
pip install -e .
```

4. **Install SnpEff (if not already installed):**
```bash
# Download SnpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip

# Set environment variables
export SNPEFF_HOME="$(pwd)/snpEff"
export SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
export SNPEFF_DATA="${SNPEFF_HOME}/data"
```

5. **Configure environment variables:**
```bash
# Copy and edit the configuration template
cp vicast_config.template.sh ~/.vicast/config.sh

# Edit with your paths
nano ~/.vicast/config.sh

# Add to your ~/.bashrc or ~/.zshrc:
source ~/.vicast/config.sh
```

6. **Validate installation:**
```bash
python -c "from vicast.config import get_config; get_config().print_status()"
```

### Option 3: Singularity (for HPC)

```bash
# Build Singularity image from Docker
singularity build vicast.sif docker://ghcr.io/mihinduk/vicast:latest

# Run with bind mounts
singularity exec --bind /data:/data vicast.sif \
    python /opt/vicast/vicast-annotate/step1_parse_viral_genome.py NC_001477
```

## Configuration

VICAST uses environment variables for flexible configuration across different systems:

| Variable | Description | Example |
|----------|-------------|---------|
| `SNPEFF_HOME` | SnpEff installation directory | `/opt/snpEff` |
| `SNPEFF_JAR` | Path to snpEff.jar | `${SNPEFF_HOME}/snpEff.jar` |
| `SNPEFF_DATA` | SnpEff data directory | `${SNPEFF_HOME}/data` |
| `JAVA_HOME` | Java installation (Java 21+ required) | `/usr/lib/jvm/java-21` |
| `NCBI_EMAIL` | Email for NCBI Entrez queries | `user@institution.edu` |
| `VICAST_SCRATCH` | Scratch directory for temp files | `/scratch/$USER/vicast` |
| `VICAST_THREADS` | Default thread count | `4` |

See `vicast_config.template.sh` for a complete configuration template.

## Usage

### VICAST-annotate Pipeline

```bash
# Pathway detection - check if genome is in SnpEff
python vicast-annotate/step0_check_snpeff.py NC_001477

# Pathway 2: Parse well-annotated NCBI genome
python vicast-annotate/step1_parse_viral_genome.py NC_001477
# Edit the TSV file as needed, then:
python vicast-annotate/step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv

# Pathway 3: BLASTx for poorly annotated genomes
python vicast-annotate/step1_blastx_annotate.py genome.fasta --blast-db nr
# Edit the TSV file, then add to SnpEff

# Pathway 4: Segmented viruses (e.g., influenza)
python vicast-annotate/vicast_annotate_segmented.py influenza \
    --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
    --names PB2,PB1,PA,HA,NP,NA,M,NS
```

### VICAST-analyze Pipeline

```bash
# QC workflow (Steps 1-6) - run first
./vicast-analyze/run_vicast_analyze_qc_only.sh R1.fastq.gz R2.fastq.gz NC_001477 4

# Review QC outputs, then continue with annotation (Steps 7-9)
./vicast-analyze/run_vicast_analyze_annotate_only.sh R1.fastq.gz R2.fastq.gz NC_001477

# Or run full pipeline
./vicast-analyze/run_vicast_analyze_full.sh R1.fastq.gz R2.fastq.gz NC_001477 4
```

## Documentation

- [VICAST-annotate README](vicast-annotate/README.md) - Detailed annotation pipeline documentation
- [VICAST-analyze Workflow](vicast-analyze/WORKFLOW_ARCHITECTURE.md) - Analysis pipeline architecture
- [Configuration Guide](vicast_config.template.sh) - Environment configuration template
- [Environment Setup](ENVIRONMENT_README.md) - Detailed environment setup guide

## Requirements

### Software Dependencies
- Python 3.9+
- Java 21+ (for SnpEff 5.2+)
- SnpEff 5.2+
- Biopython
- BLAST+ (for BLASTx annotation pathway)

For analysis pipeline:
- bwa
- samtools
- bcftools
- lofreq
- fastp

All dependencies are specified in the conda environment files.

## Testing

```bash
# Run Python tests
pytest tests/ -v

# Validate configuration
python -c "from vicast.config import get_config; get_config().print_status()"
```

## Project Structure

```
VICAST/
├── src/vicast/               # Python package
│   ├── config.py             # Centralized configuration
│   ├── annotate/             # Annotation module
│   └── analyze/              # Analysis module
├── vicast-annotate/          # Annotation pipeline scripts
│   ├── step0_check_snpeff.py
│   ├── step1_parse_viral_genome.py
│   ├── step1_blastx_annotate.py
│   ├── step2_add_to_snpeff.py
│   └── vicast_validation.py
├── vicast-analyze/           # Analysis pipeline scripts
│   ├── run_vicast_analyze_qc_only.sh
│   ├── run_vicast_analyze_annotate_only.sh
│   └── viral_pipeline.py
├── tests/                    # Test suite
├── docs/                     # Documentation
├── environment*.yml          # Conda environment files
├── pyproject.toml            # Python package configuration
├── Dockerfile                # Docker build file
└── vicast_config.template.sh # Configuration template
```

## Citation

If you use VICAST in your research, please cite:

```
Mihindukulasuriya KA, Handley SA. VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit.
GitHub: https://github.com/mihinduk/VICAST (2024)
```

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

### Development Setup

```bash
git clone https://github.com/mihinduk/VICAST.git
cd VICAST
pip install -e ".[dev]"
pytest tests/
```

## Contact

For questions or support:
- GitHub Issues: https://github.com/mihinduk/VICAST/issues
- Maintainer: Kathie A. Mihindukulasuriya

## Acknowledgments

Developed in the Handley Lab for viral genomics research with a focus on cultured virus systems and passage studies.

## Version History

### v2.2.0 (2025) - Unified Release
- **Architectural milestone:** Unified VICAST-annotate and VICAST-analyze into single integrated toolkit
- Portable configuration system across Docker/Conda/HPC
- Docker and Singularity support for reproducibility
- Python package structure with pyproject.toml
- CI/CD with GitHub Actions
- Comprehensive test suite (5 test modules)
- 8 integrated user guides (>3,200 lines of documentation)

**Development history:** VICAST originated as two independent repositories ([VICAST-analyze](https://github.com/mihinduk/VICAST-analyze): July 2025, 89 commits; [VICAST-annotate](https://github.com/mihinduk/VICAST-annotate): September 2025, 19 commits) that were unified in October 2025 into this integrated system. See [DEVELOPMENT_HISTORY.md](DEVELOPMENT_HISTORY.md) for complete repository evolution documentation.

### v2.1.0 (2024) - Initial Standalone Releases
- VICAST-annotate: Complete annotation curation pipeline
- VICAST-analyze: Variant calling pipeline
- Automated setup and validation tools
- Initial documentation

---

**VICAST** - Systematic viral genomics for cultured virus research

---

**Note:** Conda installation files are preserved in the `conda-legacy` branch for reference.
