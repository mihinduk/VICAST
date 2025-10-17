# VICAST

**V**iral **C**ultured-virus **A**nnotation and **S**npEff **T**oolkit

![VICAST Logo](vicast-annotate/VICAST_logo.png)

A comprehensive suite of semi-automated pipelines for cultured virus genomic analysis, specializing in annotation curation and variant calling for viral passage studies.

## Overview

VICAST provides systematic workflows for:
- Adding poorly annotated viral genomes to SnpEff databases with quality curation
- Comprehensive variant calling pipelines optimized for tissue culture passage analysis
- Integrated validation and testing functions
- Automated project initialization and environment setup

## Components

### VICAST-annotate
Pipeline for curating viral genome annotations and integrating them into SnpEff databases. Handles NCBI genome downloads, VADR annotation validation, and custom database creation.

**Key Features:**
- Automated NCBI genome retrieval
- VADR-based annotation validation
- Semi-automated curation with manual checkpoints
- SnpEff database integration
- Built-in testing and validation functions

### VICAST-analyze (Coming Soon)
Variant calling pipeline for cultured virus passage studies, optimized for identifying low-frequency variants and tracking evolutionary changes across passages.

## Installation

### Prerequisites
- Linux/Unix environment (tested on CentOS 7)
- Conda/Mamba package manager
- Git

### Quick Start

1. **Clone the repository:**
```bash
git clone https://github.com/mihinduk/VICAST.git
cd VICAST
```

2. **Create conda environment(s):**

**For most users (Pathways 1, 2, 4):**
```bash
conda env create -f environment.yml
conda activate vicast
```

**For Pathway 3 (VADR) - separate environment:**
```bash
conda env create -f environment_vadr.yml
conda activate vicast-vadr
```

**Why two environments?** VADR requires BLAST >=2.15.0, while other pathways work with BLAST 2.13.0. See [ENVIRONMENT_SETUP.md](docs/ENVIRONMENT_SETUP.md) for details.

3. **Setup SnpEff and dependencies:**
```bash
# Define your software installation directory
export SOFTWARE_DIR=/path/to/software

# Run setup script
source setup/setup_snpeff_custom_paths.sh $SOFTWARE_DIR
```

4. **Install VADR (required for annotation):**
```bash
bash setup/install_vadr.sh $SOFTWARE_DIR
```

5. **Validate installation:**
```bash
validate_vicast_setup
```

## Usage

### Initialize a new project
```bash
init_vicast_project my_virus_project
cd my_virus_project
```

### Run VICAST-annotate pipeline
```bash
# Test with a sample virus
test_vicast_annotate NC_001477

# Run on your virus
python vicast-annotate/step1_parse_viral_genome.py NC_XXXXXX
python vicast-annotate/step2_add_to_snpeff.py NC_XXXXXX
```

### Environment Management
```bash
# Load VICAST environment variables
source setup/snpeff_env.sh

# Verify paths
echo $SNPEFF_JAR
echo $SNPEFF_DATA
```

## Documentation

- [VICAST-annotate README](vicast-annotate/README.md) - Detailed annotation pipeline documentation
- [Setup Guide](docs/SETUP.md) - Installation and configuration guide
- [Examples](examples/) - Example workflows and test data

## Requirements

### Software Dependencies
- Python 3.8+
- Java 21 (for SnpEff 5.2+)
- SnpEff
- VADR (Viral Annotation DefineR)
- Biopython
- Standard bioinformatics tools (samtools, bcftools, etc.)

All Python dependencies are specified in `environment.yml`.

## Configuration

VICAST uses environment variables for flexible configuration:

- `SNPEFF_JAR` - Path to SnpEff executable
- `SNPEFF_DATA` - Path to SnpEff data directory
- `VADR_DIR` - Path to VADR installation
- `SCRATCH_DIR` - Working directory for temporary files

Set these via `setup/snpeff_env.sh` or export them directly.

## Testing

VICAST includes comprehensive validation functions:

```bash
# Validate complete setup
validate_vicast_setup

# Test annotation pipeline
test_vicast_annotate NC_001477

# Check specific components
validate_snpeff
validate_vadr
validate_python_packages
```

## Project Structure

```
VICAST/
├── vicast-annotate/          # Annotation pipeline scripts
│   ├── step1_parse_viral_genome.py
│   ├── step2_add_to_snpeff.py
│   └── vicast_validation.py
├── setup/                    # Installation and setup scripts
│   ├── setup_snpeff_custom_paths.sh
│   ├── install_vadr.sh
│   ├── install_java21_conda.sh
│   └── snpeff_env.sh
├── docs/                     # Documentation
├── examples/                 # Example data and workflows
├── environment.yml           # Conda environment specification
└── README.md                # This file
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

## Contact

For questions or support:
- GitHub Issues: https://github.com/mihinduk/VICAST/issues
- Maintainer: Kathie A. Mihindukulasuriya

## Acknowledgments

Developed in the Handley Lab for viral genomics research with a focus on cultured virus systems and passage studies.

## Version History

### v2.1.0 (2024)
- Initial public release
- Complete VICAST-annotate pipeline
- Automated setup and validation tools
- Comprehensive documentation

---

**VICAST** - Systematic viral genomics for cultured virus research
