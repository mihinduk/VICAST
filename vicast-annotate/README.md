# VICAST-annotate

Semi-automated pipeline for curating viral genome annotations and integrating them into SnpEff databases.

## Overview

VICAST-annotate addresses the challenge of poorly annotated viral genomes in public databases by providing a systematic workflow to:
1. Download viral genomes from NCBI
2. Validate and improve annotations using VADR
3. Create custom SnpEff databases for variant effect prediction
4. Provide manual curation checkpoints for quality control

## Workflow

### Step 1: Parse and Validate Viral Genome
`step1_parse_viral_genome.py` handles genome retrieval and annotation validation.

**Input:** NCBI RefSeq accession (e.g., NC_001477)

**Process:**
- Downloads genome sequence and GenBank annotation from NCBI
- Converts GenBank to GFF3 format
- Runs VADR annotation validation
- Generates curated GFF3 file for manual review

**Output:**
- `{accession}.fasta` - Genome sequence
- `{accession}.gff3` - Initial annotation
- `{accession}_curated.gff3` - VADR-validated annotation
- `{accession}_vadr/` - VADR validation results

**Usage:**
```bash
python step1_parse_viral_genome.py NC_001477
```

**Manual Curation Checkpoint:**
After Step 1, review `{accession}_curated.gff3` and edit as needed. This is your opportunity to:
- Correct gene boundaries
- Add missing features
- Remove spurious annotations
- Validate gene names and products

### Step 2: Add to SnpEff Database
`step2_add_to_snpeff.py` integrates the curated genome into SnpEff.

**Input:** 
- Curated GFF3 file from Step 1
- Genome FASTA file
- NCBI accession

**Process:**
- Creates SnpEff database structure
- Generates genes.gbk and protein.fa files
- Updates SnpEff configuration
- Builds and validates the database
- Tests variant annotation

**Output:**
- Custom SnpEff database in `$SNPEFF_DATA/{accession}/`
- Updated `snpEff.config`
- Test variant annotations

**Usage:**
```bash
python step2_add_to_snpeff.py NC_001477
```

## Validation and Testing

`vicast_validation.py` provides comprehensive testing functions:

### Setup Validation
```bash
validate_vicast_setup
```
Checks:
- SnpEff installation and configuration
- VADR installation
- Python package dependencies
- Environment variables
- File permissions

### Pipeline Testing
```bash
test_vicast_annotate NC_001477
```
Runs complete annotation pipeline on a test virus with validation at each step.

### Individual Component Tests
```bash
validate_snpeff      # Test SnpEff installation
validate_vadr        # Test VADR installation
validate_python_packages  # Check Python dependencies
```

## Requirements

### Software
- Python 3.8+
- SnpEff 5.2+ (with Java 21)
- VADR 1.6+
- Biopython

### Environment Variables
Must be set (via `source setup/snpeff_env.sh`):
- `SNPEFF_JAR` - Path to SnpEff executable
- `SNPEFF_DATA` - Path to SnpEff data directory
- `VADR_DIR` - Path to VADR installation

## Input Files

### Genome Accession
- Must be a valid NCBI RefSeq accession (NC_XXXXXX format)
- Genome will be automatically downloaded from NCBI

### Manual Curation (Optional)
After Step 1, you can manually edit `{accession}_curated.gff3` to:
- Refine gene annotations
- Add biological knowledge not captured by automated tools
- Correct any VADR errors or warnings

## Output Files

### Step 1 Outputs
```
{accession}.fasta                # Genome sequence
{accession}.gff3                 # Initial GenBank annotation
{accession}_curated.gff3         # VADR-validated annotation (edit this!)
{accession}_vadr/                # VADR results directory
```

### Step 2 Outputs
```
$SNPEFF_DATA/{accession}/
├── genes.gbk                    # Gene annotations
├── protein.fa                   # Protein sequences
└── snpEffectPredictor.bin       # SnpEff database

Updated: $SNPEFF_CONFIG           # SnpEff configuration
```

## Usage Examples

### Basic Workflow
```bash
# 1. Setup environment
conda activate vicast
source setup/snpeff_env.sh

# 2. Run annotation pipeline
python vicast-annotate/step1_parse_viral_genome.py NC_001477

# 3. Review and edit NC_001477_curated.gff3 if needed

# 4. Create SnpEff database
python vicast-annotate/step2_add_to_snpeff.py NC_001477

# 5. Test the database
snpEff ann NC_001477 test.vcf > annotated.vcf
```

### Batch Processing
```bash
# Process multiple viruses
for accession in NC_001477 NC_001802 NC_003977; do
    python vicast-annotate/step1_parse_viral_genome.py $accession
    # Manual curation checkpoint
    read -p "Review ${accession}_curated.gff3 and press enter to continue..."
    python vicast-annotate/step2_add_to_snpeff.py $accession
done
```

### Testing Before Production Use
```bash
# Validate your setup
validate_vicast_setup

# Test with Dengue virus
test_vicast_annotate NC_001477

# If all tests pass, proceed with your viruses
```

## Troubleshooting

### Common Issues

**VADR fails with "model not found":**
- Ensure VADR models are installed: `bash setup/install_vadr.sh $SOFTWARE_DIR`
- Check `$VADR_DIR` is set correctly

**SnpEff build fails:**
- Verify Java 21 is installed: `java -version`
- Check SnpEff paths: `echo $SNPEFF_JAR $SNPEFF_DATA`
- Ensure write permissions to `$SNPEFF_DATA`

**NCBI download fails:**
- Check internet connectivity
- Verify accession is valid RefSeq ID (NC_XXXXXX)
- Try manual download and place files in current directory

**GFF3 parsing errors:**
- Review `{accession}_curated.gff3` for formatting issues
- Ensure all features have required attributes (ID, Name)
- Check for special characters in gene names

## Tips for Manual Curation

1. **Review VADR warnings:** Check `{accession}_vadr/` output for annotation issues
2. **Validate gene boundaries:** Ensure start/stop codons are correct
3. **Check overlapping features:** Resolve conflicts between genes
4. **Verify protein translations:** Use ORF prediction tools to confirm
5. **Reference literature:** Compare with published genome annotations

## Advanced Options

### Custom VADR Models
```bash
# Use specific VADR model
export VADR_MODEL_DIR=/path/to/custom/models
python step1_parse_viral_genome.py NC_XXXXXX --vadr-model custom_model
```

### Custom SnpEff Configuration
```bash
# Use alternate SnpEff config
export SNPEFF_CONFIG=/path/to/snpEff.config
python step2_add_to_snpeff.py NC_XXXXXX
```

## Version History

### v2.1.0
- Initial release with complete annotation pipeline
- VADR integration for validation
- Automated SnpEff database creation
- Comprehensive testing suite

## Contributing

Report issues or suggest improvements at: https://github.com/mihinduk/VICAST/issues

## Citation

```
Mihindukulasuriya KA, Handley SA. VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit. 
GitHub: https://github.com/mihinduk/VICAST (2024)
```

---

For more information, see the [main VICAST README](../README.md).
