# VICAST-annotate

Semi-automated pipeline for curating viral genome annotations and integrating them into SnpEff databases.

## Overview

VICAST-annotate addresses the challenge of poorly annotated viral genomes by providing **four systematic pathways** based on genome annotation quality and availability:

1. **Pathway 1**: Already in SnpEff → Use existing database
2. **Pathway 2**: Well-annotated (NCBI) → Standard pipeline with manual curation
3. **Pathway 3**: VADR annotation → Enhanced validation for poorly annotated genomes
4. **Pathway 4**: BLASTx annotation → Homology-based annotation for novel/unannotated genomes

## Quick Start

### Automatic Pathway Selection (Recommended)

```bash
# VICAST automatically determines the best pathway
python3 vicast_annotate.py NC_001477
```

The master controller will:
1. Check if genome is already in SnpEff
2. Assess NCBI annotation quality
3. Select and execute the appropriate pathway
4. Provide manual curation checkpoint
5. Add genome to SnpEff database

### Manual Pathway Selection

```bash
# Force specific pathway
python3 vicast_annotate.py NC_001477 --pathway 2

# Auto-proceed without manual curation (not recommended)
python3 vicast_annotate.py NC_001477 --auto
```

## The Four Pathways

### Pathway 1: Already in SnpEff

**When to use**: Genome is already in SnpEff database (built-in or custom)

**Command**:
```bash
python3 step0_check_snpeff.py NC_001477
```

**Result**: No action needed, genome ready to use

**Usage**:
```bash
snpeff NC_001477 variants.vcf > annotated.vcf
```

---

### Pathway 2: Well-Annotated (NCBI)

**When to use**: NCBI has good quality annotations (most CDS have products and gene names)

**Workflow**:
```bash
# Step 0: Verify pathway
python3 step0_check_snpeff.py NC_001477

# Step 1: Parse GenBank, skip polyproteins
python3 step1_parse_viral_genome.py NC_001477

# Manual curation: Review and edit TSV file
# Edit: NC_001477_no_polyprotein.tsv

# Step 2: Add to SnpEff
python3 step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv
```

**Key Features**:
- Automatically skips polyprotein features
- Creates editable TSV for manual curation
- Preserves gene names and products from NCBI
- Handles frameshifts

---

### Pathway 3: VADR-Enhanced Annotation

**When to use**: NCBI annotations are incomplete but genome has recognizable viral features

**Workflow**:
```bash
# Step 1: Parse with VADR validation
python3 step1_parse_viral_genome.py NC_009942.1 --use-vadr

# VADR will:
# - Validate existing annotations
# - Identify annotation errors
# - Generate quality reports

# Manual curation: Review VADR alerts and edit TSV
# Files: NC_009942.1_vadr_curated.tsv
#        NC_009942.1_vadr/ (VADR results)

# Step 2: Add to SnpEff
python3 step2_add_to_snpeff.py NC_009942.1 NC_009942.1_vadr_curated.tsv
```

**VADR Models**:
- `flavi` - Flaviviruses (default)
- `corona` - Coronaviruses
- `calici` - Caliciviruses
- `flua` - Influenza A
- Custom models in `$VADR_DIR/vadr-models/`

**Custom model**:
```bash
python3 step1_parse_viral_genome.py NC_XXXXXX --use-vadr --vadr-model corona
```

---

### Pathway 4: BLASTx-Based Annotation

**When to use**: No quality annotations available, need homology-based annotation

**Workflow**:
```bash
# Ensure you have a FASTA file
# If from NCBI, download first

# Step 1: BLASTx annotation
python3 step1_blastx_annotate.py NC_XXXXXX.fasta --blast-db nr

# BLASTx will:
# - Search against protein database
# - Identify homologous proteins
# - Merge overlapping hits
# - Assign gene names and products

# Manual curation: Review and refine BLASTx results
# Edit: NC_XXXXXX_blastx.tsv
# - Verify gene boundaries
# - Check for missed ORFs
# - Refine annotations

# Step 2: Add to SnpEff
python3 step2_add_to_snpeff.py NC_XXXXXX NC_XXXXXX_blastx.tsv
```

**BLAST Databases**:
- `nr` - NCBI non-redundant protein database (comprehensive but slow)
- `refseq_protein` - RefSeq proteins (faster)
- `/path/to/viral_db` - Custom viral protein database (recommended)

**Tips**:
- Use viral-specific database for faster, more relevant results
- Adjust E-value for sensitivity: `--evalue 1e-10` (stricter) or `--evalue 1e-3` (relaxed)
- BLASTx results need careful manual curation

## Requirements

### Software
- Python 3.8+
- SnpEff 5.2+ (with Java 21)
- Biopython

### Pathway-Specific Requirements

**Pathway 3 (VADR)**:
- VADR 1.6+ installed
- VADR models for your virus family
- Install: `bash setup/install_vadr.sh /path/to/software`

**Pathway 4 (BLASTx)**:
- BLAST+ installed: `conda install -c bioconda blast`
- Access to BLAST database (nr, refseq_protein, or custom)

### Environment Variables

Must be set before running:
```bash
# Load VICAST environment
source setup/snpeff_env.sh

# Or set manually:
export SNPEFF_JAR=/path/to/snpEff.jar
export SNPEFF_DATA=/path/to/snpEff/data
export VADR_DIR=/path/to/vadr  # For pathway 3
```

## Decision Tree

```
Start: Check genome ID
    │
    ├──> Already in SnpEff? ──YES──> Pathway 1: Done!
    │         │
    │        NO
    │         │
    │         v
    ├──> NCBI has good annotations? ──YES──> Pathway 2: Standard
    │         │
    │        NO
    │         │
    │         v
    ├──> VADR model available? ──YES──> Pathway 3: VADR
    │         │
    │        NO
    │         │
    │         v
    └──> Use BLASTx ──────────────────────> Pathway 4: BLASTx
```

## Manual Curation Guidelines

All pathways include a manual curation checkpoint. Use this to:

1. **Verify gene boundaries**
   - Check start/stop codons
   - Ensure no frame shifts (unless expected)

2. **Check gene names**
   - Use standard nomenclature
   - Be consistent across similar viruses

3. **Refine products**
   - Use descriptive names
   - Avoid vague terms like "hypothetical protein"

4. **Handle special features**
   - Polyproteins (usually removed in pathways 2-3)
   - Programmed frameshifts
   - Overlapping genes
   - Ribosomal slippage sites

5. **Remove artifacts**
   - False positive ORFs (pathway 4)
   - Annotation errors
   - Redundant features

## Validation and Testing

VICAST includes comprehensive validation:

```bash
# Validate complete setup
validate_vicast_setup

# Test pathway
test_vicast_annotate NC_001477

# Check specific components
validate_snpeff
validate_vadr
validate_python_packages
```

## Output Files

### All Pathways
- `{genome_id}.fasta` - Genome sequence
- `{genome_id}_*.tsv` - Editable annotation table
- `{genome_id}_*.gff3` - GFF3 annotation file

### Pathway-Specific
**Pathway 2**:
- `{genome_id}_no_polyprotein.tsv`
- `{genome_id}_no_polyprotein.gff3`

**Pathway 3**:
- `{genome_id}_vadr/` - VADR validation results
- `{genome_id}_vadr_curated.tsv`
- `{genome_id}_vadr.gff3`

**Pathway 4**:
- `{genome_id}_blastx.txt` - BLAST results
- `{genome_id}_blastx.tsv`

### SnpEff Database
```
$SNPEFF_DATA/{genome_id}/
├── sequences.fa
├── genes.gff
└── snpEffectPredictor.bin
```

## Troubleshooting

### Common Issues

**"Genome not found in SnpEff"** (Pathway 1)
- Run pathway 2, 3, or 4 to add it

**"No CDS features found"** (Pathway 2)
- Try pathway 3 (VADR) or pathway 4 (BLASTx)

**"VADR model not found"** (Pathway 3)
- Install VADR models: `bash setup/install_vadr.sh`
- Check `$VADR_DIR/vadr-models/`
- Use different model: `--vadr-model corona`

**"BLASTx no hits"** (Pathway 4)
- Use viral-specific database
- Relax E-value: `--evalue 1e-3`
- Check if genome is truly novel

**"SnpEff build failed"**
- Check Java version: `java -version` (need Java 21)
- Verify paths: `echo $SNPEFF_JAR $SNPEFF_DATA`
- Check GFF3 format validity

## Examples

### Example 1: Dengue Virus (Well-Annotated)
```bash
python3 vicast_annotate.py NC_001477
# → Pathway 2: Standard pipeline
# → Manual curation of polyprotein cleavage sites
# → Added to SnpEff in ~5 minutes
```

### Example 2: Novel Flavivirus (Poor Annotations)
```bash
python3 vicast_annotate.py NC_XXXXXX --pathway 3
# → Pathway 3: VADR with flavi model
# → VADR identifies and corrects annotation errors
# → Manual refinement of gene boundaries
```

### Example 3: Uncharacterized Virus (No Annotations)
```bash
python3 vicast_annotate.py novel_virus.fasta --pathway 4 --blast-db /data/viral_proteins
# → Pathway 4: BLASTx against viral database
# → Identifies homologous proteins
# → Extensive manual curation needed
```

### Example 4: Already Available
```bash
python3 step0_check_snpeff.py NC_045512
# → Pathway 1: SARS-CoV-2 already in SnpEff
# → Ready to use immediately
```

## Version History

### v2.1.0
- Initial release with 4-pathway system
- Automatic pathway detection
- VADR integration
- BLASTx homology-based annotation
- Comprehensive testing suite

## Citation

```
Mihindukulasuriya KA, Handley SA. VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit.
GitHub: https://github.com/mihinduk/VICAST (2024)
```

## Contributing

Report issues or suggest improvements at: https://github.com/mihinduk/VICAST/issues

---

For main VICAST documentation, see: [../README.md](../README.md)

---

## Special Case: Segmented Viruses

### Overview

Some viruses have segmented genomes (multiple chromosomes), such as:
- **Influenza A/B** (8 segments: PB2, PB1, PA, HA, NP, NA, M, NS)
- **Rotavirus** (11 segments)
- **Bunyaviruses** (3 segments)
- **Reoviruses** (10-12 segments)

VICAST provides a specialized script to handle segmented viruses by combining all segments into a single SnpEff database.

### Standard Approach (Separate Databases)

Run VICAST-annotate on each segment individually:
```bash
# Influenza example - creates 8 separate databases
python3 vicast_annotate.py CY121680  # PB2
python3 vicast_annotate.py CY121681  # PB1
python3 vicast_annotate.py CY121682  # PA
# ... etc for all 8 segments
```

**Pros**: Simple, follows NCBI structure
**Cons**: Need to run 8 times, manage 8 databases

### Recommended Approach (Combined Database)

Use `vicast_annotate_segmented.py` to combine all segments into one database:

```bash
# Influenza A (H1N1) 2009 pandemic strain
python3 vicast_annotate_segmented.py influenza_h1n1_2009 \
  --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
  --names PB2,PB1,PA,HA,NP,NA,M,NS
```

**What it does:**
1. Downloads all 8 segments from NCBI
2. Combines into single `sequences.fa` (multi-FASTA)
3. Combines annotations into single `genes.gff`
4. Builds unified SnpEff database `influenza_h1n1_2009`

**Usage:**
```bash
# Annotate variants (VCF CHROM must match segment names)
snpeff influenza_h1n1_2009 variants.vcf > annotated.vcf
```

### VCF Chromosome Naming

**Important**: Your VCF file's CHROM field must match the segment names:

**Using segment names (recommended):**
```vcf
#CHROM    POS    ID    REF    ALT    QUAL    FILTER    INFO
PB2       100    .     A      G      .       .         .
HA        500    .     C      T      .       .         .
```

**Using accessions (alternative):**
```vcf
#CHROM       POS    ID    REF    ALT    QUAL    FILTER    INFO
CY121680     100    .     A      G      .       .         .
CY121683     500    .     C      T      .       .         .
```

Set this when combining segments with `--names`:
```bash
# Use segment names (PB2, PB1, PA...)
python3 vicast_annotate_segmented.py influenza_h1n1_2009 \
  --segments CY121680,CY121681,... \
  --names PB2,PB1,PA,...

# Or use accessions (CY121680, CY121681...)
python3 vicast_annotate_segmented.py influenza_h1n1_2009 \
  --segments CY121680,CY121681,...
  # No --names flag = uses accessions
```

### Using Local Files

If you already have segment files downloaded:

```bash
python3 vicast_annotate_segmented.py my_influenza \
  --fasta-files PB2.fasta,PB1.fasta,PA.fasta,HA.fasta,NP.fasta,NA.fasta,M.fasta,NS.fasta \
  --gb-files PB2.gb,PB1.gb,PA.gb,HA.gb,NP.gb,NA.gb,M.gb,NS.gb \
  --names PB2,PB1,PA,HA,NP,NA,M,NS
```

### Common Segmented Viruses

#### Influenza A
```bash
python3 vicast_annotate_segmented.py influenza_a_strain_X \
  --segments ACC1,ACC2,ACC3,ACC4,ACC5,ACC6,ACC7,ACC8 \
  --names PB2,PB1,PA,HA,NP,NA,M,NS
```

#### Rotavirus
```bash
python3 vicast_annotate_segmented.py rotavirus_strain_X \
  --segments ACC1,ACC2,ACC3,ACC4,ACC5,ACC6,ACC7,ACC8,ACC9,ACC10,ACC11 \
  --names VP1,VP2,VP3,VP4,NSP1,VP6,NSP3,NSP2,VP7,NSP4,NSP5
```

#### Bunyavirus (3 segments)
```bash
python3 vicast_annotate_segmented.py bunyavirus_strain_X \
  --segments ACC1,ACC2,ACC3 \
  --names L,M,S
```

### Troubleshooting

**"Chromosome not found in database"**
- Check VCF CHROM field matches segment names used in `--names`
- Use `snpeff dump genome_id | head` to see chromosome names

**"No sequences combined"**
- Verify all accessions are valid
- Check NCBI downloads completed successfully
- Ensure files exist if using `--fasta-files`

**Polyprotein handling**
- By default, polyproteins are skipped (use `--skip-polyprotein`)
- To keep polyproteins: use `--keep-polyprotein`
- Individual mature peptides are preserved

### Benefits of Combined Database

✓ Single database for entire segmented genome
✓ One command to annotate all segments
✓ Easier to manage and share
✓ Consistent with SnpEff's multi-chromosome support
✓ Works with standard variant calling pipelines

---

For questions about segmented virus annotation, see: [GitHub Issues](https://github.com/mihinduk/VICAST/issues)
