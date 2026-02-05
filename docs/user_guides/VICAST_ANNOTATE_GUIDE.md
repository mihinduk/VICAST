# VICAST-Annotate: Genome Annotation Guide

**Module:** vicast-annotate
**Version:** 2.2.0
**Last Updated:** 2026-02-05

---

## Table of Contents

1. [Overview](#overview)
2. [The Four Annotation Pathways](#the-four-annotation-pathways)
3. [Step-by-Step: Pathway 2 (Standard)](#step-by-step-pathway-2-standard)
4. [Pathway 3: BLASTx Homology Annotation](#pathway-3-blastx-homology-annotation)
5. [Pathway 4: Segmented Viruses](#pathway-4-segmented-viruses)
6. [Manual Curation Guidelines](#manual-curation-guidelines)
7. [Troubleshooting](#troubleshooting)
8. [Examples for Different Genome Types](#examples-for-different-genome-types)

---

## Overview

VICAST-annotate prepares viral genomes for variant annotation by creating curated SnpEff databases. It addresses the common challenge of **poorly annotated viral genomes** that lack quality gene annotations.

### Why Use VICAST-annotate?

**Problem:** Many viral genomes in NCBI have:
- Missing or incomplete gene annotations
- Polyprotein features that obscure individual proteins
- Generic names like "hypothetical protein"
- No SnpEff database available

**Solution:** VICAST-annotate provides systematic workflows to:
- Download genomes from NCBI
- Parse existing annotations or create new ones via BLASTx
- Curate annotations with manual checkpoints
- Build SnpEff databases with viral-friendly settings
- Support segmented genomes (influenza, rotavirus, etc.)

### Key Features

- ✅ **Four pathways** for different annotation scenarios
- ✅ **Automatic pathway detection** based on genome quality
- ✅ **Manual curation checkpoints** for quality control
- ✅ **Polyprotein handling** - skip or process cleavage products
- ✅ **Segmented virus support** - combine multiple chromosomes
- ✅ **BLASTx annotation** for poorly characterized genomes

---

## The Four Annotation Pathways

VICAST-annotate provides **four systematic pathways** based on genome availability and annotation quality:

| Pathway | Scenario | Time | Manual Work |
|---------|----------|------|-------------|
| **1** | Already in SnpEff | 1 min | None - ready to use |
| **2** | Well-annotated (NCBI) | 5-10 min | Light curation |
| **3** | BLASTx annotation | 30-60 min | Heavy curation |
| **4** | Segmented virus | 10-20 min | Moderate curation |

### Decision Tree

```
Start: Do you have a genome accession?
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
    ├──> Multiple segments? ──YES──> Pathway 4: Segmented
    │         │
    │        NO
    │         │
    │         v
    └──> Poorly annotated ────────────────> Pathway 3: BLASTx
```

### Pathway Selection

**Use automatic detection (recommended):**

```bash
python vicast_annotate.py NC_001477
```

**Or manually select:**

```bash
python vicast_annotate.py NC_001477 --pathway 2
```

---

## Step-by-Step: Pathway 2 (Standard)

**Best for:** Well-annotated NCBI genomes (most common viruses)

**Examples:** Dengue virus, Zika virus, HIV, Poliovirus, many others

### Prerequisites

```bash
# Activate environment
conda activate vicast_annotate

# Verify SnpEff variables
echo $SNPEFF_JAR
echo $SNPEFF_DATA

# Set NCBI email (required for downloads)
export NCBI_EMAIL="your_email@institution.edu"
```

### Step 0: Check if Already in SnpEff

```bash
cd /path/to/VICAST/vicast-annotate

python step0_check_snpeff.py NC_001477
```

**Output:**
```
Checking if NC_001477 is in SnpEff database...
❌ NC_001477 NOT found in SnpEff
Recommended: Use Pathway 2 (standard annotation)
```

### Step 1: Parse Genome and Annotations

```bash
python step1_parse_viral_genome.py NC_001477
```

**What this does:**
1. Downloads FASTA from NCBI
2. Downloads GenBank annotations
3. Parses CDS features
4. **Automatically skips polyprotein features**
5. Creates editable TSV file

**Output files:**
```
NC_001477.fasta                    # Genome sequence
NC_001477_no_polyprotein.tsv       # Editable annotation table
```

**Example TSV content:**

```
seqname  source  feature  start  end    score  strand  frame  gene_name      product
NC_001477  NCBI  CDS      95     1491   .      +       0      Protein_C      capsid protein C
NC_001477  NCBI  CDS      95     794    .      +       0      Protein_prM-M  precursor membrane protein prM
NC_001477  NCBI  CDS      795    2393   .      +       0      Protein_E      envelope protein E
...
```

### Step 2: Manual Curation (Critical!)

**Open the TSV file:**

```bash
nano NC_001477_no_polyprotein.tsv
# Or use Excel, LibreOffice, etc.
```

**What to check:**

1. **Gene names** - Use standard nomenclature
   - ✅ Good: `Protein_E`, `NS5`, `capsid`
   - ❌ Bad: `hypothetical_protein`, `unnamed`, `gene1`

2. **Product descriptions** - Be specific
   - ✅ Good: `envelope protein E`, `RNA-dependent RNA polymerase`
   - ❌ Bad: `unknown function`, `protein`

3. **Start/end coordinates** - Verify in-frame
   - Check that (end - start + 1) is divisible by 3
   - Verify no frameshifts (unless biological)

4. **Overlapping features** - Decide which to keep
   - Some viruses have legitimate overlapping genes
   - Remove redundant annotations

**Example curation:**

```
Before:
gene_name: hypothetical_protein
product: unknown

After (search literature):
gene_name: NS4B
product: non-structural protein 4B
```

### Step 3: Add to SnpEff Database

```bash
python step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv
```

**What this does:**
1. Converts TSV to GFF3 format
2. Adds genome to SnpEff configuration
3. Builds SnpEff database
4. Validates database structure

**Output:**
```
✓ Created NC_001477.gff3
✓ Added NC_001477 to snpEff.config
✓ Building SnpEff database...
✓ Database built successfully!
✓ Validation: 10 genes, 10 transcripts, 10 exons

SnpEff database ready: NC_001477
```

**Database location:**
```
$SNPEFF_DATA/NC_001477/
├── sequences.fa
├── genes.gff
└── snpEffectPredictor.bin
```

### Step 4: Verify Database

```bash
# Check database exists
java -jar $SNPEFF_JAR databases | grep NC_001477

# Test annotation (if you have a VCF)
java -jar $SNPEFF_JAR NC_001477 test.vcf > annotated.vcf
```

---

## Pathway 3: BLASTx Homology Annotation

**Best for:** Poorly annotated or novel viral genomes

**When to use:**
- NCBI annotations are missing or low-quality
- Novel virus with limited characterization
- Need to annotate based on homology to known viruses

### Prerequisites

```bash
# Requires BLAST+ tools
conda install -c bioconda blast

# Need BLAST database
# Option 1: NCBI nr (comprehensive but slow)
# Option 2: Custom viral protein database (recommended)
```

### Step 1: Prepare FASTA File

```bash
# If from NCBI, download:
python -c "
from Bio import Entrez
Entrez.email = 'your@email.com'
handle = Entrez.efetch(db='nucleotide', id='NC_XXXXXX', rettype='fasta')
with open('NC_XXXXXX.fasta', 'w') as f:
    f.write(handle.read())
"

# Or use existing FASTA file
```

### Step 2: Run BLASTx Annotation

```bash
python step1_blastx_annotate.py NC_XXXXXX.fasta \
    --blast-db nr \
    --evalue 1e-10 \
    --threads 4
```

**Parameters:**
- `--blast-db` - BLAST database (nr, refseq_protein, or custom)
- `--evalue` - E-value cutoff (default: 1e-5, stricter: 1e-10)
- `--threads` - Parallel threads for BLAST

**What this does:**
1. Translates genome in all 6 reading frames
2. BLASTx search against protein database
3. Identifies ORFs with homology hits
4. Merges overlapping hits
5. Assigns gene names from best hits
6. Creates TSV for curation

**Output files:**
```
NC_XXXXXX_blastx.txt    # Raw BLAST results
NC_XXXXXX_blastx.tsv    # Curated annotation table
```

### Step 3: Manual Curation (Heavy)

**BLASTx results need careful review:**

```bash
nano NC_XXXXXX_blastx.tsv
```

**Common issues:**

1. **Overlapping ORFs** - BLASTx may predict multiple ORFs in same region
   - Choose the longest or best-scoring
   - Remove shorter overlapping hits

2. **Missed ORFs** - Small proteins may be missed
   - Manually check genome for short ORFs
   - Use ORF finder tools if needed

3. **Wrong reading frame** - Ensure correct frame
   - Check that start codon is ATG (or viral alternative)
   - Verify in-frame stop codon at end

4. **Gene name conflicts** - Multiple hits may suggest different names
   - Research virus family conventions
   - Choose most specific, accurate name

5. **False positives** - Non-coding regions may show weak hits
   - Remove hits with low identity (<30%)
   - Remove very short alignments (<100 aa)

### Step 4: Add to SnpEff

```bash
python step2_add_to_snpeff.py NC_XXXXXX NC_XXXXXX_blastx.tsv
```

---

## Pathway 4: Segmented Viruses

**Best for:** Multi-segment genomes (influenza, rotavirus, bunyaviruses)

**Examples:**
- **Influenza A/B:** 8 segments (PB2, PB1, PA, HA, NP, NA, M, NS)
- **Rotavirus:** 11 segments
- **Bunyavirus:** 3 segments (L, M, S)

### Approach: Combined Database

Instead of creating 8 separate databases, combine all segments into one:

```bash
python vicast_annotate_segmented.py influenza_h1n1_2009 \
    --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
    --names PB2,PB1,PA,HA,NP,NA,M,NS
```

**What this does:**
1. Downloads all 8 segments from NCBI
2. Combines sequences into single multi-FASTA
3. Combines annotations into single GFF3
4. Builds unified SnpEff database

**Output:**
```
$SNPEFF_DATA/influenza_h1n1_2009/
├── sequences.fa      # All 8 segments
├── genes.gff         # All annotations
└── snpEffectPredictor.bin
```

### Segment Naming

**Option 1: Use functional names (recommended)**

```bash
--names PB2,PB1,PA,HA,NP,NA,M,NS
```

VCF must use these names:
```vcf
#CHROM  POS   REF  ALT
PB2     100   A    G
HA      500   C    T
```

**Option 2: Use accessions**

```bash
# Don't provide --names flag
```

VCF must use accessions:
```vcf
#CHROM      POS   REF  ALT
CY121680    100   A    G
CY121683    500   C    T
```

### Using Local Files

If you already have segment files:

```bash
python vicast_annotate_segmented.py my_influenza \
    --fasta-files seg1.fa,seg2.fa,seg3.fa,seg4.fa,seg5.fa,seg6.fa,seg7.fa,seg8.fa \
    --gb-files seg1.gb,seg2.gb,seg3.gb,seg4.gb,seg5.gb,seg6.gb,seg7.gb,seg8.gb \
    --names PB2,PB1,PA,HA,NP,NA,M,NS
```

### Common Segmented Viruses

**Influenza A:**
```bash
python vicast_annotate_segmented.py my_influenza_a \
    --segments ACC1,ACC2,ACC3,ACC4,ACC5,ACC6,ACC7,ACC8 \
    --names PB2,PB1,PA,HA,NP,NA,M,NS
```

**Rotavirus:**
```bash
python vicast_annotate_segmented.py my_rotavirus \
    --segments ACC1,ACC2,ACC3,ACC4,ACC5,ACC6,ACC7,ACC8,ACC9,ACC10,ACC11 \
    --names VP1,VP2,VP3,VP4,NSP1,VP6,NSP3,NSP2,VP7,NSP4,NSP5
```

**Bunyavirus:**
```bash
python vicast_annotate_segmented.py my_bunyavirus \
    --segments ACC1,ACC2,ACC3 \
    --names L,M,S
```

---

## Manual Curation Guidelines

All pathways include manual curation checkpoints. Follow these best practices:

### 1. Verify Gene Boundaries

**Check start codons:**
- Standard: ATG
- Viral alternatives: CTG, GTG (some viruses)

**Check stop codons:**
- TAA, TAG, TGA

**Verify in-frame:**
```python
# Length should be divisible by 3
length = end - start + 1
assert length % 3 == 0, "Frame shift detected!"
```

### 2. Use Standard Nomenclature

**Follow virus family conventions:**

| Virus Family | Naming Convention | Example |
|--------------|-------------------|---------|
| Flavivirus | Protein_X or X | `Protein_E`, `NS5` |
| Coronavirus | Gene names | `S`, `N`, `ORF1ab` |
| HIV | Standard genes | `gag`, `pol`, `env` |
| Influenza | Segment name | `HA`, `NA`, `PB2` |

**Resources:**
- NCBI Virus database
- ViralZone (ExPASy)
- Family-specific nomenclature papers

### 3. Handle Polyproteins

**Decision: Keep or skip?**

**Skip polyproteins (default):**
- Individual mature peptides are more useful for annotation
- Pathway 2 automatically skips with `_no_polyprotein.tsv`

**Keep polyproteins (if needed):**
```bash
python step1_parse_viral_genome.py NC_001477 --keep-polyprotein
```

### 4. Handle Special Features

**Programmed frameshifts:**
- Some viruses use ribosomal frameshifting
- Annotate as two separate features: `gene_FS1` and `gene_FS2`
- Document in product description

**Ribosomal slippage:**
- Common in HIV gag-pol, coronavirus ORF1ab
- Can annotate as single feature or split

**Overlapping genes:**
- Some viruses have legitimate overlapping ORFs
- Keep both if both are functional
- Document overlap in description

### 5. Remove Artifacts

**Remove these:**
- Duplicate features (same coordinates)
- Generic hypothetical proteins (if better annotation available)
- False positive ORFs from BLASTx (very short, low identity)
- Non-coding RNA mistakenly annotated as CDS

---

## Troubleshooting

### Issue: "Genome not found in NCBI"

**Cause:** Accession doesn't exist or is incorrect

**Solution:**
```bash
# Verify accession at NCBI
# https://www.ncbi.nlm.nih.gov/nuccore/NC_001477

# Check for version number
python step1_parse_viral_genome.py NC_001477.1  # include version
```

### Issue: "No CDS features found"

**Cause:** Genome lacks annotations in GenBank

**Solution:**
```bash
# Use Pathway 3 (BLASTx)
python step1_blastx_annotate.py NC_XXXXXX.fasta --blast-db nr
```

### Issue: "BLASTx finds no hits"

**Cause:** Novel virus or wrong database

**Solutions:**
1. **Relax E-value:**
   ```bash
   --evalue 1e-3  # instead of default 1e-5
   ```

2. **Use viral-specific database:**
   ```bash
   --blast-db /path/to/viral_proteins
   ```

3. **Try different databases:**
   ```bash
   --blast-db refseq_protein  # faster than nr
   ```

### Issue: "SnpEff build failed"

**Cause:** GFF3 format errors or Java version

**Solution:**
```bash
# Check Java version (need 21+)
java -version

# Install correct Java
conda install -c conda-forge openjdk=21

# Verify GFF3 format
grep -v "^#" NC_001477.gff3 | head

# Check SnpEff logs
cat snpEff_build.log
```

### Issue: "Polyprotein cleavage sites unknown"

**Cause:** Polyprotein not annotated with mature peptides

**Solution:**
1. Search literature for cleavage sites
2. Use Pathway 3 BLASTx to identify regions
3. Manually annotate mature peptides in TSV
4. Or keep polyprotein and skip individual peptides

### Issue: "Database builds but annotation fails"

**Cause:** Chromosome name mismatch

**Solution:**
```bash
# VCF CHROM field must match FASTA sequence ID
# Check FASTA:
grep ">" NC_001477.fasta
# Output: >NC_001477

# VCF must use:
# CHROM = NC_001477 (not NC_001477.1 or just "1")
```

---

## Examples for Different Genome Types

### Example 1: Simple Virus (Poliovirus)

**Characteristics:** Single ORF, polyprotein
**Strategy:** Pathway 2, skip polyprotein, annotate mature peptides

```bash
python vicast_annotate.py NC_002058
# Auto-detects Pathway 2
# Skips polyprotein
# Creates annotation for mature proteins: VP1, VP2, VP3, VP4, 2A, 2B, 2C, 3A, 3B, 3C, 3D
```

### Example 2: Segmented Virus (Influenza)

**Characteristics:** 8 segments
**Strategy:** Pathway 4, combine segments

```bash
python vicast_annotate_segmented.py influenza_h1n1_2009 \
    --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
    --names PB2,PB1,PA,HA,NP,NA,M,NS
```

### Example 3: Complex Virus (HIV)

**Characteristics:** Multiple overlapping ORFs, frameshifts
**Strategy:** Pathway 2, manual curation of overlaps

```bash
python step1_parse_viral_genome.py NC_001802

# Manual curation in TSV:
# - Separate gag and pol (frameshift)
# - Annotate tat, rev (spliced, complex)
# - Keep all overlapping ORFs

python step2_add_to_snpeff.py NC_001802 NC_001802_curated.tsv
```

### Example 4: Poorly Annotated Virus

**Characteristics:** Novel or poorly studied
**Strategy:** Pathway 3, BLASTx with curation

```bash
# Download FASTA
# Run BLASTx against viral database
python step1_blastx_annotate.py novel_virus.fasta \
    --blast-db /ref/viral_proteins \
    --evalue 1e-10

# Heavy manual curation needed
# Research homologous viruses
# Refine gene boundaries
# Verify reading frames

python step2_add_to_snpeff.py novel_virus novel_virus_blastx_curated.tsv
```

---

## Best Practices

### For Publication-Quality Annotations

1. **Document curation decisions**
   - Keep notes on changes made
   - Reference papers for nomenclature
   - Explain non-standard features

2. **Validate against reference strains**
   - Compare to well-characterized isolates
   - Check gene order is conserved
   - Verify protein lengths match homologs

3. **Test with known variants**
   - Run SnpEff on reference VCF if available
   - Check that annotations make biological sense
   - Verify synonymous vs. non-synonymous calls

4. **Share your database**
   - Consider submitting to SnpEff repository
   - Include in supplementary materials
   - Document in methods section

### Optimization Tips

**Faster BLAST:**
```bash
# Use threads
--threads 8

# Use smaller database
--blast-db refseq_protein  # instead of nr

# Pre-filter by taxonomy
makeblastdb -in viral_proteins.fa -dbtype prot -taxid 10239
```

**Memory efficiency:**
```bash
# For large genomes or many segments
# Process in batches
# Use HPC job submission
```

---

## Quick Reference

### Common Commands

```bash
# Check if genome exists in SnpEff
python step0_check_snpeff.py NC_001477

# Pathway 2: Standard annotation
python step1_parse_viral_genome.py NC_001477
python step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv

# Pathway 3: BLASTx annotation
python step1_blastx_annotate.py genome.fasta --blast-db nr
python step2_add_to_snpeff.py genome genome_blastx.tsv

# Pathway 4: Segmented virus
python vicast_annotate_segmented.py virus_name \
    --segments ACC1,ACC2,ACC3 \
    --names seg1,seg2,seg3

# Verify database
java -jar $SNPEFF_JAR databases | grep NC_001477
```

### File Outputs

| File | Description |
|------|-------------|
| `{acc}.fasta` | Genome sequence |
| `{acc}_no_polyprotein.tsv` | Editable annotation table |
| `{acc}.gff3` | GFF3 annotation file |
| `$SNPEFF_DATA/{acc}/` | SnpEff database directory |

---

## Next Steps

After creating your SnpEff database:

1. **[Variant Calling Guide](VICAST_ANALYZE_GUIDE.md)** - Use your database for variant annotation
2. **Test with sample data** - Verify annotations are correct
3. **Document in methods** - Describe curation process for publications

---

**Need Help?**
- [Troubleshooting Guide](TROUBLESHOOTING.md)
- [GitHub Issues](https://github.com/mihinduk/VICAST/issues)
- Check existing annotations in SnpEff for similar viruses

---

**Last Updated:** 2026-02-05
**VICAST Version:** 2.2.0
