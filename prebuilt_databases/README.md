# VICAST Pre-built SnpEff Databases

This directory contains pre-built SnpEff databases for common viral genomes, ready to use with VICAST.

## Available Databases

See `manifest.json` for the complete list of available databases.

Currently available:
- **DENV-2** (NC_001474.2) - Dengue virus type 2
- More coming soon!

## For Users: Installing Databases

### Using the installer script (recommended)

```bash
# List available databases
install_prebuilt_database.sh --list

# Install specific database
install_prebuilt_database.sh --install NC_001474.2

# Install all dengue databases
install_prebuilt_database.sh --tag dengue

# See what's installed
install_prebuilt_database.sh --installed
```

### Manual installation

1. Download the `.tar.gz` file for your genome
2. Extract to your SnpEff data directory:
   ```bash
   tar -xzf NC_001474.2.tar.gz -C $SNPEFF_DATA/
   ```
3. Verify installation:
   ```bash
   snpeff dump NC_001474.2 | head
   ```

## For Contributors: Adding New Databases

If you've built a SnpEff database for a viral genome and want to contribute it:

### 1. Package your database

```bash
# From VICAST root directory
bash prebuilt_databases/package_database.sh NC_001474.2

# This creates:
#   - NC_001474.2.tar.gz  (packaged database)
#   - NC_001474.2.md5     (checksum)
```

### 2. Upload to GitHub

Add the `.tar.gz` file to this directory and commit:

```bash
git add prebuilt_databases/NC_001474.2.tar.gz
git commit -m "Add pre-built database for DENV-2 (NC_001474.2)"
```

### 3. Update manifest.json

Add an entry for your database:

```json
{
  "accession": "NC_001474.2",
  "name": "DENV-2",
  "description": "Dengue virus type 2",
  "size_mb": 2.1,
  "features": 10,
  "tags": ["dengue", "flavivirus", "arbovirus"],
  "url": "https://github.com/shandley/VICAST/raw/main/prebuilt_databases/NC_001474.2.tar.gz",
  "md5": "abc123def456..."
}
```

**Fields:**
- `accession`: GenBank accession (genome identifier)
- `name`: Short name (e.g., "DENV-2", "ZIKV")
- `description`: Full organism name
- `size_mb`: Package size in MB
- `features`: Number of gene features in annotation
- `tags`: Categories for filtering (e.g., ["dengue", "flavivirus"])
- `url`: Direct download URL
- `md5`: MD5 checksum from `.md5` file

### 4. Submit Pull Request

Submit your changes via GitHub pull request with:
- The packaged `.tar.gz` file
- Updated `manifest.json`
- Brief description of the genome

## Inspecting Pre-built Databases

Before using a pre-built database, you can inspect its annotation quality:

```bash
# View all genes/features in the database
snpeff dump NC_001474.2 | head -50

# See the GFF3 annotation
cat $SNPEFF_DATA/NC_001474.2/genes.gff

# Count how many features are annotated
grep -c "CDS" $SNPEFF_DATA/NC_001474.2/genes.gff

# Check protein sequences
head $SNPEFF_DATA/NC_001474.2/protein.fa
```

## Overwriting Pre-built Databases

**Don't like our annotation?** You can always replace it with your own:

### Option 1: Build from scratch (overwrites existing)

```bash
# Download genome and build custom database
step1_parse_viral_genome.py NC_001474.2

# This creates NC_001474.2_features.tsv - edit it to your liking
# Then build the database (overwrites pre-built)
step2_add_to_snpeff.py NC_001474.2 NC_001474.2_features.tsv --force
```

### Option 2: Reinstall pre-built (revert changes)

```bash
# Reinstall original pre-built database
install_prebuilt_database.sh --install NC_001474.2 --force
```

**Key points:**
- ✅ Pre-built databases are **starting points**, not requirements
- ✅ You have **full control** to modify or replace them
- ✅ Custom annotations always **overwrite** pre-built ones
- ✅ Use `--force` flag to overwrite existing databases

## Building from Scratch

If you want to build your own database instead of using pre-built ones:

```bash
# Download genome from NCBI and build database
step1_parse_viral_genome.py NC_001474.2
step2_add_to_snpeff.py NC_001474.2 NC_001474.2_features.tsv

# Then package it to share
bash prebuilt_databases/package_database.sh NC_001474.2
```

## Database Quality Standards

Pre-built databases should meet these criteria:

- ✅ Complete genome annotation (all CDSs annotated)
- ✅ Verified with `snpeff dump` command
- ✅ Tested with at least one VCF file
- ✅ Clear documentation of source (GenBank accession)
- ✅ Tagged appropriately (virus family, type)

## Priority Genomes

We welcome contributions of these commonly-used viral genomes:

**Flaviviruses:**
- DENV-1, DENV-3, DENV-4 (DENV-2 already available)
- Zika virus (ZIKV)
- Yellow fever virus (YFV)
- West Nile virus (WNV)
- Japanese encephalitis virus (JEV)

**Other Arboviruses:**
- Chikungunya virus (CHIKV)
- Ross River virus (RRV)
- Venezuelan equine encephalitis virus (VEEV)

**Respiratory Viruses:**
- SARS-CoV-2
- Influenza A (H1N1, H3N2)
- RSV

**Others:**
- Norovirus
- Enterovirus
- Poliovirus

## File Size Guidelines

- Keep individual databases under 10 MB
- Use `tar -czf` for maximum compression
- Large genome families (e.g., all Influenza segments) can be bundled

## Contact

Questions about contributing databases? Open an issue on GitHub or contact the VICAST maintainers.
