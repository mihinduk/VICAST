# Dengue Virus 2 Validation Dataset

This dataset validates VICAST's ability to annotate **flavivirus genomes** with per-gene resolution across structural and nonstructural proteins.

## Why This Matters

Dengue virus 2 encodes a single polyprotein (~3,390 amino acids) that is co- and post-translationally cleaved into 14 mature proteins. Accurate per-gene annotation is critical for:

- Tracking mutations in the E protein (neutralization, tropism, virulence)
- Monitoring changes in NS1 (diagnostic target, complement activation)
- Identifying resistance mutations in NS3 (protease) and NS5 (RdRp) drug targets

VICAST resolves the polyprotein into individual genes, enabling variant annotation at the per-protein level rather than reporting all mutations relative to a single polyprotein ORF.

## Dataset Contents

| File | Description |
|------|-------------|
| `NC_001474.gb` | Original GenBank file from NCBI RefSeq |
| `NC_001474.fasta` | Reference genome sequence (10,723 bp) |
| `NC_001474.gff3` | GFF3 annotation with per-gene CDS features |
| `test_variants.vcf` | Test VCF with variants in key functional domains |
| `setup_and_validate.sh` | One-command setup and validation script |

## Test Variants

The test VCF includes variants targeting functionally important domains in structural and nonstructural proteins:

| Position | Mutation | Gene | AA Change | Functional Region |
|----------|----------|------|-----------|-------------------|
| 1231 | A>G | E | R99G | Fusion loop (neutralization determinant) |
| 1834 | T>C | E | S300P | Domain III (primary neutralizing epitope) |
| 2800 | G>A | NS1 | E127K | Diagnostic antigen, complement activation |
| 4700 | A>G | NS3 | H60R | Serine protease domain (drug target) |
| 5500 | C>A | NS3 | Q327K | Helicase domain (NTPase activity) |
| 9000 | G>A | NS5 | M477I | RNA-dependent RNA polymerase (drug target) |

## Genome Structure

| Gene | Position | Size (aa) | Function |
|------|----------|-----------|----------|
| C | 97-396 | 100 | Capsid protein |
| prM | 439-936 | 166 | Pre-membrane protein |
| E | 937-2421 | 495 | Envelope protein (receptor binding, fusion) |
| NS1 | 2422-3477 | 352 | Secreted glycoprotein (immune evasion, diagnostics) |
| NS2A | 3478-4131 | 218 | Membrane protein (assembly) |
| NS2B | 4132-4521 | 130 | NS3 protease cofactor |
| NS3 | 4522-6375 | 618 | Serine protease + RNA helicase |
| NS4A | 6376-6756 | 127 | Membrane rearrangement |
| 2K | 6757-6825 | 23 | Signal peptide |
| NS4B | 6826-7569 | 248 | Replication complex scaffold |
| NS5 | 7570-10269 | 900 | Methyltransferase + RNA-dependent RNA polymerase |

## Running Validation

```bash
cd examples/dengue2_validation
./setup_and_validate.sh
```

This script will:
1. Install SnpEff (if needed)
2. Validate VCF REF bases against reference
3. Build the SnpEff database with per-gene annotations
4. Annotate test variants
5. Verify all variants are correctly assigned to their respective genes

## Data Source

Reference genome downloaded from NCBI RefSeq:
- Strain: Dengue virus 2, New Guinea C (NGC)
- Accession: NC_001474.2
- Downloaded: February 2026

## Key Features Demonstrated

- **Flavivirus polyprotein resolution**: Single ORF correctly split into 14 mature proteins
- **Structural protein annotation**: C, prM, E with domain-level variant context
- **Nonstructural protein annotation**: NS1-NS5 with functional domain awareness
- **VCF REF base validation**: Prevents silent annotation errors
