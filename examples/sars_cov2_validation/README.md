# SARS-CoV-2 Polyprotein Validation Dataset

This dataset validates VICAST's ability to handle **polyprotein annotation** using SARS-CoV-2 as a test case.

## Why This Matters

SARS-CoV-2 ORF1ab encodes a large polyprotein (~7,000 amino acids) that is cleaved into **16 mature peptides** (nsp1-16). Proper annotation of variants within these mature peptides is critical for:

- Identifying drug resistance mutations (e.g., in nsp5/3CLpro or nsp12/RdRp)
- Understanding viral fitness
- Tracking evolutionary changes

VICAST correctly annotates the ORF1ab polyprotein including:
- Ribosomal frameshift site
- All 16 mature peptides as distinct features
- Proper coordinate mapping for variants

## Dataset Contents

| File | Description |
|------|-------------|
| `NC_045512.gb` | Original GenBank file from NCBI RefSeq |
| `NC_045512.fasta` | Reference genome sequence |
| `NC_045512.gff3` | GFF3 annotation with polyprotein/mature peptides |
| `test_variants.vcf` | Test VCF with known SARS-CoV-2 mutations |
| `genbank_to_gff3_polyprotein.py` | GenBank to GFF3 converter (handles mature peptides) |
| `setup_and_validate.sh` | One-command setup and validation script |

## Test Variants

The test VCF includes clinically significant SARS-CoV-2 mutations:

| Position | Mutation | Gene | AA Change | Significance |
|----------|----------|------|-----------|--------------|
| 23403 | A>G | S (Spike) | D614G | Early fitness mutation, now fixed |
| 23063 | A>T | S (Spike) | N501Y | Receptor binding, Alpha/Omicron |
| 23012 | G>A | S (Spike) | E484K | Immune escape, Beta/Gamma |
| 28881 | G>A | N | R203K | Common in Alpha and later variants |
| 10500 | G>A | ORF1ab | G3412D | nsp5 (3CLpro) region |
| 15000 | T>C | ORF1ab | V4912A | nsp12 (RdRp) region |

## Genome Structure

| Gene | Position | Product | Notes |
|------|----------|---------|-------|
| ORF1ab | 266-21555 | Polyprotein | 16 mature peptides (nsp1-16) |
| S | 21563-25384 | Spike glycoprotein | Receptor binding |
| ORF3a | 25393-26220 | Accessory protein | |
| E | 26245-26472 | Envelope protein | |
| M | 26523-27191 | Membrane glycoprotein | |
| ORF6 | 27202-27387 | Accessory protein | |
| ORF7a | 27394-27759 | Accessory protein | |
| ORF7b | 27756-27887 | Accessory protein | |
| ORF8 | 27894-28259 | Accessory protein | |
| N | 28274-29533 | Nucleocapsid phosphoprotein | |
| ORF10 | 29558-29674 | Accessory protein | |

## ORF1ab Mature Peptides

The polyprotein is cleaved into 16 non-structural proteins:

| nsp | Position | Product | Function |
|-----|----------|---------|----------|
| nsp1 | 266-805 | Leader protein | Host shutoff |
| nsp2 | 806-2719 | nsp2 | Unknown |
| nsp3 | 2720-8554 | Papain-like protease | Polyprotein processing |
| nsp4 | 8555-10054 | nsp4 | Membrane scaffold |
| nsp5 | 10055-10972 | 3C-like proteinase | Main protease (drug target) |
| nsp6 | 10973-11842 | nsp6 | Autophagy inhibition |
| nsp7 | 11843-12091 | nsp7 | RdRp cofactor |
| nsp8 | 12092-12685 | nsp8 | Primase |
| nsp9 | 12686-13024 | nsp9 | RNA binding |
| nsp10 | 13025-13441 | nsp10 | Cofactor |
| nsp11 | 13442-13480 | nsp11 | Short peptide |
| nsp12 | 13442-16236 | RdRp | RNA-dependent RNA polymerase (drug target) |
| nsp13 | 16237-18039 | Helicase | RNA unwinding |
| nsp14 | 18040-19620 | ExoN | Proofreading |
| nsp15 | 19621-20658 | EndoRNase | Immune evasion |
| nsp16 | 20659-21552 | 2'-O-MTase | Cap methylation |

## Running Validation

```bash
cd examples/sars_cov2_validation
./setup_and_validate.sh
```

This script will:
1. Install SnpEff (if needed)
2. Validate VCF REF bases against reference
3. Build the SnpEff database
4. Annotate test variants
5. Verify all expected mutations are correctly annotated

## Building SnpEff Database Manually

```bash
# 1. Add genome to snpEff.config
echo "sars_cov2_wuhan.genome : SARS-CoV-2 Wuhan-Hu-1" >> $SNPEFF_HOME/snpEff.config

# 2. Create data directory
mkdir -p $SNPEFF_HOME/data/sars_cov2_wuhan

# 3. Copy files
cp NC_045512.fasta $SNPEFF_HOME/data/sars_cov2_wuhan/sequences.fa
cp NC_045512.gff3 $SNPEFF_HOME/data/sars_cov2_wuhan/genes.gff

# 4. Build database
java -jar $SNPEFF_HOME/snpEff.jar build -gff3 -v sars_cov2_wuhan
```

## Expected Output

After annotation, variants should show:
- Spike mutations with correct amino acid changes (D614G, N501Y, E484K)
- N protein mutations (R203K)
- ORF1ab mutations mapped to polyprotein coordinates

## Data Source

Reference genome downloaded from NCBI RefSeq:
- Strain: SARS-CoV-2 Wuhan-Hu-1
- Accession: NC_045512.2
- Downloaded: February 2026

## Key Validation Points

1. **Polyprotein annotation**: ORF1ab correctly represented with frameshift
2. **Mature peptides**: All 16 nsps annotated in GFF3
3. **Spike mutations**: Known variants correctly positioned
4. **VCF validation**: REF bases match reference genome

## Comparison with Other Tools

| Feature | VICAST | VADR | Other |
|---------|--------|------|-------|
| Polyprotein/mature peptides | Yes | Partial | Varies |
| Custom database building | Yes | Model-based | Varies |
| VCF REF validation | Yes | No | Rare |
| SnpEff integration | Native | Requires conversion | Varies |

## Citation

If using this dataset, please cite:
- VICAST: Mihindukulasuriya KA, Handley SA. GitHub: https://github.com/shandley/VICAST
- SARS-CoV-2 reference: Wu et al. (2020). Nature 579(7798):265-269
