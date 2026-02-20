#!/usr/bin/env python3
"""
Generate consensus genome and per-gene protein sequences.

Tier 1 of the two-tier VICAST post-processing pipeline:
  - Applies all high-AF variants to the reference genome
  - AF threshold is genome-type-aware:
      ssRNA/ssDNA >= 0.95 (haploid)
      dsRNA       >= 0.45 (diploid-like)
  - Translates per-gene CDS regions from the consensus
  - Produces a single consensus genome + protein FASTA + report
"""

import argparse
import sys
import os
import json
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from viral_translator import viral_translate

# ── Default AF thresholds by genome type ──────────────────────────────
GENOME_TYPE_THRESHOLDS = {
    'ssRNA': 0.95,
    'dsRNA': 0.45,
    'ssDNA': 0.95,
    'dsDNA': 0.95,
    'unknown': 0.95,
}


def read_virus_config(accession):
    """Read virus configuration from known_viruses.json."""
    search_paths = [
        Path(__file__).parent / "known_viruses.json",
        Path(__file__).parent.parent / "virus_configs" / "known_viruses.json",
        Path(__file__).parent.parent / "visualization" / "known_viruses.json",
        Path("known_viruses.json"),
    ]
    for path in search_paths:
        if path.exists():
            with open(path, 'r') as f:
                viruses = json.load(f)
            return viruses.get(accession)
    return None


def get_consensus_threshold(virus_config, genome_type_override=None, af_override=None):
    """Determine the consensus AF threshold.

    Priority: explicit --consensus-af > --genome-type override > config genome_type > default.
    """
    if af_override is not None:
        return af_override

    genome_type = genome_type_override
    if genome_type is None and virus_config:
        genome_type = virus_config.get('genome_type', 'unknown')
    if genome_type is None:
        genome_type = 'unknown'

    # Normalise: accept 'ss', 'ds' shorthand
    gt = genome_type.lower().replace('-', '').replace('_', '')
    if gt in ('ss', 'ssrna'):
        return GENOME_TYPE_THRESHOLDS['ssRNA']
    elif gt in ('ds', 'dsrna'):
        return GENOME_TYPE_THRESHOLDS['dsRNA']
    elif gt in ('ssdna',):
        return GENOME_TYPE_THRESHOLDS['ssDNA']
    elif gt in ('dsdna',):
        return GENOME_TYPE_THRESHOLDS['dsDNA']

    return GENOME_TYPE_THRESHOLDS.get(genome_type, 0.95)


def load_filtered_tsv(tsv_path):
    """Load the filtered mutations TSV produced by parse_snpeff_tsv.py."""
    df = pd.read_csv(tsv_path, sep='\t')
    return df


def apply_variants_to_sequence(ref_seq, variants_df):
    """Apply SNPs and simple substitutions to reference sequence.

    Returns (consensus_sequence, list_of_applied_mutations).
    Indels are recorded but only simple SNPs are applied for now.
    """
    seq_list = list(ref_seq)
    applied = []
    skipped = []

    # Sort by position descending for safe indel handling (future)
    for _, row in variants_df.sort_values('POS').iterrows():
        pos = int(row['POS']) - 1  # 0-based
        ref_base = str(row['REF'])
        alt_base = str(row['ALT'])

        # SNP
        if len(ref_base) == 1 and len(alt_base) == 1:
            if pos < len(seq_list):
                if seq_list[pos].upper() == ref_base.upper():
                    seq_list[pos] = alt_base
                    applied.append(row)
                else:
                    # Mismatch — apply anyway but warn
                    seq_list[pos] = alt_base
                    applied.append(row)
                    print(f"  Warning: REF mismatch at {pos+1}: expected {ref_base}, found {seq_list[pos]}")
        else:
            # Indel — record but skip for now
            skipped.append(row)
            print(f"  Skipping indel at {pos+1}: {ref_base}>{alt_base}")

    return ''.join(seq_list), applied, skipped


def translate_genes(consensus_seq, gene_coords):
    """Translate per-gene CDS from consensus sequence.

    Returns list of (gene_name, protein_seq) tuples.
    """
    proteins = []
    for gene, (start, end) in gene_coords.items():
        gene_seq = consensus_seq[start - 1:end]
        protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
        proteins.append((gene, protein))
    return proteins


def write_report(report_path, args, virus_config, threshold, variants_df, applied, skipped, proteins):
    """Write a human-readable consensus report."""
    genome_type = 'unknown'
    if virus_config:
        genome_type = virus_config.get('genome_type', 'unknown')
    if args.genome_type:
        genome_type = args.genome_type

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("VICAST CONSENSUS GENOME REPORT\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Accession:          {args.accession}\n")
        if virus_config:
            f.write(f"Virus name:         {virus_config.get('name', 'N/A')}\n")
        f.write(f"Genome type:        {genome_type}\n")
        f.write(f"Consensus AF:       >= {threshold}\n")
        f.write(f"Input TSV:          {args.vcf}\n")
        f.write(f"Reference FASTA:    {args.reference}\n\n")

        f.write(f"Total variants in TSV:       {len(variants_df)}\n")
        consensus_variants = variants_df[variants_df['Allele_Frequency'] >= threshold]
        f.write(f"Variants >= AF threshold:    {len(consensus_variants)}\n")
        f.write(f"Variants applied (SNPs):     {len(applied)}\n")
        f.write(f"Variants skipped (indels):   {len(skipped)}\n\n")

        f.write("-" * 80 + "\n")
        f.write("APPLIED MUTATIONS\n")
        f.write("-" * 80 + "\n")
        if applied:
            for row in applied:
                gene = row.get('GENE_NAME', '')
                hgvsp = row.get('HGVSp', '')
                effect = row.get('EFFECT', '')
                af = row['Allele_Frequency']
                f.write(f"  {row['POS']} {row['REF']}>{row['ALT']}  AF={af:.4f}  "
                        f"gene={gene}  effect={effect}  HGVSp={hgvsp}\n")
        else:
            f.write("  (none)\n")

        f.write("\n")
        f.write("-" * 80 + "\n")
        f.write("TRANSLATED PROTEINS\n")
        f.write("-" * 80 + "\n")
        for gene_name, protein in proteins:
            f.write(f"  {gene_name}: {len(protein)} aa\n")

    print(f"Report written to: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate consensus genome and proteins from high-AF variants')
    parser.add_argument('--vcf', required=True,
                        help='Filtered TSV from parse_snpeff_tsv.py')
    parser.add_argument('--reference', required=True,
                        help='Reference genome FASTA')
    parser.add_argument('--accession', required=True,
                        help='Virus accession')
    parser.add_argument('--consensus-af', type=float, default=None,
                        help='Override consensus AF threshold (default: auto from genome_type)')
    parser.add_argument('--genome-type', default=None,
                        help='Override genome type: ss|ds|ssRNA|dsRNA (default: auto from known_viruses.json)')
    parser.add_argument('--output-prefix', required=True,
                        help='Output prefix for consensus files')
    args = parser.parse_args()

    print("=" * 80)
    print("VICAST CONSENSUS GENOME GENERATION (Tier 1)")
    print("=" * 80)

    # ── Load reference ────────────────────────────────────────────────
    ref_record = next(SeqIO.parse(args.reference, "fasta"))
    ref_seq = str(ref_record.seq)
    print(f"Reference genome: {len(ref_seq)} bp")

    # ── Load virus config ─────────────────────────────────────────────
    virus_config = read_virus_config(args.accession)
    if virus_config:
        print(f"Virus config: {virus_config.get('name', 'N/A')} "
              f"({virus_config.get('genome_type', 'unknown')})")
    else:
        print(f"Warning: No config for {args.accession} — will generate polyprotein only")

    # ── Determine AF threshold ────────────────────────────────────────
    threshold = get_consensus_threshold(virus_config, args.genome_type, args.consensus_af)
    print(f"Consensus AF threshold: >= {threshold}")

    # ── Load variants ─────────────────────────────────────────────────
    variants_df = load_filtered_tsv(args.vcf)
    print(f"Loaded {len(variants_df)} filtered variants")

    consensus_variants = variants_df[variants_df['Allele_Frequency'] >= threshold].copy()
    print(f"Variants at AF >= {threshold}: {len(consensus_variants)}")

    # ── Build consensus ───────────────────────────────────────────────
    if len(consensus_variants) == 0:
        print("No variants above consensus threshold — consensus equals reference")
        consensus_seq = ref_seq
        applied = []
        skipped = []
    else:
        consensus_seq, applied, skipped = apply_variants_to_sequence(ref_seq, consensus_variants)
        print(f"Applied {len(applied)} mutations, skipped {len(skipped)} indels")

    # ── Write consensus FASTA ─────────────────────────────────────────
    sample_name = Path(args.output_prefix).name
    consensus_fasta = f"{args.output_prefix}_consensus.fasta"

    desc = f"consensus genome ({len(applied)} mutations applied, AF >= {threshold})"
    record = SeqRecord(Seq(consensus_seq), id=f"{sample_name}", description=desc)
    SeqIO.write([record], consensus_fasta, "fasta")
    print(f"Consensus genome written to: {consensus_fasta}")

    # ── Translate proteins ────────────────────────────────────────────
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    if gene_coords:
        proteins = translate_genes(consensus_seq, gene_coords)
        print(f"Translated {len(proteins)} gene products")
    else:
        protein = viral_translate(consensus_seq, coordinates=None, stop_at_stop_codon=True)
        proteins = [('Polyprotein', protein)]
        print("No gene coordinates — generated single polyprotein")

    # ── Write protein FASTA ───────────────────────────────────────────
    protein_fasta = f"{args.output_prefix}_consensus_proteins.fasta"
    protein_records = []
    for gene_name, protein_seq in proteins:
        rec = SeqRecord(
            Seq(protein_seq),
            id=f"{sample_name}_{gene_name}",
            description=f"{gene_name} ({len(protein_seq)} aa)"
        )
        protein_records.append(rec)
    SeqIO.write(protein_records, protein_fasta, "fasta")
    print(f"Consensus proteins written to: {protein_fasta}")

    # ── Write report ──────────────────────────────────────────────────
    report_path = f"{args.output_prefix}_consensus_report.txt"
    write_report(report_path, args, virus_config, threshold,
                 variants_df, applied, skipped, proteins)

    print("")
    print("=" * 80)
    print("Consensus generation complete")
    print("=" * 80)


if __name__ == "__main__":
    main()
