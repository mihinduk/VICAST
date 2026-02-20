#!/usr/bin/env python3
"""
Generate consensus genome and per-gene protein sequences.

Tier 1 of the two-tier VICAST post-processing pipeline:
  - Applies all high-AF variants to the reference genome
  - AF threshold is genome-type-aware:
      ssRNA/ssDNA >= 0.95 (haploid)
      dsRNA       >= 0.45 (diploid-like)
  - Masks positions below --min-depth with N
  - Translates per-gene CDS regions from the consensus (skips UTRs)
  - Produces a single consensus genome + protein FASTA + report
"""

import argparse
import sys
import os
import json
import re
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from viral_translator import viral_translate

import subprocess

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
    return pd.read_csv(tsv_path, sep='\t')


def apply_variants_to_sequence(ref_seq, variants_df):
    """Apply SNPs and simple substitutions to reference sequence.

    Returns (consensus_sequence, list_of_applied_mutation_rows, list_of_skipped_rows).
    """
    seq_list = list(ref_seq)
    applied = []
    skipped = []

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
                    seq_list[pos] = alt_base
                    applied.append(row)
                    print(f"  Warning: REF mismatch at {pos+1}: expected {ref_base}, "
                          f"found {seq_list[pos]}")
        else:
            skipped.append(row)
            print(f"  Skipping indel at {pos+1}: {ref_base}>{alt_base}")

    return ''.join(seq_list), applied, skipped


def _load_depth_from_file(depth_file, genome_len):
    """Parse a samtools-depth TSV (chrom, pos, depth) into an array."""
    depth_array = [0] * genome_len
    with open(depth_file, 'r') as f:
        for line in f:
            if line.startswith('chrom') or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                try:
                    pos = int(parts[1]) - 1  # 1-based → 0-based
                    depth = int(parts[2])
                    if 0 <= pos < genome_len:
                        depth_array[pos] = depth
                except ValueError:
                    continue
    return depth_array


def _load_depth_from_bam(bam_path, genome_len):
    """Run samtools depth -a on a BAM and return per-base depth array."""
    depth_array = [0] * genome_len
    try:
        result = subprocess.run(
            ['samtools', 'depth', '-a', str(bam_path)],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 3:
                pos = int(parts[1]) - 1
                depth = int(parts[2])
                if 0 <= pos < genome_len:
                    depth_array[pos] = depth
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  Warning: samtools depth failed ({e}) — skipping masking")
        return None
    return depth_array


def mask_low_coverage(consensus_seq, min_depth, bam_path=None, depth_file=None):
    """Replace positions with coverage < min_depth with N.

    Uses depth_file if provided, otherwise runs samtools depth on the BAM.
    Returns (masked_seq, n_count, list of masked regions as (start, end) 1-based).
    """
    seq_list = list(consensus_seq)
    genome_len = len(seq_list)

    # Load depth from file or BAM
    if depth_file and Path(depth_file).exists():
        print(f"  Reading depth from: {depth_file}")
        depth_array = _load_depth_from_file(depth_file, genome_len)
    elif bam_path and Path(bam_path).exists():
        print(f"  Computing depth from BAM: {bam_path}")
        depth_array = _load_depth_from_bam(bam_path, genome_len)
        if depth_array is None:
            return consensus_seq, 0, []
    else:
        print("  Warning: No depth source available — skipping masking")
        return consensus_seq, 0, []

    # Mask low-coverage positions
    n_count = 0
    masked_regions = []
    in_region = False
    region_start = None

    for i in range(genome_len):
        if depth_array[i] < min_depth:
            seq_list[i] = 'N'
            n_count += 1
            if not in_region:
                region_start = i + 1  # 1-based
                in_region = True
        else:
            if in_region:
                masked_regions.append((region_start, i))  # end is 1-based
                in_region = False

    if in_region:
        masked_regions.append((region_start, genome_len))

    return ''.join(seq_list), n_count, masked_regions


def is_utr(gene_name):
    """Return True for UTR entries that should not be translated."""
    name = gene_name.lower().replace("'", "").replace("_", "")
    return 'utr' in name


def translate_genes(consensus_seq, gene_coords):
    """Translate per-gene CDS from consensus sequence (skips UTRs).

    Returns list of (gene_name, protein_seq) tuples.
    """
    proteins = []
    for gene, (start, end) in gene_coords.items():
        if is_utr(gene):
            continue
        gene_seq = consensus_seq[start - 1:end]
        protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
        proteins.append((gene, protein))
    return proteins


def get_mutations_for_gene(applied, gene_name, gene_start, gene_end):
    """Get list of mutations affecting a specific gene.

    Returns list of dicts with 'nt_change' and 'aa_change' keys.
    """
    mutations = []
    for row in applied:
        pos = int(row['POS'])
        if gene_start <= pos <= gene_end:
            nt_change = f"{pos}{row['REF']}>{row['ALT']}"
            # Parse HGVSp for amino acid change
            hgvsp = str(row.get('HGVSp', ''))
            aa_change = ''
            if hgvsp and hgvsp != 'nan':
                # Extract the p. notation
                match = re.search(r'(p\.\S+)', hgvsp)
                if match:
                    aa_change = match.group(1)
            mutations.append({'nt_change': nt_change, 'aa_change': aa_change})
    return mutations


def format_mutations_for_header(mutations, mode='nt'):
    """Format mutation list for FASTA header.

    mode='nt': nucleotide changes (e.g., 241C>T,3037C>T)
    mode='aa': amino acid changes (e.g., p.Asp614Gly,p.Asn501Tyr)
    mode='both': both nt and aa
    """
    if not mutations:
        return ''
    if mode == 'nt':
        return ','.join(m['nt_change'] for m in mutations)
    elif mode == 'aa':
        aa_parts = [m['aa_change'] for m in mutations if m['aa_change']]
        return ','.join(aa_parts) if aa_parts else ''
    elif mode == 'both':
        parts = []
        for m in mutations:
            if m['aa_change']:
                parts.append(f"{m['nt_change']}({m['aa_change']})")
            else:
                parts.append(m['nt_change'])
        return ','.join(parts)
    return ''


def write_report(report_path, args, virus_config, threshold, variants_df,
                 applied, skipped, proteins, n_count, masked_regions, min_depth):
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
        f.write(f"Min depth:          {min_depth}X\n")
        f.write(f"Input TSV:          {args.vcf}\n")
        f.write(f"Reference FASTA:    {args.reference}\n")
        if args.depth_file:
            f.write(f"Depth file:         {args.depth_file}\n")
        elif args.bam:
            f.write(f"BAM file:           {args.bam}\n")
        f.write("\n")

        f.write(f"Total variants in TSV:       {len(variants_df)}\n")
        consensus_variants = variants_df[variants_df['Allele_Frequency'] >= threshold]
        f.write(f"Variants >= AF threshold:    {len(consensus_variants)}\n")
        f.write(f"Variants applied (SNPs):     {len(applied)}\n")
        f.write(f"Variants skipped (indels):   {len(skipped)}\n")
        f.write(f"Positions masked (N):        {n_count}\n\n")

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

        if masked_regions:
            f.write("\n")
            f.write("-" * 80 + "\n")
            f.write(f"LOW COVERAGE REGIONS (< {min_depth}X) — masked with N\n")
            f.write("-" * 80 + "\n")
            for start, end in masked_regions:
                f.write(f"  {start}-{end} ({end - start + 1} bp)\n")

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
    parser.add_argument('--bam', default=None,
                        help='BAM file for depth masking (fallback if no --depth-file)')
    parser.add_argument('--depth-file', default=None,
                        help='Pre-computed samtools depth file (preferred over --bam)')
    parser.add_argument('--min-depth', type=int, default=20,
                        help='Minimum depth — positions below this are masked with N (default: 20)')
    parser.add_argument('--consensus-af', type=float, default=None,
                        help='Override consensus AF threshold (default: auto from genome_type)')
    parser.add_argument('--genome-type', default=None,
                        help='Override genome type: ss|ds|ssRNA|dsRNA')
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

    # ── Mask low-coverage positions ───────────────────────────────────
    n_count = 0
    masked_regions = []
    if args.depth_file or args.bam:
        print(f"Masking positions with depth < {args.min_depth}X ...")
        consensus_seq, n_count, masked_regions = mask_low_coverage(
            consensus_seq, args.min_depth,
            bam_path=args.bam, depth_file=args.depth_file)
        print(f"  Masked {n_count} positions with N ({len(masked_regions)} regions)")
    else:
        print("No BAM or depth file provided — skipping depth masking")

    # ── Build genome FASTA header ─────────────────────────────────────
    sample_name = Path(args.output_prefix).name
    consensus_fasta = f"{args.output_prefix}_consensus.fasta"

    # Header includes: mutation count, nt mutations, N count
    header_parts = [f"{len(applied)} mutations"]
    if applied:
        nt_muts = ','.join(f"{int(r['POS'])}{r['REF']}>{r['ALT']}" for r in applied)
        header_parts.append(f"variants={nt_muts}")
    header_parts.append(f"AF>={threshold}")
    if n_count > 0:
        header_parts.append(f"{n_count} Ns (depth<{args.min_depth}X)")
    desc = ' '.join(header_parts)

    record = SeqRecord(Seq(consensus_seq), id=f"{sample_name}", description=desc)
    SeqIO.write([record], consensus_fasta, "fasta")
    print(f"Consensus genome written to: {consensus_fasta}")

    # ── Translate proteins ────────────────────────────────────────────
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    if gene_coords:
        proteins = translate_genes(consensus_seq, gene_coords)
        print(f"Translated {len(proteins)} gene products (UTRs excluded)")
    else:
        protein = viral_translate(consensus_seq, coordinates=None, stop_at_stop_codon=True)
        proteins = [('Polyprotein', protein)]
        print("No gene coordinates — generated single polyprotein")

    # ── Write protein FASTA with mutation info in headers ─────────────
    protein_fasta = f"{args.output_prefix}_consensus_proteins.fasta"
    protein_records = []
    for gene_name, protein_seq in proteins:
        # Find mutations in this gene
        gene_start, gene_end = gene_coords.get(gene_name, (0, 0))
        gene_muts = get_mutations_for_gene(applied, gene_name, gene_start, gene_end)

        # Build header
        desc_parts = [f"{gene_name} ({len(protein_seq)} aa)"]
        if gene_muts:
            nt_str = format_mutations_for_header(gene_muts, mode='nt')
            aa_str = format_mutations_for_header(gene_muts, mode='aa')
            desc_parts.append(f"nt={nt_str}")
            if aa_str:
                desc_parts.append(f"aa={aa_str}")
        else:
            desc_parts.append("no mutations")

        rec = SeqRecord(
            Seq(protein_seq),
            id=f"{sample_name}_{gene_name}",
            description=' '.join(desc_parts)
        )
        protein_records.append(rec)
    SeqIO.write(protein_records, protein_fasta, "fasta")
    print(f"Consensus proteins written to: {protein_fasta}")

    # ── Write report ──────────────────────────────────────────────────
    report_path = f"{args.output_prefix}_consensus_report.txt"
    write_report(report_path, args, virus_config, threshold,
                 variants_df, applied, skipped, proteins,
                 n_count, masked_regions, args.min_depth)

    print("")
    print("=" * 80)
    print("Consensus generation complete")
    print("=" * 80)


if __name__ == "__main__":
    main()
