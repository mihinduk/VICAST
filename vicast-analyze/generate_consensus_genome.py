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
  - Supports multi-segment viruses (e.g. influenza): auto-detected from reference FASTA
"""

import argparse
import sys
import os
import json
import re
from pathlib import Path
from collections import defaultdict, OrderedDict

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
    # keep_default_na=False: preserve literal "NA" values (e.g. influenza NA segment)
    return pd.read_csv(tsv_path, sep='\t', keep_default_na=False)


def apply_variants_to_sequence(ref_seq, variants_df):
    """Apply SNPs and indels to reference sequence.

    Indels are applied in reverse position order to preserve upstream coordinates.
    Returns (consensus_sequence, applied_rows, skipped_rows, indel_offsets).
    indel_offsets is a sorted list of (1-based pos, size_change) for coordinate adjustment.
    """
    seq_list = list(ref_seq)
    applied = []
    skipped = []
    indel_offsets = []  # (pos_1based, size_change) for downstream coordinate adjustment

    # Sort by position descending so indels don't shift upstream coordinates
    sorted_variants = variants_df.sort_values('POS', ascending=False)

    for _, row in sorted_variants.iterrows():
        pos = int(row['POS']) - 1  # 0-based
        ref_base = str(row['REF'])
        alt_base = str(row['ALT'])

        if pos >= len(seq_list):
            skipped.append(row)
            continue

        if len(ref_base) == 1 and len(alt_base) == 1:
            # SNP
            seq_list[pos] = alt_base
            applied.append(row)
        elif len(ref_base) < len(alt_base) and alt_base.startswith(ref_base):
            # Insertion: REF=A, ALT=ATCG → insert TCG after pos
            insert_seq = alt_base[len(ref_base):]
            insert_pos = pos + len(ref_base)
            seq_list[insert_pos:insert_pos] = list(insert_seq)
            indel_offsets.append((int(row['POS']), len(insert_seq)))
            applied.append(row)
            print(f"  Applied insertion at {int(row['POS'])}: +{len(insert_seq)}bp")
        elif len(ref_base) > len(alt_base) and ref_base.startswith(alt_base):
            # Deletion: REF=ATCG, ALT=A → delete TCG after anchor
            del_start = pos + len(alt_base)
            del_len = len(ref_base) - len(alt_base)
            del seq_list[del_start:del_start + del_len]
            indel_offsets.append((int(row['POS']), -del_len))
            applied.append(row)
            print(f"  Applied deletion at {int(row['POS'])}: -{del_len}bp")
        else:
            # Complex substitution (MNV) — apply as replacement
            seq_list[pos:pos + len(ref_base)] = list(alt_base)
            size_change = len(alt_base) - len(ref_base)
            if size_change != 0:
                indel_offsets.append((int(row['POS']), size_change))
            applied.append(row)
            print(f"  Applied complex variant at {int(row['POS'])}: "
                  f"{ref_base}>{alt_base}")

    # Re-sort applied list by position (ascending) for reporting
    applied.sort(key=lambda r: int(r['POS']))
    # Sort offsets by position ascending for coordinate adjustment
    indel_offsets.sort(key=lambda x: x[0])

    return ''.join(seq_list), applied, skipped, indel_offsets


def adjust_gene_coords(gene_coords, indel_offsets):
    """Adjust gene coordinates based on indel offsets.

    Returns new dict of gene_name -> (adjusted_start, adjusted_end).
    """
    if not indel_offsets:
        return gene_coords

    adjusted = {}
    for gene, (start, end) in gene_coords.items():
        cum_offset = 0
        for indel_pos, size_change in indel_offsets:
            if indel_pos < start:
                cum_offset += size_change
            elif indel_pos <= end:
                # Indel within the gene — adjust end only
                end += size_change
        adjusted[gene] = (start + cum_offset, end + cum_offset)
    return adjusted


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


def _load_depth_from_file_segmented(depth_file, ref_records):
    """Parse a samtools-depth TSV into per-segment depth arrays.

    Returns dict[segment_id -> depth_array].
    """
    depth_by_seg = {}
    for seg_id, seg_seq in ref_records.items():
        depth_by_seg[seg_id] = [0] * len(seg_seq)

    with open(depth_file, 'r') as f:
        for line in f:
            if line.startswith('chrom') or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 3:
                try:
                    chrom = parts[0]
                    pos = int(parts[1]) - 1  # 1-based → 0-based
                    depth = int(parts[2])
                    if chrom in depth_by_seg and 0 <= pos < len(depth_by_seg[chrom]):
                        depth_by_seg[chrom][pos] = depth
                except ValueError:
                    continue
    return depth_by_seg


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


def _load_depth_from_bam_segmented(bam_path, ref_records):
    """Run samtools depth -a on a BAM and return per-segment depth arrays.

    Returns dict[segment_id -> depth_array], or None on failure.
    """
    depth_by_seg = {}
    for seg_id, seg_seq in ref_records.items():
        depth_by_seg[seg_id] = [0] * len(seg_seq)

    try:
        result = subprocess.run(
            ['samtools', 'depth', '-a', str(bam_path)],
            capture_output=True, text=True, check=True
        )
        for line in result.stdout.splitlines():
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                pos = int(parts[1]) - 1
                depth = int(parts[2])
                if chrom in depth_by_seg and 0 <= pos < len(depth_by_seg[chrom]):
                    depth_by_seg[chrom][pos] = depth
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"  Warning: samtools depth failed ({e}) — skipping masking")
        return None
    return depth_by_seg


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


def mask_low_coverage_from_array(consensus_seq, min_depth, depth_array):
    """Mask positions with coverage < min_depth using a pre-loaded depth array.

    Returns (masked_seq, n_count, list of masked regions as (start, end) 1-based).
    """
    seq_list = list(consensus_seq)
    genome_len = len(seq_list)
    n_count = 0
    masked_regions = []
    in_region = False
    region_start = None

    for i in range(genome_len):
        if i < len(depth_array) and depth_array[i] < min_depth:
            seq_list[i] = 'N'
            n_count += 1
            if not in_region:
                region_start = i + 1
                in_region = True
        else:
            if in_region:
                masked_regions.append((region_start, i))
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


def get_gene_coords_for_segment(virus_config, segment_id):
    """Get gene coordinates for a specific segment.

    For segmented viruses, returns gene_coords_by_segment[segment_id].
    For non-segmented viruses, returns flat gene_coords.
    """
    if not virus_config:
        return {}
    by_segment = virus_config.get('gene_coords_by_segment', {})
    if by_segment and segment_id in by_segment:
        return by_segment[segment_id]
    return virus_config.get('gene_coords', {})


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
        n_snps_applied = sum(1 for r in applied if len(str(r['REF'])) == 1 and len(str(r['ALT'])) == 1)
        n_indels_applied = len(applied) - n_snps_applied
        f.write(f"Variants applied:            {len(applied)} ({n_snps_applied} SNPs, {n_indels_applied} indels)\n")
        f.write(f"Variants skipped:            {len(skipped)}\n")
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


def write_report_segmented(report_path, args, virus_config, threshold, variants_df,
                           segment_results, min_depth):
    """Write a human-readable consensus report for segmented viruses."""
    genome_type = 'unknown'
    if virus_config:
        genome_type = virus_config.get('genome_type', 'unknown')
    if args.genome_type:
        genome_type = args.genome_type

    # Aggregate totals
    total_applied = sum(len(sr['applied']) for sr in segment_results.values())
    total_skipped = sum(len(sr['skipped']) for sr in segment_results.values())
    total_n_count = sum(sr['n_count'] for sr in segment_results.values())
    total_proteins = []
    for sr in segment_results.values():
        total_proteins.extend(sr['proteins'])

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("VICAST CONSENSUS GENOME REPORT (SEGMENTED)\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Accession:          {args.accession}\n")
        if virus_config:
            f.write(f"Virus name:         {virus_config.get('name', 'N/A')}\n")
        f.write(f"Genome type:        {genome_type}\n")
        f.write(f"Segments:           {len(segment_results)}\n")
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
        n_snps_total = sum(1 for sr in segment_results.values()
                          for r in sr['applied']
                          if len(str(r['REF'])) == 1 and len(str(r['ALT'])) == 1)
        n_indels_total = total_applied - n_snps_total
        f.write(f"Variants applied:            {total_applied} ({n_snps_total} SNPs, {n_indels_total} indels)\n")
        f.write(f"Variants skipped:            {total_skipped}\n")
        f.write(f"Positions masked (N):        {total_n_count}\n\n")

        # Per-segment details
        for seg_id, sr in segment_results.items():
            f.write("-" * 80 + "\n")
            f.write(f"SEGMENT: {seg_id} ({sr['ref_len']} bp)\n")
            f.write("-" * 80 + "\n")
            n_snps = sum(1 for r in sr['applied']
                         if len(str(r['REF'])) == 1 and len(str(r['ALT'])) == 1)
            n_indels = len(sr['applied']) - n_snps
            f.write(f"  Variants applied: {len(sr['applied'])} ({n_snps} SNPs, {n_indels} indels)\n")
            f.write(f"  Positions masked: {sr['n_count']}\n")

            if sr['applied']:
                f.write("  Mutations:\n")
                for row in sr['applied']:
                    gene = row.get('GENE_NAME', '')
                    hgvsp = row.get('HGVSp', '')
                    effect = row.get('EFFECT', '')
                    af = row['Allele_Frequency']
                    f.write(f"    {row['POS']} {row['REF']}>{row['ALT']}  AF={af:.4f}  "
                            f"gene={gene}  effect={effect}  HGVSp={hgvsp}\n")

            if sr['masked_regions']:
                f.write(f"  Low coverage regions (< {min_depth}X):\n")
                for start, end in sr['masked_regions']:
                    f.write(f"    {start}-{end} ({end - start + 1} bp)\n")

            if sr['proteins']:
                f.write("  Proteins:\n")
                for gene_name, protein in sr['proteins']:
                    f.write(f"    {gene_name}: {len(protein)} aa\n")
            f.write("\n")

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

    # ── Load reference (auto-detect segmented) ────────────────────────
    ref_records = OrderedDict()
    for record in SeqIO.parse(args.reference, "fasta"):
        ref_records[record.id] = str(record.seq)
    is_segmented = len(ref_records) > 1

    if is_segmented:
        total_bp = sum(len(s) for s in ref_records.values())
        print(f"Reference genome: {len(ref_records)} segments, {total_bp} bp total")
        for seg_id, seg_seq in ref_records.items():
            print(f"  {seg_id}: {len(seg_seq)} bp")
    else:
        ref_seq = list(ref_records.values())[0]
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

    # ── Branch: segmented vs single-sequence ──────────────────────────
    if is_segmented:
        _run_segmented(args, ref_records, virus_config, threshold,
                       variants_df, consensus_variants)
    else:
        _run_single(args, ref_seq, virus_config, threshold,
                    variants_df, consensus_variants)

    print("")
    print("=" * 80)
    print("Consensus generation complete")
    print("=" * 80)


def _run_single(args, ref_seq, virus_config, threshold, variants_df, consensus_variants):
    """Original single-sequence consensus pipeline."""
    # ── Build consensus ───────────────────────────────────────────────
    indel_offsets = []
    if len(consensus_variants) == 0:
        print("No variants above consensus threshold — consensus equals reference")
        consensus_seq = ref_seq
        applied = []
        skipped = []
    else:
        consensus_seq, applied, skipped, indel_offsets = apply_variants_to_sequence(
            ref_seq, consensus_variants)
        n_snps = sum(1 for r in applied if len(str(r['REF'])) == 1 and len(str(r['ALT'])) == 1)
        n_indels = len(applied) - n_snps
        print(f"Applied {len(applied)} mutations ({n_snps} SNPs, {n_indels} indels), "
              f"skipped {len(skipped)}")
        if indel_offsets:
            total_offset = sum(s for _, s in indel_offsets)
            print(f"  Net genome size change from indels: {total_offset:+d} bp "
                  f"(ref {len(ref_seq)} → consensus {len(consensus_seq)})")

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

    # Adjust gene coordinates if indels changed the genome length
    adjusted_coords = adjust_gene_coords(gene_coords, indel_offsets) if gene_coords else {}
    if indel_offsets and gene_coords:
        print(f"  Adjusted {len(adjusted_coords)} gene coordinates for indels")

    if adjusted_coords:
        proteins = translate_genes(consensus_seq, adjusted_coords)
        print(f"Translated {len(proteins)} gene products (UTRs excluded)")
    elif gene_coords:
        # No indels — use original coords
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
        # Find mutations in this gene (use original coords for POS matching)
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


def _run_segmented(args, ref_records, virus_config, threshold,
                   variants_df, consensus_variants):
    """Segmented virus consensus pipeline — per-segment processing."""
    sample_name = Path(args.output_prefix).name

    # ── Group variants by CHROM ───────────────────────────────────────
    if 'CHROM' in consensus_variants.columns:
        variants_by_seg = dict(tuple(consensus_variants.groupby('CHROM')))
    else:
        # Fallback: all variants go to the first segment
        first_seg = list(ref_records.keys())[0]
        variants_by_seg = {first_seg: consensus_variants} if len(consensus_variants) > 0 else {}

    # ── Load depth data if available ──────────────────────────────────
    depth_by_seg = None
    if args.depth_file and Path(args.depth_file).exists():
        print(f"Reading depth from: {args.depth_file}")
        depth_by_seg = _load_depth_from_file_segmented(args.depth_file, ref_records)
    elif args.bam and Path(args.bam).exists():
        print(f"Computing depth from BAM: {args.bam}")
        depth_by_seg = _load_depth_from_bam_segmented(args.bam, ref_records)

    # ── Process each segment ──────────────────────────────────────────
    consensus_records = []
    protein_records = []
    segment_results = OrderedDict()

    for seg_id, seg_ref_seq in ref_records.items():
        print(f"\n--- Segment: {seg_id} ({len(seg_ref_seq)} bp) ---")

        seg_variants = variants_by_seg.get(seg_id, pd.DataFrame())
        seg_indel_offsets = []

        if len(seg_variants) == 0:
            print(f"  No variants — consensus equals reference")
            seg_consensus = seg_ref_seq
            seg_applied = []
            seg_skipped = []
        else:
            seg_consensus, seg_applied, seg_skipped, seg_indel_offsets = \
                apply_variants_to_sequence(seg_ref_seq, seg_variants)
            n_snps = sum(1 for r in seg_applied
                         if len(str(r['REF'])) == 1 and len(str(r['ALT'])) == 1)
            n_indels = len(seg_applied) - n_snps
            print(f"  Applied {len(seg_applied)} mutations ({n_snps} SNPs, {n_indels} indels), "
                  f"skipped {len(seg_skipped)}")
            if seg_indel_offsets:
                total_offset = sum(s for _, s in seg_indel_offsets)
                print(f"  Net size change: {total_offset:+d} bp")

        # ── Per-segment depth masking ─────────────────────────────────
        seg_n_count = 0
        seg_masked_regions = []
        if depth_by_seg and seg_id in depth_by_seg:
            print(f"  Masking positions with depth < {args.min_depth}X ...")
            seg_consensus, seg_n_count, seg_masked_regions = \
                mask_low_coverage_from_array(seg_consensus, args.min_depth,
                                             depth_by_seg[seg_id])
            print(f"  Masked {seg_n_count} positions with N ({len(seg_masked_regions)} regions)")
        elif not (args.depth_file or args.bam):
            pass  # No depth source — already reported above

        # ── Build FASTA record for this segment ──────────────────────
        header_parts = [f"segment={seg_id}", f"{len(seg_applied)} mutations"]
        if seg_applied:
            nt_muts = ','.join(f"{int(r['POS'])}{r['REF']}>{r['ALT']}" for r in seg_applied)
            header_parts.append(f"variants={nt_muts}")
        header_parts.append(f"AF>={threshold}")
        if seg_n_count > 0:
            header_parts.append(f"{seg_n_count} Ns (depth<{args.min_depth}X)")
        desc = ' '.join(header_parts)

        rec = SeqRecord(Seq(seg_consensus), id=f"{sample_name}_{seg_id}",
                        description=desc)
        consensus_records.append(rec)

        # ── Per-segment protein translation ───────────────────────────
        seg_gene_coords = get_gene_coords_for_segment(virus_config, seg_id)
        seg_adjusted_coords = adjust_gene_coords(seg_gene_coords, seg_indel_offsets) \
            if seg_indel_offsets and seg_gene_coords else seg_gene_coords

        if seg_adjusted_coords:
            seg_proteins = translate_genes(seg_consensus, seg_adjusted_coords)
        elif seg_gene_coords:
            seg_proteins = translate_genes(seg_consensus, seg_gene_coords)
        else:
            protein = viral_translate(seg_consensus, coordinates=None, stop_at_stop_codon=True)
            seg_proteins = [('Polyprotein', protein)]

        if seg_proteins:
            print(f"  Translated {len(seg_proteins)} proteins: "
                  f"{', '.join(g for g, _ in seg_proteins)}")

        for gene_name, protein_seq in seg_proteins:
            gene_start, gene_end = seg_gene_coords.get(gene_name, (0, 0))
            gene_muts = get_mutations_for_gene(seg_applied, gene_name, gene_start, gene_end)

            desc_parts = [f"{seg_id}:{gene_name} ({len(protein_seq)} aa)"]
            if gene_muts:
                nt_str = format_mutations_for_header(gene_muts, mode='nt')
                aa_str = format_mutations_for_header(gene_muts, mode='aa')
                desc_parts.append(f"nt={nt_str}")
                if aa_str:
                    desc_parts.append(f"aa={aa_str}")
            else:
                desc_parts.append("no mutations")

            prot_rec = SeqRecord(
                Seq(protein_seq),
                id=f"{sample_name}_{seg_id}_{gene_name}",
                description=' '.join(desc_parts)
            )
            protein_records.append(prot_rec)

        segment_results[seg_id] = {
            'ref_len': len(seg_ref_seq),
            'consensus_len': len(seg_consensus),
            'applied': seg_applied,
            'skipped': seg_skipped,
            'indel_offsets': seg_indel_offsets,
            'n_count': seg_n_count,
            'masked_regions': seg_masked_regions,
            'proteins': seg_proteins,
        }

    # ── Write multi-record consensus FASTA ────────────────────────────
    consensus_fasta = f"{args.output_prefix}_consensus.fasta"
    SeqIO.write(consensus_records, consensus_fasta, "fasta")
    print(f"\nConsensus genome written to: {consensus_fasta} "
          f"({len(consensus_records)} segments)")

    # ── Write protein FASTA ───────────────────────────────────────────
    protein_fasta = f"{args.output_prefix}_consensus_proteins.fasta"
    SeqIO.write(protein_records, protein_fasta, "fasta")
    print(f"Consensus proteins written to: {protein_fasta} "
          f"({len(protein_records)} proteins)")

    # ── Write segmented report ────────────────────────────────────────
    report_path = f"{args.output_prefix}_consensus_report.txt"
    write_report_segmented(report_path, args, virus_config, threshold,
                           variants_df, segment_results, args.min_depth)


if __name__ == "__main__":
    main()
