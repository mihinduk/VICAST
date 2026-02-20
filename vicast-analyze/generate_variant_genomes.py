#!/usr/bin/env python3
"""
Generate BAM-validated minor variant genomes and per-gene proteins.

Tier 2 of the two-tier VICAST post-processing pipeline:
  - Selects minor variants within a user-specified AF range (default 0.5%-5%)
  - Uses BAM co-occurrence data to identify linked variant groups
  - Builds one variant genome per linked group (applied on top of consensus)
  - Unlinked variants become singleton variant genomes
  - Translates per-gene proteins for each variant genome
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


def read_virus_config(accession):
    """Read virus configuration from known_viruses.json."""
    search_paths = [
        Path(__file__).parent / "known_viruses.json",
        Path(__file__).parent.parent / "virus_configs" / "known_viruses.json",
        Path("known_viruses.json"),
    ]
    for path in search_paths:
        if path.exists():
            with open(path, 'r') as f:
                viruses = json.load(f)
            return viruses.get(accession)
    return None


def load_filtered_tsv(tsv_path):
    """Load the filtered mutations TSV."""
    return pd.read_csv(tsv_path, sep='\t')


def load_cooccurrence(cooccurrence_path):
    """Load co-occurrence TSV from check_read_cooccurrence.py."""
    if not Path(cooccurrence_path).exists():
        print(f"Warning: co-occurrence file not found: {cooccurrence_path}")
        return pd.DataFrame()
    df = pd.read_csv(cooccurrence_path, sep='\t')
    if len(df) == 0:
        print("Co-occurrence file is empty — no linked pairs found")
    return df


def select_minor_variants(variants_df, min_af, max_af):
    """Select variants within the minor AF range."""
    minor = variants_df[
        (variants_df['Allele_Frequency'] >= min_af) &
        (variants_df['Allele_Frequency'] <= max_af)
    ].copy()
    return minor


def build_linkage_graph(minor_variants, cooccurrence_df):
    """Build an adjacency list from co-occurrence data.

    A variant pair is linked if both_alt > 0 in the co-occurrence table.
    Variants are identified by their 1-based POS.

    Returns dict mapping POS -> set of linked POS values.
    """
    graph = defaultdict(set)
    minor_positions = set(minor_variants['POS'].astype(int))

    if cooccurrence_df.empty:
        return graph

    for _, row in cooccurrence_df.iterrows():
        pos1 = int(row.get('variant1_pos', 0))
        pos2 = int(row.get('variant2_pos', 0))

        # Parse both_alt — might be int or string
        both_alt = row.get('both_alt', 0)
        try:
            both_alt = int(both_alt)
        except (ValueError, TypeError):
            both_alt = 0

        if both_alt > 0 and pos1 in minor_positions and pos2 in minor_positions:
            graph[pos1].add(pos2)
            graph[pos2].add(pos1)

    return graph


def find_connected_components(graph, all_positions):
    """Find connected components in the linkage graph.

    Returns list of sets, each set being a group of linked positions.
    Positions not in the graph become singleton groups.
    """
    visited = set()
    components = []

    def bfs(start):
        component = set()
        queue = [start]
        while queue:
            node = queue.pop(0)
            if node in visited:
                continue
            visited.add(node)
            component.add(node)
            for neighbour in graph.get(node, set()):
                if neighbour not in visited:
                    queue.append(neighbour)
        return component

    # Process nodes in the graph first
    for node in sorted(graph.keys()):
        if node not in visited:
            comp = bfs(node)
            components.append(comp)

    # Add singletons for positions not in any graph edge
    for pos in sorted(all_positions):
        if pos not in visited:
            components.append({pos})

    return components


def apply_mutations(base_seq, mutations_df):
    """Apply SNP mutations to a sequence.

    Returns (mutated_seq, applied_count).
    """
    seq_list = list(base_seq)
    applied = 0

    for _, row in mutations_df.sort_values('POS').iterrows():
        pos = int(row['POS']) - 1  # 0-based
        ref_base = str(row['REF'])
        alt_base = str(row['ALT'])

        if len(ref_base) == 1 and len(alt_base) == 1 and pos < len(seq_list):
            seq_list[pos] = alt_base
            applied += 1

    return ''.join(seq_list), applied


def translate_genes(seq, gene_coords):
    """Translate per-gene CDS from a genome sequence."""
    proteins = []
    for gene, (start, end) in gene_coords.items():
        gene_seq = seq[start - 1:end]
        protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
        proteins.append((gene, protein))
    return proteins


def write_report(report_path, args, groups, minor_variants, cooccurrence_df,
                 virus_config, variant_genomes_info):
    """Write a human-readable variant genomes report."""
    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("VICAST VARIANT GENOMES REPORT (Tier 2)\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Accession:            {args.accession}\n")
        if virus_config:
            f.write(f"Virus name:           {virus_config.get('name', 'N/A')}\n")
        f.write(f"Minor AF range:       {args.min_af}-{args.max_af}\n")
        f.write(f"Consensus FASTA:      {args.consensus}\n")
        f.write(f"Co-occurrence file:   {args.cooccurrence}\n\n")

        f.write(f"Minor variants found:           {len(minor_variants)}\n")
        n_pairs = len(cooccurrence_df) if not cooccurrence_df.empty else 0
        f.write(f"Co-occurrence pairs tested:     {n_pairs}\n")
        linked_pairs = 0
        if not cooccurrence_df.empty and 'both_alt' in cooccurrence_df.columns:
            linked_pairs = int((cooccurrence_df['both_alt'].astype(int) > 0).sum())
        f.write(f"Linked pairs (both_alt > 0):    {linked_pairs}\n")
        f.write(f"Variant groups (genomes):       {len(groups)}\n\n")

        f.write("-" * 80 + "\n")
        f.write("VARIANT GROUPS\n")
        f.write("-" * 80 + "\n\n")

        for info in variant_genomes_info:
            f.write(f"Group: {info['name']}\n")
            f.write(f"  Positions:  {', '.join(str(p) for p in sorted(info['positions']))}\n")
            f.write(f"  Mutations:  {info['applied_count']} applied\n")
            f.write(f"  Linkage:    {info['linkage_type']}\n")
            for _, row in info['mutations_df'].iterrows():
                gene = row.get('GENE_NAME', '')
                hgvsp = row.get('HGVSp', '')
                effect = row.get('EFFECT', '')
                af = row['Allele_Frequency']
                f.write(f"    {row['POS']} {row['REF']}>{row['ALT']}  "
                        f"AF={af:.4f}  gene={gene}  effect={effect}  HGVSp={hgvsp}\n")
            f.write("\n")

    print(f"Report written to: {report_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate BAM-linked minor variant genomes and proteins')
    parser.add_argument('--vcf', required=True,
                        help='Filtered TSV from parse_snpeff_tsv.py')
    parser.add_argument('--cooccurrence', required=True,
                        help='Co-occurrence TSV from check_read_cooccurrence.py')
    parser.add_argument('--consensus', required=True,
                        help='Consensus FASTA from generate_consensus_genome.py')
    parser.add_argument('--accession', required=True,
                        help='Virus accession')
    parser.add_argument('--min-af', type=float, default=0.005,
                        help='Minimum AF for minor variants (default: 0.005 = 0.5%%)')
    parser.add_argument('--max-af', type=float, default=0.05,
                        help='Maximum AF for minor variants (default: 0.05 = 5%%)')
    parser.add_argument('--output-prefix', required=True,
                        help='Output prefix for variant genome files')
    args = parser.parse_args()

    print("=" * 80)
    print("VICAST VARIANT GENOME GENERATION (Tier 2)")
    print("=" * 80)

    # ── Load consensus genome ─────────────────────────────────────────
    consensus_record = next(SeqIO.parse(args.consensus, "fasta"))
    consensus_seq = str(consensus_record.seq)
    print(f"Consensus genome: {len(consensus_seq)} bp")

    # ── Load virus config ─────────────────────────────────────────────
    virus_config = read_virus_config(args.accession)
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    if virus_config:
        print(f"Virus config: {virus_config.get('name', 'N/A')} ({len(gene_coords)} genes)")

    # ── Load variants ─────────────────────────────────────────────────
    variants_df = load_filtered_tsv(args.vcf)
    print(f"Total filtered variants: {len(variants_df)}")

    minor_variants = select_minor_variants(variants_df, args.min_af, args.max_af)
    print(f"Minor variants (AF {args.min_af}-{args.max_af}): {len(minor_variants)}")

    if len(minor_variants) == 0:
        print("\nNo minor variants in the specified AF range.")
        print("No variant genomes to generate.")

        # Write empty outputs
        genome_fasta = f"{args.output_prefix}_variant_genomes.fasta"
        protein_fasta = f"{args.output_prefix}_variant_proteins.fasta"
        report_path = f"{args.output_prefix}_variant_report.txt"

        with open(genome_fasta, 'w') as f:
            pass
        with open(protein_fasta, 'w') as f:
            pass
        with open(report_path, 'w') as f:
            f.write("No minor variants found in the specified AF range.\n")

        print(f"Empty outputs written to {args.output_prefix}_variant_*")
        return

    # ── Load co-occurrence ────────────────────────────────────────────
    cooccurrence_df = load_cooccurrence(args.cooccurrence)

    # ── Build linkage graph & find groups ─────────────────────────────
    graph = build_linkage_graph(minor_variants, cooccurrence_df)
    all_positions = set(minor_variants['POS'].astype(int))
    groups = find_connected_components(graph, all_positions)
    print(f"Variant groups: {len(groups)} (from {len(all_positions)} positions)")

    # ── Generate variant genomes ──────────────────────────────────────
    sample_name = Path(args.output_prefix).name
    genome_records = []
    protein_records = []
    variant_genomes_info = []

    for i, group_positions in enumerate(sorted(groups, key=lambda g: min(g))):
        group_df = minor_variants[minor_variants['POS'].isin(group_positions)].copy()
        n_muts = len(group_df)

        # Name the group
        if n_muts == 1:
            row = group_df.iloc[0]
            group_name = f"var_{row['POS']}_{row['REF']}_{row['ALT']}"
        else:
            group_name = f"linked_group_{i+1}"

        # Determine linkage type
        if n_muts > 1 and any(pos in graph for pos in group_positions):
            linkage_type = "BAM-validated co-occurrence"
        elif n_muts > 1:
            linkage_type = "unvalidated — beyond insert size"
        else:
            linkage_type = "singleton"

        # Apply mutations on top of consensus
        variant_seq, applied_count = apply_mutations(consensus_seq, group_df)

        # Create genome record
        positions_str = ','.join(str(p) for p in sorted(group_positions))
        desc = (f"variant genome ({applied_count} mutations at pos {positions_str}; "
                f"linkage: {linkage_type})")
        genome_rec = SeqRecord(
            Seq(variant_seq),
            id=f"{sample_name}_{group_name}",
            description=desc
        )
        genome_records.append(genome_rec)

        # Translate proteins
        if gene_coords:
            proteins = translate_genes(variant_seq, gene_coords)
        else:
            protein = viral_translate(variant_seq, coordinates=None, stop_at_stop_codon=True)
            proteins = [('Polyprotein', protein)]

        for gene_name, protein_seq in proteins:
            prot_rec = SeqRecord(
                Seq(protein_seq),
                id=f"{sample_name}_{group_name}_{gene_name}",
                description=f"{gene_name} from {group_name} ({len(protein_seq)} aa)"
            )
            protein_records.append(prot_rec)

        variant_genomes_info.append({
            'name': group_name,
            'positions': group_positions,
            'applied_count': applied_count,
            'linkage_type': linkage_type,
            'mutations_df': group_df,
        })

        print(f"  {group_name}: {applied_count} mutations, linkage={linkage_type}")

    # ── Write outputs ─────────────────────────────────────────────────
    genome_fasta = f"{args.output_prefix}_variant_genomes.fasta"
    SeqIO.write(genome_records, genome_fasta, "fasta")
    print(f"Variant genomes written to: {genome_fasta} ({len(genome_records)} genomes)")

    protein_fasta = f"{args.output_prefix}_variant_proteins.fasta"
    SeqIO.write(protein_records, protein_fasta, "fasta")
    print(f"Variant proteins written to: {protein_fasta} ({len(protein_records)} proteins)")

    report_path = f"{args.output_prefix}_variant_report.txt"
    write_report(report_path, args, groups, minor_variants, cooccurrence_df,
                 virus_config, variant_genomes_info)

    print("")
    print("=" * 80)
    print(f"Variant genome generation complete: {len(genome_records)} genomes")
    print("=" * 80)


if __name__ == "__main__":
    main()
