#!/usr/bin/env python3
"""
Generate biologically realistic viral haplotypes based on mutation frequency patterns
Uses frequency-based grouping and realistic population genetics modeling
Correctly handles high-frequency mutations that dominate viral populations
"""

import argparse
import sys; import os; sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "visualization")); from viral_translator import viral_translate
from pathlib import Path
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
from collections import defaultdict
import itertools
from math import prod
import numpy as np

def read_virus_config(accession):
    """Read virus configuration from known_viruses.json"""
    possible_paths = [
        Path(__file__).parent.parent.parent / "virus_configs" / "known_viruses.json",
        Path(__file__).parent / "known_viruses.json",
        Path(__file__).parent.parent / "visualization" / "known_viruses.json",
        Path(__file__).parent.parent.parent / "known_viruses.json",
        Path("known_viruses.json"),
    ]
    
    config_file = None
    for path in possible_paths:
        if path.exists():
            config_file = path
            print(f"Using virus config from: {config_file}")
            break
    
    if not config_file:
        print(f"Warning: known_viruses.json not found in any expected location")
        return None
    
    with open(config_file, 'r') as f:
        viruses = json.load(f)
    
    return viruses.get(accession, None)

def filter_vcf(vcf_file, quality_cutoff, depth_cutoff, freq_cutoff):
    """Filter VCF/TSV file by quality, depth, and frequency"""
    print(f"Reading variants from {vcf_file}")
    
    with open(vcf_file, 'r') as f:
        first_line = f.readline()
        if first_line.startswith('#CHROM'):
            df = pd.read_csv(vcf_file, sep='\t')
        else:
            df = pd.read_csv(vcf_file, sep='\t')
    
    print(f"Found {len(df)} total variants")
    
    filtered = df[
        (df['QUAL'] >= quality_cutoff) & 
        (df['Total_Depth'] >= depth_cutoff) &
        (df['Allele_Frequency'] >= freq_cutoff)
    ].copy()
    
    print(f"After filtering: {len(filtered)} variants pass criteria")
    return filtered

def find_gene_for_position(pos, gene_coords):
    """Find all genes that contain a given position, showing overlaps"""
    matching_genes = []
    for gene, (start, end) in gene_coords.items():
        if start <= pos <= end:
            matching_genes.append((gene, start, end))
    
    if not matching_genes:
        return None, None, None
    
    # Sort by gene length (longest first) for consistent ordering
    matching_genes.sort(key=lambda x: x[2] - x[1] + 1, reverse=True)
    
    # Return combined gene name and coordinates of the longest gene
    gene_names = [g[0] for g in matching_genes]
    combined_gene = "/".join(gene_names)
    longest_start = matching_genes[0][1]
    longest_end = matching_genes[0][2]
    
    return combined_gene, longest_start, longest_end

def parse_aa_change(hgvsp):
    """Parse amino acid change from HGVSp annotation"""
    if not hgvsp or pd.isna(hgvsp):
        return None, None, None
    
    # Handle various formats: p.Leu287Ile, p.L287I, etc
    import re
    pattern = r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)'
    match = re.search(pattern, hgvsp)
    if match:
        ref_aa = match.group(1)
        position = int(match.group(2))
        alt_aa = match.group(3)
        
        # Convert 3-letter to 1-letter if needed
        aa_3to1 = {
            'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
            'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
            'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
            'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
            'Ter': '*', 'TER': '*'
        }
        
        if len(ref_aa) > 1:
            ref_aa = aa_3to1.get(ref_aa, ref_aa)
        if len(alt_aa) > 1:
            alt_aa = aa_3to1.get(alt_aa, alt_aa)
            
        return ref_aa, position, alt_aa
    
    return None, None, None

def reconstruct_realistic_haplotypes(mutations_df, gene_coords, 
                                   high_freq_threshold=0.80, 
                                   medium_freq_threshold=0.40,
                                   min_haplotype_freq=0.01):
    """
    Reconstruct viral haplotypes using biologically realistic frequency-based grouping
    """
    print(f"\n🧬 Reconstructing realistic viral haplotypes...")
    
    # Prepare mutation data
    mutations = []
    for _, mut in mutations_df.iterrows():
        pos = int(mut['POS'])
        gene, start, end = find_gene_for_position(pos, gene_coords)
        
        ref_aa, aa_pos, alt_aa = parse_aa_change(mut.get('HGVSp', ''))
        aa_change = f"{ref_aa}{aa_pos}{alt_aa}" if ref_aa and alt_aa else "Unknown"
        
        mutations.append({
            'pos': pos,
            'gene': gene,
            'freq': mut['Allele_Frequency'],
            'mutation': mut,
            'aa_change': aa_change,
            'nucleotide_id': f"{mut['POS']}{mut['REF']}>{mut['ALT']}",
            'effect': mut.get('EFFECT', 'unknown')
        })
    
    # Sort mutations by frequency (highest first)
    mutations.sort(key=lambda x: x['freq'], reverse=True)
    
    print(f"  📊 Mutation frequency analysis:")
    high_freq_muts = [m for m in mutations if m['freq'] >= high_freq_threshold]
    medium_freq_muts = [m for m in mutations if medium_freq_threshold <= m['freq'] < high_freq_threshold]
    low_freq_muts = [m for m in mutations if m['freq'] < medium_freq_threshold]
    
    print(f"    High frequency (≥{high_freq_threshold*100:.0f}%): {len(high_freq_muts)} mutations")
    print(f"    Medium frequency ({medium_freq_threshold*100:.0f}%-{high_freq_threshold*100:.0f}%): {len(medium_freq_muts)} mutations")
    print(f"    Low frequency (<{medium_freq_threshold*100:.0f}%): {len(low_freq_muts)} mutations")
    
    for mut in mutations:
        freq_category = "HIGH" if mut['freq'] >= high_freq_threshold else \
                       "MED" if mut['freq'] >= medium_freq_threshold else "LOW"
        print(f"    - {mut['gene']}: {mut['aa_change']} ({mut['freq']:.2%}) [{freq_category}]")
    
    # Reconstruct haplotypes using realistic frequency modeling
    haplotypes = []
    
    # Step 1: Create the dominant haplotype with high-frequency mutations
    if high_freq_muts:
        # Calculate realistic frequency for dominant haplotype
        # Use the minimum of high-frequency mutations (most restrictive)
        dominant_freq = min(m['freq'] for m in high_freq_muts)
        
        # But also account for the fact that some high-freq mutations might not co-occur
        # Use a more realistic estimate: average of top frequencies
        top_freqs = [m['freq'] for m in high_freq_muts]
        if len(top_freqs) > 3:
            # For many high-freq mutations, use 95th percentile to be conservative
            dominant_freq = np.percentile(top_freqs, 85)
        else:
            # For few mutations, use minimum
            dominant_freq = min(top_freqs)
        
        haplotypes.append({
            'id': 'Dominant',
            'mutations': high_freq_muts,
            'estimated_freq': dominant_freq,
            'description': f'Dominant strain with {len(high_freq_muts)} high-frequency mutations'
        })
        
        print(f"  🎯 Dominant haplotype: {dominant_freq:.2%} frequency, {len(high_freq_muts)} mutations")
    
    # Step 2: Create medium-frequency variant haplotypes
    if medium_freq_muts:
        # Group medium-frequency mutations by frequency similarity
        freq_groups = []
        current_group = []
        
        for mut in medium_freq_muts:
            if not current_group:
                current_group = [mut]
            else:
                # Group mutations within 10% frequency difference
                freq_diff = abs(current_group[0]['freq'] - mut['freq'])
                if freq_diff <= 0.10:
                    current_group.append(mut)
                else:
                    freq_groups.append(current_group)
                    current_group = [mut]
        
        if current_group:
            freq_groups.append(current_group)
        
        # Create haplotypes for each frequency group
        for i, group in enumerate(freq_groups):
            group_freq = np.mean([m['freq'] for m in group])
            
            # Medium-frequency variants might also include some high-frequency mutations
            # Model as: medium-freq mutations + subset of high-freq mutations
            variant_mutations = group.copy()
            
            # Add high-frequency mutations that are compatible
            for high_mut in high_freq_muts:
                # Compatible if high-freq mutation frequency >= group frequency
                if high_mut['freq'] >= group_freq:
                    variant_mutations.append(high_mut)
            
            # But adjust frequency to be realistic
            # Medium variants can't have higher frequency than their defining mutations
            variant_freq = min(group_freq, group_freq * 0.9)  # Slightly lower due to combination
            
            if variant_freq >= min_haplotype_freq:
                haplotypes.append({
                    'id': f'Medium_variant_{i+1}',
                    'mutations': variant_mutations,
                    'estimated_freq': variant_freq,
                    'description': f'Medium-frequency variant with {len(variant_mutations)} mutations'
                })
                
                print(f"  📊 Medium variant {i+1}: {variant_freq:.2%} frequency, {len(variant_mutations)} mutations")
    
    # Step 3: Create low-frequency variants
    significant_low_freq = [m for m in low_freq_muts if m['freq'] >= min_haplotype_freq * 2]
    
    for mut in significant_low_freq:
        # Low-frequency variants are typically single mutations or small combinations
        variant_mutations = [mut]
        
        # Might co-occur with some high-frequency mutations
        for high_mut in high_freq_muts:
            # Lower probability of co-occurrence
            if high_mut['freq'] >= 0.90:  # Only very high frequency mutations
                variant_mutations.append(high_mut)
        
        variant_freq = mut['freq'] * 0.8  # Adjust for combination uncertainty
        
        if variant_freq >= min_haplotype_freq:
            haplotypes.append({
                'id': f'Low_freq_{mut["aa_change"]}',
                'mutations': variant_mutations,
                'estimated_freq': variant_freq,
                'description': f'Low-frequency variant: {mut["aa_change"]}'
            })
            
            print(f"  📉 Low-freq variant {mut['aa_change']}: {variant_freq:.2%} frequency")
    
    # Step 4: Calculate wild-type frequency
    total_mutant_freq = sum(h['estimated_freq'] for h in haplotypes)
    
    # Realistic wild-type calculation
    if high_freq_muts:
        # If we have high-frequency mutations (>80%), wild-type should be very low
        max_high_freq = max(m['freq'] for m in high_freq_muts)
        wildtype_freq = max(0.01, 1.0 - max_high_freq)  # At least 1%, but based on highest mutation
    else:
        # If no high-frequency mutations, use traditional calculation
        wildtype_freq = max(0.0, 1.0 - total_mutant_freq)
    
    if wildtype_freq >= min_haplotype_freq:
        haplotypes.append({
            'id': 'Wildtype',
            'mutations': [],
            'estimated_freq': wildtype_freq,
            'description': 'Reference sequence (no mutations)'
        })
        print(f"  🧬 Wild-type: {wildtype_freq:.2%} frequency")
    
    # Step 5: Normalize frequencies to sum to 100%
    total_freq = sum(h['estimated_freq'] for h in haplotypes)
    if total_freq > 1.0:
        print(f"  ⚖️  Normalizing frequencies (total was {total_freq:.2%})")
        for haplotype in haplotypes:
            haplotype['estimated_freq'] = haplotype['estimated_freq'] / total_freq
    
    # Sort haplotypes by frequency
    haplotypes.sort(key=lambda x: x['estimated_freq'], reverse=True)
    
    print(f"\n  🎯 Reconstructed {len(haplotypes)} realistic viral haplotypes:")
    for h in haplotypes:
        print(f"    - {h['id']}: {h['estimated_freq']:.2%} frequency, {len(h['mutations'])} mutations")
    
    return haplotypes

def apply_haplotype_mutations(ref_seq, haplotype, haplotype_id=""):
    """Apply mutations from a specific haplotype to reference sequence"""
    seq_list = list(str(ref_seq))
    
    if not haplotype['mutations']:
        # Wild-type haplotype
        print(f"Generated wild-type sequence [{haplotype_id}]")
        return ''.join(seq_list)
    
    # Sort mutations by position
    mutations = sorted(haplotype['mutations'], key=lambda x: x['pos'])
    
    applied_count = 0
    for mut_info in mutations:
        mut = mut_info['mutation']
        pos = int(mut['POS']) - 1
        ref_base = mut['REF']
        alt_base = mut['ALT']
        
        if pos < len(seq_list) and seq_list[pos] == ref_base:
            seq_list[pos] = alt_base
            applied_count += 1
        else:
            # Check if this is a known issue (like position 1793 T/C conflict)
            if seq_list[pos] != ref_base:
                # This might be due to multiple mutations at the same position
                # Apply anyway but warn
                seq_list[pos] = alt_base
                applied_count += 1
                print(f"Applied mutation at position {pos+1} despite reference mismatch: {ref_base}>{alt_base} (found {seq_list[pos]})")
    
    haplotype_suffix = f" [{haplotype_id}]" if haplotype_id else ""
    print(f"Applied {applied_count} mutations to sequence{haplotype_suffix}")
    return ''.join(seq_list)

def generate_haplotype_proteins(haplotypes, ref_seq, virus_config, summary_data):
    """Generate protein sequences for each viral haplotype"""
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    
    if not gene_coords:
        print("No gene coordinates found - generating polyproteins for each haplotype")
        haplotype_proteins = {}
        
        for haplotype in haplotypes:
            haplotype_seq = apply_haplotype_mutations(ref_seq, haplotype, haplotype['id'])
            protein = viral_translate(haplotype_seq, coordinates=None, stop_at_stop_codon=True)
            
            key = (haplotype['id'], 'Polyprotein', protein)
            haplotype_proteins[key] = {
                'sequence': protein,
                'haplotype_info': haplotype,
                'gene': 'Polyprotein',
                'mutations': [m['aa_change'] for m in haplotype['mutations']]
            }
        
        return haplotype_proteins
    
    # Generate proteins for each gene in each haplotype
    haplotype_proteins = {}
    protein_mutation_summary = defaultdict(lambda: defaultdict(list))
    
    for haplotype in haplotypes:
        # Generate genomic sequence for this haplotype
        haplotype_seq = apply_haplotype_mutations(ref_seq, haplotype, haplotype['id'])
        
        # Group mutations by gene for this haplotype
        gene_mutations = defaultdict(list)
        for mut in haplotype['mutations']:
            if mut['gene']:
                gene_mutations[mut['gene']].append(mut)
        
        # Generate proteins for each gene
        for gene, (start, end) in gene_coords.items():
            # Extract gene sequence from haplotype
            gene_seq = haplotype_seq[start-1:end]
            protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
            
            # Get mutations affecting this gene in this haplotype
            gene_muts = gene_mutations.get(gene, [])
            aa_changes = [m['aa_change'] for m in gene_muts if m['aa_change'] != 'Unknown']
            nucleotide_ids = [m['nucleotide_id'] for m in gene_muts]
            
            # Track mutation summary
            for mut in gene_muts:
                if mut['aa_change'] != 'Unknown':
                    protein_mutation_summary[gene][mut['aa_change']].append({
                        'haplotype': haplotype['id'],
                        'nucleotide': mut['nucleotide_id'],
                        'haplotype_freq': haplotype['estimated_freq']
                    })
            
            # Create unique key for this protein
            mutations_tuple = tuple(sorted(aa_changes))
            key = (haplotype['id'], gene, protein, mutations_tuple)
            
            haplotype_proteins[key] = {
                'sequence': protein,
                'haplotype_info': haplotype,
                'gene': gene,
                'mutations': aa_changes,
                'nucleotide_ids': nucleotide_ids
            }
    
    # Update summary data
    summary_data['protein_mutation_summary'] = dict(protein_mutation_summary)
    summary_data['haplotype_analysis'] = haplotypes
    
    return haplotype_proteins

def generate_summary_report(summary_data, output_prefix):
    """Generate comprehensive realistic haplotype reconstruction report"""
    report_file = f"{output_prefix}_realistic_haplotype_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("REALISTIC VIRAL HAPLOTYPE RECONSTRUCTION SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        # Basic statistics
        f.write(f"Total mutations passing filters: {summary_data['total_mutations']}\n")
        f.write(f"Reconstructed haplotypes: {len(summary_data.get('haplotype_analysis', []))}\n")
        f.write(f"Unique protein sequences: {summary_data['unique_proteins']}\n\n")
        
        # Haplotype analysis
        if summary_data.get('haplotype_analysis'):
            f.write("-" * 80 + "\n")
            f.write("REALISTIC VIRAL HAPLOTYPE RECONSTRUCTION\n")
            f.write("-" * 80 + "\n\n")
            
            total_freq = 0
            for haplotype in summary_data['haplotype_analysis']:
                f.write(f"{haplotype['id']}:\n")
                f.write(f"  Estimated frequency: {haplotype['estimated_freq']:.2%}\n")
                f.write(f"  Number of mutations: {len(haplotype['mutations'])}\n")
                f.write(f"  Description: {haplotype['description']}\n")
                
                if haplotype['mutations']:
                    f.write(f"  Mutations:\n")
                    for mut in haplotype['mutations']:
                        f.write(f"    - {mut['gene']}: {mut['aa_change']} ({mut['nucleotide_id']}, {mut['freq']:.2%})\n")
                else:
                    f.write(f"  Wild-type (no mutations)\n")
                
                total_freq += haplotype['estimated_freq']
                f.write("\n")
            
            f.write(f"Total frequency accounted for: {total_freq:.2%}\n\n")
        
        # Protein summary by haplotype
        if summary_data['protein_mutation_summary']:
            f.write("-" * 80 + "\n")
            f.write("PROTEIN VARIANTS BY HAPLOTYPE\n")
            f.write("-" * 80 + "\n\n")
            
            for gene, mutations in summary_data['protein_mutation_summary'].items():
                f.write(f"{gene}:\n")
                if mutations:
                    for aa_change, occurrences in mutations.items():
                        f.write(f"  - {aa_change}:\n")
                        for occ in occurrences:
                            f.write(f"    {occ['haplotype']}: {occ['nucleotide']} (haplotype freq: {occ['haplotype_freq']:.2%})\n")
                else:
                    f.write("  No mutations across haplotypes\n")
                f.write("\n")
    
    print(f"\n📄 Realistic haplotype reconstruction report saved to: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate realistic viral haplotypes based on mutation frequencies')
    parser.add_argument('--vcf', required=True, help='Input VCF or filtered TSV file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--accession', required=True, help='Virus accession number')
    parser.add_argument('--quality', type=float, default=1000, help='Minimum quality score (default: 1000)')
    parser.add_argument('--depth', type=int, default=200, help='Minimum depth (default: 200)')
    parser.add_argument('--freq', type=float, default=0.01, help='Minimum allele frequency (default: 0.01)')
    parser.add_argument('--high-freq-threshold', type=float, default=0.80, help='High frequency threshold (default: 0.80)')
    parser.add_argument('--medium-freq-threshold', type=float, default=0.40, help='Medium frequency threshold (default: 0.40)')
    parser.add_argument('--min-haplotype-freq', type=float, default=0.01, help='Minimum haplotype frequency (default: 0.01)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("REALISTIC VIRAL HAPLOTYPE RECONSTRUCTION")
    print("=" * 80)
    print(f"Input VCF/TSV: {args.vcf}")
    print(f"Reference: {args.reference}")
    print(f"Accession: {args.accession}")
    print(f"Quality filter: >= {args.quality}")
    print(f"Depth filter: >= {args.depth}")
    print(f"Frequency filter: >= {args.freq}")
    print(f"High frequency threshold: >= {args.high_freq_threshold}")
    print(f"Medium frequency threshold: >= {args.medium_freq_threshold}")
    print(f"Min haplotype frequency: >= {args.min_haplotype_freq}")
    print(f"Output prefix: {args.output_prefix}")
    print("=" * 80)
    
    # Read reference sequence
    ref_record = next(SeqIO.parse(args.reference, "fasta"))
    ref_seq = str(ref_record.seq)
    print(f"Reference genome length: {len(ref_seq)} bp")
    
    # Read virus configuration
    virus_config = read_virus_config(args.accession)
    if not virus_config:
        print(f"Warning: No configuration found for {args.accession}, will generate polyproteins only")
    else:
        print(f"Found configuration for {virus_config['name']}")
        gene_count = len(virus_config.get('gene_coords', {}))
        print(f"Genes defined: {gene_count}")
    
    # Filter mutations
    mutations_df = filter_vcf(args.vcf, args.quality, args.depth, args.freq)
    
    if len(mutations_df) == 0:
        print("\n🧬 No mutations pass the filtering criteria")
        print("📝 Generating wild-type only analysis...")
        
        # Initialize summary data for wild-type only
        summary_data = {
            'total_mutations': 0,
            'unique_proteins': 0,
            'protein_mutation_summary': {},
            'haplotype_analysis': [{
                'id': 'Wildtype',
                'mutations': [],
                'estimated_freq': 1.0,
                'description': 'Reference sequence - 100% wild-type population'
            }]
        }
        
        # Create wild-type haplotype
        wildtype_haplotype = {
            'id': 'Wildtype',
            'mutations': [],
            'estimated_freq': 1.0,
            'description': 'Reference sequence - 100% wild-type population'
        }
        
        print(f"\n🎯 Analysis Results:")
        print(f"  Wild-type: 100.00% frequency, 0 mutations")
        print(f"  Population: Homogeneous wild-type virus")
        
        # Generate wild-type proteins
        print(f"\n🧬 Generating wild-type proteins...")
        gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
        
        if gene_coords:
            protein_records = []
            sample_name = Path(args.output_prefix).name
            
            # Generate proteins for each gene using reference sequence
            for gene, (start, end) in gene_coords.items():
                gene_seq = ref_seq[start-1:end]
                protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
                
                desc = f"{gene} protein (Wild-type, 100.00%)"
                
                record = SeqRecord(
                    Seq(protein),
                    id=f"{sample_name}_Wildtype_{gene}",
                    description=desc
                )
                protein_records.append(record)
            
            summary_data['unique_proteins'] = len(protein_records)
            
            # Write wild-type proteins
            protein_file = f"{args.output_prefix}_realistic_proteins.fasta"
            SeqIO.write(protein_records, protein_file, "fasta")
            print(f"Generated {len(protein_records)} wild-type protein sequences")
            print(f"Wrote proteins to: {protein_file}")
        else:
            # Generate single polyprotein
            protein = viral_translate(ref_seq, coordinates=None, stop_at_stop_codon=True)
            
            protein_record = SeqRecord(
                Seq(protein),
                id=f"{sample_name}_Wildtype_Polyprotein",
                description="Polyprotein (Wild-type, 100.00%)"
            )
            
            protein_file = f"{args.output_prefix}_realistic_proteins.fasta"
            SeqIO.write([protein_record], protein_file, "fasta")
            print(f"Generated 1 wild-type polyprotein sequence")
            print(f"Wrote proteins to: {protein_file}")
            
            summary_data['unique_proteins'] = 1
        
        # Write wild-type consensus sequence
        print(f"\n📝 Writing wild-type consensus sequence...")
        sample_name = Path(args.output_prefix).name
        
        consensus_record = SeqRecord(
            Seq(ref_seq),
            id=f"{sample_name}_Wildtype",
            description="Wild-type reference sequence - 100% frequency"
        )
        
        consensus_file = f"{args.output_prefix}_realistic_haplotypes.fasta"
        SeqIO.write([consensus_record], consensus_file, "fasta")
        print(f"Wrote wild-type sequence to: {consensus_file}")
        
        # Generate wild-type summary report
        report_file = f"{args.output_prefix}_realistic_haplotype_report.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("REALISTIC VIRAL HAPLOTYPE RECONSTRUCTION SUMMARY REPORT\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("ANALYSIS RESULT: 100% WILD-TYPE POPULATION\n\n")
            
            f.write(f"Total mutations passing filters: 0\n")
            f.write(f"Quality filter: >= {args.quality}\n")
            f.write(f"Depth filter: >= {args.depth}\n")
            f.write(f"Frequency filter: >= {args.freq}\n\n")
            
            f.write(f"Reconstructed haplotypes: 1\n")
            f.write(f"Unique protein sequences: {summary_data['unique_proteins']}\n\n")
            
            f.write("-" * 80 + "\n")
            f.write("VIRAL POPULATION STRUCTURE\n")
            f.write("-" * 80 + "\n\n")
            
            f.write("Wildtype:\n")
            f.write("  Estimated frequency: 100.00%\n")
            f.write("  Number of mutations: 0\n")
            f.write("  Description: Reference sequence - homogeneous wild-type population\n")
            f.write("  Interpretation: No significant mutations detected above threshold\n\n")
            
            f.write("Total frequency accounted for: 100.00%\n\n")
            
            f.write("-" * 80 + "\n")
            f.write("BIOLOGICAL INTERPRETATION\n")
            f.write("-" * 80 + "\n\n")
            
            f.write("Population Status: Homogeneous wild-type virus\n")
            f.write("Genetic Diversity: Minimal - no mutations above detection threshold\n")
            f.write("Clinical Significance: Standard reference strain characteristics\n")
            f.write("Evolutionary Status: Stable, no detectable adaptation\n\n")
        
        print(f"📄 Wild-type analysis report saved to: {report_file}")
        
        print(f"\n✅ Wild-type analysis complete!")
        print(f"Result: 100% wild-type viral population with no significant mutations")
        
        return
    
    # Initialize summary data
    summary_data = {
        'total_mutations': len(mutations_df),
        'unique_proteins': 0,
        'protein_mutation_summary': {},
        'haplotype_analysis': []
    }
    
    # Reconstruct realistic viral haplotypes
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    haplotypes = reconstruct_realistic_haplotypes(
        mutations_df, gene_coords, 
        args.high_freq_threshold, 
        args.medium_freq_threshold, 
        args.min_haplotype_freq
    )
    
    # Generate proteins for each haplotype
    print(f"\n🧬 Generating haplotype-specific proteins...")
    haplotype_proteins = generate_haplotype_proteins(haplotypes, ref_seq, virus_config, summary_data)
    
    print(f"Generated {len(haplotype_proteins)} haplotype-specific protein sequences")
    summary_data['unique_proteins'] = len(haplotype_proteins)
    
    # Write consensus sequences for each haplotype
    print(f"\n📝 Writing haplotype consensus sequences...")
    sample_name = Path(args.output_prefix).name
    
    consensus_records = []
    for haplotype in haplotypes:
        haplotype_seq = apply_haplotype_mutations(ref_seq, haplotype, haplotype['id'])
        
        # Create description
        if haplotype['mutations']:
            mut_summary = ", ".join([f"{m['aa_change']}" for m in haplotype['mutations'][:3]])
            if len(haplotype['mutations']) > 3:
                mut_summary += f" (+{len(haplotype['mutations'])-3} more)"
            desc = f"{haplotype['id']} with {len(haplotype['mutations'])} mutations [{mut_summary}] (est. {haplotype['estimated_freq']:.2%})"
        else:
            desc = f"{haplotype['id']} - Wild-type reference (est. {haplotype['estimated_freq']:.2%})"
        
        record = SeqRecord(
            Seq(haplotype_seq),
            id=f"{sample_name}_{haplotype['id']}",
            description=desc
        )
        consensus_records.append(record)
    
    consensus_file = f"{args.output_prefix}_realistic_haplotypes.fasta"
    SeqIO.write(consensus_records, consensus_file, "fasta")
    print(f"Wrote {len(consensus_records)} haplotype sequences to: {consensus_file}")
    
    # Write protein sequences
    print(f"\n🧬 Writing haplotype-specific proteins...")
    protein_records = []
    
    # Group proteins by haplotype for organized output
    haplotype_groups = defaultdict(list)
    for key, protein_data in haplotype_proteins.items():
        if len(key) == 3:
            haplotype_id, gene, protein = key
            mutations = tuple()
        elif len(key) == 4:
            haplotype_id, gene, protein, mutations = key
        else:
            print(f"Warning: Unexpected key format: {key}")
            continue
        
        haplotype_groups[haplotype_id].append((key, protein_data))
    
    # Sort haplotypes by frequency
    sorted_haplotypes = sorted(haplotype_groups.keys(), 
                              key=lambda h: next(hp['haplotype_info']['estimated_freq'] 
                                                 for hp in [v[1] for v in haplotype_groups[h]] if hp), 
                              reverse=True)
    
    for haplotype_id in sorted_haplotypes:
        proteins = haplotype_groups[haplotype_id]
        
        for key, protein_data in proteins:
            gene = protein_data['gene']
            haplotype_info = protein_data['haplotype_info']
            
            # Create description
            if protein_data['mutations']:
                mutation_str = ", ".join(protein_data['mutations'])
                desc = f"{gene} protein [{mutation_str}]"
            else:
                desc = f"{gene} protein"
            
            # Add haplotype and frequency info
            desc += f" ({haplotype_info['id']}, est. {haplotype_info['estimated_freq']:.2%})"
            
            # Add nucleotide changes if available
            if protein_data.get('nucleotide_ids'):
                nuc_changes = ", ".join(protein_data['nucleotide_ids'])
                desc += f" (from {nuc_changes})"
            
            record = SeqRecord(
                Seq(protein_data['sequence']),
                id=f"{sample_name}_{haplotype_id}_{gene}",
                description=desc
            )
            protein_records.append(record)
    
    protein_file = f"{args.output_prefix}_realistic_proteins.fasta"
    SeqIO.write(protein_records, protein_file, "fasta")
    print(f"Wrote {len(protein_records)} haplotype-specific proteins to: {protein_file}")
    
    # Generate summary report
    generate_summary_report(summary_data, args.output_prefix)
    
    print(f"\n✅ Realistic viral haplotype reconstruction complete!")
    print(f"Reconstructed {len(haplotypes)} biologically realistic viral haplotypes")
    print(f"Generated {len(haplotype_proteins)} haplotype-specific protein sequences")
    
    # Print haplotype summary
    print(f"\n🎯 Realistic Haplotype Summary:")
    for haplotype in haplotypes:
        print(f"  {haplotype['id']}: {haplotype['estimated_freq']:.2%} frequency, {len(haplotype['mutations'])} mutations")

if __name__ == "__main__":
    main()