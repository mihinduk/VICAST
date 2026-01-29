#!/usr/bin/env python3
"""
Generate consensus genome and proteins with intelligent multi-allelic site handling
Creates separate sequences for each unique mutation to support quasispecies analysis
Generates comprehensive summary report
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
    """Find which gene contains a given position"""
    for gene, (start, end) in gene_coords.items():
        # Skip UTRs - not translated
        if "UTR" in gene:
            continue
        if start <= pos <= end:
            return gene, start, end
    return None, None, None

def parse_aa_change(hgvsp):
    """Parse amino acid change from HGVSp annotation"""
    if not hgvsp or pd.isna(hgvsp):
        return None, None, None
    if not hgvsp or pd.isna(hgvsp):
        return None, None, None
    
    import re
    
    # Handle insertions: p.Arg214_Asp215insGluProGlu
    ins_pattern = r'p\.([A-Za-z]+)(\d+)_([A-Za-z]+)(\d+)ins([A-Za-z]+)'
    ins_match = re.search(ins_pattern, hgvsp)
    if ins_match:
        # Return as "R214_D215insEPE" format
        ref_aa1 = convert_aa_code(ins_match.group(1))
        pos1 = int(ins_match.group(2))
        ref_aa2 = convert_aa_code(ins_match.group(3))
        pos2 = int(ins_match.group(4))
        inserted = convert_aa_code(ins_match.group(5))
        return f"{ref_aa1}{pos1}_{ref_aa2}{pos2}", None, f"ins{inserted}"
    
    # Handle deletions: p.Arg214_Asp215del or p.Arg214del
    del_range_pattern = r'p\.([A-Za-z]+)(\d+)_([A-Za-z]+)(\d+)del'
    del_match = re.search(del_range_pattern, hgvsp)
    if del_match:
        ref_aa1 = convert_aa_code(del_match.group(1))
        pos1 = int(del_match.group(2))
        ref_aa2 = convert_aa_code(del_match.group(3))
        pos2 = int(del_match.group(4))
        return f"{ref_aa1}{pos1}_{ref_aa2}{pos2}", None, "del"
    
    del_single_pattern = r'p\.([A-Za-z]+)(\d+)del'
    del_single_match = re.search(del_single_pattern, hgvsp)
    if del_single_match:
        ref_aa = convert_aa_code(del_single_match.group(1))
        pos = int(del_single_match.group(2))
        return ref_aa, pos, "del"
    
    # Handle delins: p.Arg214_Asp215delinsXxx
    delins_pattern = r'p\.([A-Za-z]+)(\d+)_([A-Za-z]+)(\d+)delins([A-Za-z]+)'
    delins_match = re.search(delins_pattern, hgvsp)
    if delins_match:
        ref_aa1 = convert_aa_code(delins_match.group(1))
        pos1 = int(delins_match.group(2))
        ref_aa2 = convert_aa_code(delins_match.group(3))
        pos2 = int(delins_match.group(4))
        inserted = convert_aa_code(delins_match.group(5))
        return f"{ref_aa1}{pos1}_{ref_aa2}{pos2}", None, f"delins{inserted}"
    
    # Handle simple substitutions: p.Leu287Ile
    sub_pattern = r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)'
    sub_match = re.search(sub_pattern, hgvsp)
    if sub_match:
        ref_aa = convert_aa_code(sub_match.group(1))
        position = int(sub_match.group(2))
        alt_aa = convert_aa_code(sub_match.group(3))
        return ref_aa, position, alt_aa
    
    return None, None, None

def convert_aa_code(aa_code):
    """Convert 3-letter or 1-letter amino acid code to 1-letter"""
    if len(aa_code) == 1:
        return aa_code
    
    aa_3to1 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Ter': '*', 'TER': '*'
    }
    
    # Handle multi-AA codes (e.g., "GluProGlu" -> "EPE")
    result = ""
    i = 0
    while i < len(aa_code):
        # Try 3-letter code
        if i + 3 <= len(aa_code):
            three_letter = aa_code[i:i+3]
            if three_letter in aa_3to1:
                result += aa_3to1[three_letter]
                i += 3
                continue
        # If not found, keep as-is and move on
        result += aa_code[i]
        i += 1
    
    return result if result else aa_code

def apply_mutations_single_sequence(ref_seq, mutations_df, allele_id=""):
    """Apply mutations to single reference sequence"""
    seq_list = list(str(ref_seq))
    mutations_df = mutations_df.sort_values('POS')
    
    applied_count = 0
    for _, mut in mutations_df.iterrows():
        pos = int(mut['POS']) - 1
        ref_base = mut['REF']
        alt_base = mut['ALT']
        
        if pos < len(seq_list) and seq_list[pos] == ref_base:
            seq_list[pos] = alt_base
            applied_count += 1
        else:
            print(f"Warning: Reference mismatch at position {pos+1}: expected {ref_base}, found {seq_list[pos] if pos < len(seq_list) else 'out of range'}")
    
    allele_suffix = f" [{allele_id}]" if allele_id else ""
    print(f"Applied {applied_count} mutations to sequence{allele_suffix}")
    return ''.join(seq_list)

def generate_individual_variant_proteins(ref_seq, mutations_df, virus_config, summary_data):
    """Generate separate protein sequences for each individual mutation for quasispecies analysis"""
    # Get gene coordinates
    gene_coords = virus_config.get('gene_coords', {}) if virus_config else {}
    
    if not gene_coords:
        print("No gene coordinates found - generating single polyprotein consensus")
        consensus_seq = apply_mutations_single_sequence(ref_seq, mutations_df)
        protein = viral_translate(consensus_seq, coordinates=None, stop_at_stop_codon=False)
        return {('Polyprotein', protein, tuple()): {
            'sequence': protein,
            'allele_ids': [''],
            'mutations': [],
            'gene': 'Polyprotein'
        }}
    
    all_proteins = {}
    protein_mutation_summary = defaultdict(lambda: defaultdict(list))
    
    # Generate consensus with all mutations first (for genes without mutations)
    full_consensus = apply_mutations_single_sequence(ref_seq, mutations_df, "full_consensus")
    
    # Group mutations by gene
    gene_mutations = defaultdict(list)
    for _, mut in mutations_df.iterrows():
        pos = int(mut['POS'])
        gene, start, end = find_gene_for_position(pos, gene_coords)
        if gene:
            gene_mutations[gene].append(mut)
    
    # For each gene, create proteins
    for gene, (start, end) in gene_coords.items():
        # Skip UTRs - they are not translated
        if "UTR" in gene:
            continue
        has_nonsynonymous = False
        if gene in gene_mutations:
            # Check if any mutations are non-synonymous
            for mut in gene_mutations[gene]:
                ref_aa_check, aa_pos_check, alt_aa_check = parse_aa_change(mut.get("HGVSp", ""))
                if ref_aa_check and alt_aa_check:
                    # Indels (ins, del, delins) are always non-synonymous
                    if alt_aa_check.startswith('ins') or alt_aa_check.startswith('del'):
                        has_nonsynonymous = True
                        break
                    # Substitutions: check if different
                    if ref_aa_check != alt_aa_check:
                        has_nonsynonymous = True
                        break

        # Add ALL mutations to report (including synonymous)
        if gene in gene_mutations:
            for mut in gene_mutations[gene]:
                ref_aa_all, aa_pos_all, alt_aa_all = parse_aa_change(mut.get("HGVSp", ""))
                if ref_aa_all and alt_aa_all:
                    # Format AA change appropriately
                    if aa_pos_all is not None:
                        aa_change_all = f"{ref_aa_all}{aa_pos_all}{alt_aa_all}"
                    else:
                        # For insertions/deletions, ref_aa already contains position info
                        aa_change_all = f"{ref_aa_all}{alt_aa_all}"
                    mutation_id_all = f"{mut['POS']}{mut['REF']}>{mut['ALT']}"
                    protein_mutation_summary[gene][aa_change_all].append({
                        'nucleotide': mutation_id_all,
                        'effect': mut.get('EFFECT', 'unknown'),
                        'freq': mut.get('Allele_Frequency', 0)
                    })
        
        if gene in gene_mutations and has_nonsynonymous:
            # This gene has mutations - create separate proteins for each mutation
            gene_muts = gene_mutations[gene]

            # Apply ALL consensus mutations to create ONE protein
            # CRITICAL: Extract gene from reference FIRST, then apply mutations
            # This avoids coordinate shift issues from upstream indels
            ref_gene_seq = ref_seq[start-1:end]
            
            # Filter mutations to only those within this gene's coordinates
            gene_specific_muts = [mut for mut in gene_muts if start <= int(mut['POS']) <= end]
            
            if gene_specific_muts:
                # Adjust mutation positions to be relative to gene start (0-indexed)
                gene_muts_df = pd.DataFrame(gene_specific_muts)
                gene_muts_df['POS_ADJUSTED'] = gene_muts_df['POS'].astype(int) - (start - 1)
                
                # Apply mutations to gene sequence
                gene_seq_list = list(ref_gene_seq)
                for _, mut in gene_muts_df.iterrows():
                    pos_in_gene = int(mut['POS_ADJUSTED']) - 1  # Convert to 0-indexed
                    ref_base = mut['REF']
                    alt_base = mut['ALT']
                    
                    # Handle SNPs, insertions, and deletions
                    ref_len = len(ref_base)
                    alt_len = len(alt_base)
                    
                    # Check if position is valid and ref matches
                    if pos_in_gene < len(gene_seq_list):
                        # For indels, we need to handle differently
                        if ref_len == 1 and alt_len == 1:
                            # Simple SNP
                            if gene_seq_list[pos_in_gene] == ref_base:
                                gene_seq_list[pos_in_gene] = alt_base
                        elif ref_len < alt_len:
                            # Insertion: replace ref and insert extra bases
                            if pos_in_gene + ref_len <= len(gene_seq_list):
                                # Replace the ref bases
                                gene_seq_list[pos_in_gene:pos_in_gene+ref_len] = list(alt_base)
                        elif ref_len > alt_len:
                            # Deletion: replace ref bases with alt
                            if pos_in_gene + ref_len <= len(gene_seq_list):
                                gene_seq_list[pos_in_gene:pos_in_gene+ref_len] = list(alt_base)
                
                gene_seq = ''.join(gene_seq_list)
            else:
                # No mutations in this gene
                gene_seq = ref_gene_seq
            
            # Translate
            protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
            
            # Collect ALL non-synonymous mutations for this protein
            all_mutations = []
            all_mutation_ids = []
            max_freq = 0
            
            for mut in gene_muts:
                ref_aa, aa_pos, alt_aa = parse_aa_change(mut.get("HGVSp", ""))
                # Include non-synonymous substitutions and all indels
                is_indel = alt_aa and (alt_aa.startswith('ins') or alt_aa.startswith('del'))
                is_nonsynonymous_sub = ref_aa and alt_aa and ref_aa != alt_aa and not is_indel
                if is_indel or is_nonsynonymous_sub:  # Only non-synonymous or indels
                    # Format AA change appropriately
                    if aa_pos is not None:
                        aa_change = f"{ref_aa}{aa_pos}{alt_aa}"
                    else:
                        # For indels, ref_aa already contains position info (e.g., "R214_D215")
                        aa_change = f"{ref_aa}{alt_aa}"
                    mutation_id = f"{mut['POS']}{mut['REF']}>{mut['ALT']}"
                    all_mutations.append(aa_change)
                    all_mutation_ids.append(mutation_id)
                    if mut.get("Allele_Frequency", 0) > max_freq:
                        max_freq = mut.get("Allele_Frequency", 0)
            
            # Create ONE protein with ALL mutations
            if all_mutations:
                key = (gene, protein, tuple(all_mutations))
                all_proteins[key] = {
                    'sequence': protein,
                    'allele_ids': all_mutation_ids,
                    'mutations': all_mutations,
                    'freq': max_freq,  # Use highest frequency
                    'effect': 'missense_variant',
                    'gene': gene
                }
        else:
            # This gene has no mutations (or only synonymous) - use reference sequence
            gene_seq = ref_seq[start-1:end]  # Use reference (coordinates are stable)
            protein = viral_translate(gene_seq, coordinates=None, stop_at_stop_codon=False)
            
            key = (gene, protein, tuple())
            all_proteins[key] = {
                'sequence': protein,
                'allele_ids': [''],
                'mutations': [],
                'gene': gene
            }
    
    # Update summary data
    summary_data['protein_mutation_summary'] = dict(protein_mutation_summary)
    
    return all_proteins

def generate_summary_report(summary_data, output_prefix):
    """Generate comprehensive summary report"""
    report_file = f"{output_prefix}_summary_report.txt"
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("INDIVIDUAL VARIANT PROTEIN GENERATION SUMMARY REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        # Basic statistics
        f.write(f"Total mutations passing filters: {summary_data['total_mutations']}\n")
        f.write(f"Multi-allelic sites found: {len(summary_data.get('multiallelic_sites', []))}\n")
        f.write(f"Unique protein sequences: {summary_data['unique_proteins']}\n\n")
        
        # Protein mutation summary
        if summary_data['protein_mutation_summary']:
            f.write("-" * 80 + "\n")
            f.write("MUTATIONS PER PROTEIN\n")
            f.write("-" * 80 + "\n\n")
            
            for gene, mutations in summary_data['protein_mutation_summary'].items():
                f.write(f"{gene}:\n")
                if mutations:
                    for aa_change, occurrences in mutations.items():
                        f.write(f"  - {aa_change}:")
                        for occ in occurrences:
                            f.write(f" {occ['nucleotide']} ({occ['effect']}, {occ['freq']:.2%})")
                        f.write("\n")
                else:
                    f.write("  No mutations\n")
                f.write("\n")
        
        # Individual protein variant summary
        if summary_data.get('protein_dedup_info'):
            f.write("-" * 80 + "\n")
            f.write("INDIVIDUAL PROTEIN VARIANTS\n")
            f.write("-" * 80 + "\n\n")
            
            for info in summary_data['protein_dedup_info']:
                f.write(f"{info['gene']}: ")
                if info['mutations']:
                    f.write(f"Variant with mutations: {', '.join(info['mutations'])}\n")
                else:
                    f.write("Reference sequence (no mutations)\n")
                if info['allele_ids'] and info['allele_ids'][0]:
                    f.write(f"  From nucleotide change: {info['allele_ids'][0]}\n")
                f.write("\n")
    
    print(f"\n📄 Summary report saved to: {report_file}")

def main():
    parser = argparse.ArgumentParser(description='Generate individual variant proteins for quasispecies analysis')
    parser.add_argument('--vcf', required=True, help='Input VCF or filtered TSV file')
    parser.add_argument('--reference', required=True, help='Reference genome FASTA file')
    parser.add_argument('--accession', required=True, help='Virus accession number')
    parser.add_argument('--quality', type=float, default=1000, help='Minimum quality score (default: 1000)')
    parser.add_argument('--depth', type=int, default=200, help='Minimum depth (default: 200)')
    parser.add_argument('--freq', type=float, default=0.01, help='Minimum allele frequency (default: 0.01)')
    parser.add_argument('--output-prefix', required=True, help='Output file prefix')
    
    args = parser.parse_args()
    
    print("=" * 80)
    print("INDIVIDUAL VARIANT PROTEIN GENERATOR")
    print("=" * 80)
    print(f"Input VCF/TSV: {args.vcf}")
    print(f"Reference: {args.reference}")
    print(f"Accession: {args.accession}")
    print(f"Quality filter: >= {args.quality}")
    print(f"Depth filter: >= {args.depth}")
    print(f"Frequency filter: >= {args.freq}")
    print(f"Output prefix: {args.output_prefix}")
    print("=" * 80)
    
    # Read reference sequence
    ref_record = next(SeqIO.parse(args.reference, "fasta"))
    ref_seq = str(ref_record.seq)
    print(f"Reference genome length: {len(ref_seq)} bp")
    
    # Read virus configuration
    virus_config = read_virus_config(args.accession)
    if not virus_config:
        print(f"Warning: No configuration found for {args.accession}, will generate polyprotein only")
    else:
        print(f"Found configuration for {virus_config['name']}")
        gene_count = len(virus_config.get('gene_coords', {}))
        print(f"Genes defined: {gene_count}")
    
    # Filter mutations
    mutations_df = filter_vcf(args.vcf, args.quality, args.depth, args.freq)
    
    if len(mutations_df) == 0:
        print("No mutations pass the filtering criteria")
        return
    
    # Initialize summary data
    summary_data = {
        'total_mutations': len(mutations_df),
        'multiallelic_sites': [],
        'unique_proteins': 0,
        'protein_mutation_summary': {},
        'protein_dedup_info': []
    }
    
    # Generate individual variant proteins
    print(f"\n🧬 Generating individual variant proteins...")
    unique_proteins = generate_individual_variant_proteins(ref_seq, mutations_df, virus_config, summary_data)
    
    print(f"Generated {len(unique_proteins)} unique protein variants")
    summary_data['unique_proteins'] = len(unique_proteins)
    
    # Write consensus genome (with all mutations applied)
    print(f"\n📝 Writing consensus genome...")
    full_consensus = apply_mutations_single_sequence(ref_seq, mutations_df, "all_mutations")
    sample_name = Path(args.output_prefix).name
    
    consensus_record = SeqRecord(
        Seq(full_consensus),
        id=f"{sample_name}_filtered_consensus",
        description=f"Consensus with {len(mutations_df)} mutations (Q>={args.quality}, D>={args.depth}, F>={args.freq}) | Mutations: " + ", ".join([f"{row['POS']}{row['REF']}>{row['ALT']} ({row['Allele_Frequency']*100:.1f}%)" for _, row in mutations_df.iterrows()])
    )
    
    consensus_file = f"{args.output_prefix}.fasta"
    SeqIO.write([consensus_record], consensus_file, "fasta")
    print(f"Wrote consensus genome to: {consensus_file}")
    
    # Write protein sequences
    print(f"\n🧬 Writing protein sequences...")
    protein_records = []
    
    for key, protein_data in unique_proteins.items():
        if len(key) == 2:
            # Polyprotein case: (gene_info, protein_seq)
            gene_info, protein_seq = key
            mutations = tuple()  # No mutations for polyprotein
        elif len(key) == 3:
            # Gene-specific case: (gene_name, protein_seq, mutations_tuple)
            gene_info, protein_seq, mutations = key
        else:
            print(f"Warning: Unexpected key format: {key}")
            continue
        gene = protein_data.get('gene', gene_info)
        
        # Create description
        if protein_data['mutations']:
            mutation_str = ", ".join(protein_data['mutations'])
            desc = f"{gene} protein [{mutation_str}]"
        else:
            desc = f"{gene} protein"
        
        # Add nucleotide change info if available
        if protein_data['allele_ids'] and any(protein_data['allele_ids']):
            nuc_changes = [nuc for nuc in protein_data['allele_ids'] if nuc]
            if nuc_changes:
                if len(nuc_changes) == 1:
                    desc += f" (from {nuc_changes[0]})"
                else:
                    desc += f" (from {', '.join(nuc_changes)})"
        
        # Add allele frequency if available
        if 'freq' in protein_data and protein_data['freq'] > 0:
            freq_pct = protein_data['freq'] * 100
            desc += f", {freq_pct:.2f}%"
        
        
        record = SeqRecord(
            Seq(protein_data['sequence']),
            id=f"{sample_name}_{gene}",
            description=desc
        )
        protein_records.append(record)
        
        # Add to dedup info
        summary_data['protein_dedup_info'].append({
            'gene': gene,
            'allele_ids': protein_data['allele_ids'],
            'mutations': protein_data['mutations']
        })
    
    protein_file = f"{args.output_prefix}_proteins.fasta"
    SeqIO.write(protein_records, protein_file, "fasta")
    print(f"Wrote {len(protein_records)} protein sequences to: {protein_file}")
    
    # Generate summary report
    generate_summary_report(summary_data, args.output_prefix)
    
    print(f"\n✅ Individual variant protein generation complete!")
    print(f"Generated {len(unique_proteins)} unique protein variants")

if __name__ == "__main__":
    main()