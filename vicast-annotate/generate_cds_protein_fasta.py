#!/usr/bin/env python3
"""
Generate CDS and protein FASTA files for SnpEff validation.
This script extracts CDS sequences from genome and translates them to proteins.
"""

import sys
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd

def generate_cds_protein_from_tsv(genome_fasta, tsv_file, output_cds, output_protein):
    """
    Generate CDS and protein FASTA files from genome and TSV annotation.

    Args:
        genome_fasta: Path to genome FASTA file
        tsv_file: Path to TSV annotation file
        output_cds: Path to output CDS FASTA
        output_protein: Path to output protein FASTA

    Returns:
        tuple: (cds_count, protein_count)
    """
    # Load genome sequence
    genome_record = SeqIO.read(genome_fasta, "fasta")
    genome_seq = genome_record.seq

    print(f"Loaded genome: {genome_record.id}")
    print(f"Genome length: {len(genome_seq)} bp")

    # Load TSV annotations
    df = pd.read_csv(tsv_file, sep='\t')

    # Filter for CDS features that are marked KEEP
    cds_df = df[(df['type'] == 'CDS') & (df.get('action', 'KEEP') != 'DELETE')]

    print(f"\nFound {len(cds_df)} CDS features to extract")

    cds_sequences = []
    protein_sequences = []

    for idx, row in cds_df.iterrows():
        start = int(row['start']) - 1  # Convert to 0-based
        end = int(row['end'])
        strand = row['strand']
        gene_name = row.get('gene_name', '') or row.get('gene', '')
        protein_id = row.get('protein_id', '')
        product = row.get('product', '')

        # Extract CDS sequence
        cds_seq = genome_seq[start:end]

        # Reverse complement if on minus strand
        if strand == '-':
            cds_seq = cds_seq.reverse_complement()

        # Generate sequence ID using gene_name|protein_id format
        if gene_name and protein_id:
            seq_id = f"{gene_name}|{protein_id}"
        elif gene_name:
            seq_id = gene_name
        elif protein_id:
            seq_id = protein_id
        else:
            seq_id = f"CDS_{idx}"

        # Description
        description = product if product else "hypothetical protein"

        # Create CDS record
        cds_record = SeqIO.SeqRecord(
            cds_seq,
            id=seq_id,
            description=description
        )
        cds_sequences.append(cds_record)

        # Translate to protein
        try:
            # Try standard translation
            protein_seq = cds_seq.translate(to_stop=True)

            # If sequence doesn't start with M, try with alternate start codon
            if not str(protein_seq).startswith('M') and len(cds_seq) >= 3:
                # Some viral genes use alternate start codons
                # Force first codon to be Met
                protein_seq = Seq('M') + cds_seq[3:].translate(to_stop=True)

            protein_record = SeqIO.SeqRecord(
                protein_seq,
                id=seq_id,
                description=description
            )
            protein_sequences.append(protein_record)

        except Exception as e:
            print(f"  Warning: Could not translate {seq_id}: {e}")
            # Create empty protein record
            protein_record = SeqIO.SeqRecord(
                Seq(''),
                id=seq_id,
                description=f"{description} [translation failed]"
            )
            protein_sequences.append(protein_record)

    # Write CDS FASTA
    SeqIO.write(cds_sequences, output_cds, "fasta")
    print(f"\nWrote {len(cds_sequences)} CDS sequences to: {output_cds}")

    # Write protein FASTA
    SeqIO.write(protein_sequences, output_protein, "fasta")
    print(f"Wrote {len(protein_sequences)} protein sequences to: {output_protein}")

    return len(cds_sequences), len(protein_sequences)

def main():
    parser = argparse.ArgumentParser(
        description='Generate CDS and protein FASTA files for SnpEff validation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate from TSV and genome
  python3 generate_cds_protein_fasta.py NC_001477.fasta NC_001477_no_polyprotein.tsv

  # Specify output names
  python3 generate_cds_protein_fasta.py genome.fasta annotations.tsv --cds cds.fa --protein protein.fa

This script should be run after step1 and before step2 to generate
the CDS and protein FASTA files that SnpEff uses for validation.

FASTA header format: gene_name|protein_id (e.g., >ancC|NP_722457.2)
        """
    )

    parser.add_argument('genome_fasta',
                       help='Path to genome FASTA file')
    parser.add_argument('tsv_file',
                       help='Path to TSV annotation file')
    parser.add_argument('--cds',
                       help='Output CDS FASTA file (default: genome_id_cds.fa)')
    parser.add_argument('--protein',
                       help='Output protein FASTA file (default: genome_id_protein.fa)')

    args = parser.parse_args()

    # Check input files exist
    if not os.path.exists(args.genome_fasta):
        print(f"Error: Genome FASTA not found: {args.genome_fasta}")
        sys.exit(1)

    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV file not found: {args.tsv_file}")
        sys.exit(1)

    # Determine output file names
    if args.cds:
        output_cds = args.cds
    else:
        base = os.path.splitext(args.genome_fasta)[0]
        output_cds = f"{base}_cds.fa"

    if args.protein:
        output_protein = args.protein
    else:
        base = os.path.splitext(args.genome_fasta)[0]
        output_protein = f"{base}_protein.fa"

    print("="*60)
    print("Generate CDS and Protein FASTA for SnpEff Validation")
    print("="*60)
    print(f"\nGenome FASTA: {args.genome_fasta}")
    print(f"TSV annotations: {args.tsv_file}")
    print(f"Output CDS: {output_cds}")
    print(f"Output protein: {output_protein}")
    print()

    # Generate files
    cds_count, protein_count = generate_cds_protein_from_tsv(
        args.genome_fasta,
        args.tsv_file,
        output_cds,
        output_protein
    )

    print("\n" + "="*60)
    print("SUCCESS")
    print("="*60)
    print(f"\nGenerated {cds_count} CDS sequences")
    print(f"Generated {protein_count} protein sequences")
    print("\nThese files can now be used by SnpEff for validation.")
    print("Copy them to the SnpEff data directory:")
    print(f"  cp {output_cds} $SNPEFF_DATA/<genome_id>/cds.fa")
    print(f"  cp {output_protein} $SNPEFF_DATA/<genome_id>/protein.fa")

if __name__ == '__main__':
    main()
