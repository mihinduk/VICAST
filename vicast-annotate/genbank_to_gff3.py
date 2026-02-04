#!/usr/bin/env python3
"""
Convert GenBank files to GFF3 format for VICAST/SnpEff.

Creates properly formatted GFF3 with gene/mRNA/CDS hierarchy that is
compatible with SnpEff database building.

Features:
- Handles spliced genes (join locations)
- Proper phase calculation for CDS features
- Unique IDs for all features
- Validates output GFF3 if VICAST is installed

Usage:
    python genbank_to_gff3.py input1.gb [input2.gb ...] -o output.gff3
    python genbank_to_gff3.py --fasta-out combined.fasta input*.gb -o output.gff3

Examples:
    # Convert single GenBank file
    python genbank_to_gff3.py NC_001477.gb -o dengue.gff3

    # Convert multiple files and create combined FASTA
    python genbank_to_gff3.py segment*.gb -o flu.gff3 --fasta-out flu.fasta

    # With validation
    python genbank_to_gff3.py NC_001477.gb -o dengue.gff3 --validate
"""

import argparse
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


def convert_genbank_to_gff3(
    gb_files: list,
    output_gff: str,
    output_fasta: str = None,
    source: str = "NCBI",
    verbose: bool = False
) -> tuple:
    """
    Convert GenBank files to GFF3 format.

    Args:
        gb_files: List of GenBank file paths
        output_gff: Output GFF3 file path
        output_fasta: Optional output FASTA file path for combined sequences
        source: Source field for GFF3 (default: NCBI)
        verbose: Print progress messages

    Returns:
        Tuple of (gene_count, feature_count)
    """
    gene_count = 0
    feature_count = 0
    sequences = []

    with open(output_gff, 'w') as out:
        out.write("##gff-version 3\n")

        for gb_file in sorted(gb_files):
            if verbose:
                print(f"Processing: {gb_file}")

            for record in SeqIO.parse(gb_file, "genbank"):
                seqid = record.id
                seq_len = len(record.seq)

                # Store sequence for FASTA output
                if output_fasta:
                    sequences.append(record)

                # Write sequence region
                out.write(f"##sequence-region {seqid} 1 {seq_len}\n")

                # Track gene IDs to avoid duplicates
                id_counts = {}

                for feature in record.features:
                    if feature.type == "CDS":
                        gene_count += 1

                        # Extract info
                        strand = '+' if feature.location.strand == 1 else '-'

                        # Get gene name and product
                        gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
                        product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                        protein_id = feature.qualifiers.get('protein_id', [''])[0]

                        # Escape special characters in product
                        product = product.replace(';', '%3B').replace('=', '%3D')

                        # Handle duplicates (e.g., M1/M2 on same segment)
                        if gene_name in id_counts:
                            id_counts[gene_name] += 1
                            unique_id = f"{gene_name}_{id_counts[gene_name]}"
                        else:
                            id_counts[gene_name] = 1
                            unique_id = gene_name

                        # Handle spliced features (join locations)
                        if hasattr(feature.location, 'parts'):
                            # Spliced gene - write each exon
                            parts = list(feature.location.parts)

                            # Gene feature spans full range
                            gene_start = min(int(p.start) + 1 for p in parts)
                            gene_end = max(int(p.end) for p in parts)

                            out.write(f"{seqid}\t{source}\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
                            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")
                            feature_count += 1

                            out.write(f"{seqid}\t{source}\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
                            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")
                            feature_count += 1

                            # Write CDS parts
                            cumulative_len = 0
                            for i, part in enumerate(parts):
                                p_start = int(part.start) + 1
                                p_end = int(part.end)

                                # Calculate phase based on cumulative length
                                phase = (3 - (cumulative_len % 3)) % 3
                                cumulative_len += (p_end - p_start + 1)

                                out.write(f"{seqid}\t{source}\tCDS\t{p_start}\t{p_end}\t.\t{strand}\t{phase}\t")
                                out.write(f"ID=cds_{unique_id}_{i+1};Parent=mRNA_{unique_id};")
                                out.write(f"Name={gene_name};gene={gene_name};product={product}")
                                if protein_id:
                                    out.write(f";protein_id={protein_id}")
                                out.write("\n")
                                feature_count += 1
                        else:
                            # Simple unspliced gene
                            start = int(feature.location.start) + 1  # 1-based
                            end = int(feature.location.end)

                            out.write(f"{seqid}\t{source}\tgene\t{start}\t{end}\t.\t{strand}\t.\t")
                            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")
                            feature_count += 1

                            out.write(f"{seqid}\t{source}\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t")
                            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")
                            feature_count += 1

                            out.write(f"{seqid}\t{source}\tCDS\t{start}\t{end}\t.\t{strand}\t0\t")
                            out.write(f"ID=cds_{unique_id};Parent=mRNA_{unique_id};")
                            out.write(f"Name={gene_name};gene={gene_name};product={product}")
                            if protein_id:
                                out.write(f";protein_id={protein_id}")
                            out.write("\n")
                            feature_count += 1

    # Write combined FASTA if requested
    if output_fasta and sequences:
        SeqIO.write(sequences, output_fasta, "fasta")
        if verbose:
            print(f"Created FASTA: {output_fasta} ({len(sequences)} sequences)")

    return gene_count, feature_count


def validate_gff3(gff_path: str, fasta_path: str = None) -> bool:
    """
    Validate GFF3 file using VICAST validation module.

    Args:
        gff_path: Path to GFF3 file
        fasta_path: Optional path to reference FASTA

    Returns:
        True if valid, False otherwise
    """
    try:
        from vicast.validation import validate_gff_for_snpeff
        is_valid, errors, warnings = validate_gff_for_snpeff(gff_path, fasta_path)

        for w in warnings:
            print(f"  WARN: {w}")
        for e in errors:
            print(f"  ERROR: {e}")

        return is_valid
    except ImportError:
        print("  VICAST validation module not available, skipping validation")
        return True


def main():
    parser = argparse.ArgumentParser(
        description="Convert GenBank files to GFF3 format for SnpEff",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Convert single GenBank file
  %(prog)s NC_001477.gb -o dengue.gff3

  # Convert multiple files and create combined FASTA
  %(prog)s segment*.gb -o flu.gff3 --fasta-out flu.fasta

  # With validation
  %(prog)s NC_001477.gb -o dengue.gff3 --validate
"""
    )

    parser.add_argument(
        'genbank_files',
        nargs='+',
        help='Input GenBank file(s)'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output GFF3 file path'
    )
    parser.add_argument(
        '--fasta-out',
        help='Also create combined FASTA file for SnpEff'
    )
    parser.add_argument(
        '--source',
        default='NCBI',
        help='Source field for GFF3 (default: NCBI)'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate output GFF3 using VICAST'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print progress messages'
    )

    args = parser.parse_args()

    # Check input files exist
    for gb_file in args.genbank_files:
        if not Path(gb_file).exists():
            print(f"Error: File not found: {gb_file}")
            sys.exit(1)

    # Convert
    gene_count, feature_count = convert_genbank_to_gff3(
        gb_files=args.genbank_files,
        output_gff=args.output,
        output_fasta=args.fasta_out,
        source=args.source,
        verbose=args.verbose
    )

    print(f"Created: {args.output}")
    print(f"  Genes: {gene_count}")
    print(f"  Features: {feature_count}")

    # Validate if requested
    if args.validate:
        print("\nValidating GFF3...")
        if validate_gff3(args.output, args.fasta_out):
            print("  Validation PASSED")
        else:
            print("  Validation FAILED")
            sys.exit(1)


if __name__ == "__main__":
    main()
