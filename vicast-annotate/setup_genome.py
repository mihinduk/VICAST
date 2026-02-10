#!/usr/bin/env python3
"""
VICAST Genome Setup Tool

Sets up a SnpEff genome database from reference files.
Supports input from GenBank files or GFF3/FASTA files.

Usage:
    # From GenBank files
    vicast-setup-genome my_virus --genbank input.gb

    # From GFF3 and FASTA
    vicast-setup-genome my_virus --gff genes.gff3 --fasta genome.fasta

    # From GenBank with auto-validation
    vicast-setup-genome my_virus --genbank *.gb --validate

Examples:
    # Set up Dengue virus from NCBI GenBank
    vicast-setup-genome dengue_1 --genbank NC_001477.gb

    # Set up Influenza with multiple segments
    vicast-setup-genome influenza_A --genbank segment*.gb

    # Use existing GFF3 and FASTA
    vicast-setup-genome my_virus --gff genes.gff3 --fasta genome.fasta

    # Specify SnpEff installation directory
    vicast-setup-genome my_virus --genbank input.gb --snpeff-dir /path/to/snpEff
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    from Bio import SeqIO
except ImportError:
    SeqIO = None


def find_snpeff_dir():
    """Find SnpEff installation directory."""
    # Check environment variable
    if os.environ.get('SNPEFF_HOME'):
        return Path(os.environ['SNPEFF_HOME'])

    # Check common locations
    locations = [
        Path.home() / 'snpEff',
        Path('/opt/snpEff'),
        Path('/usr/local/snpEff'),
        Path(__file__).parent.parent / 'tools' / 'snpEff',
    ]

    for loc in locations:
        if (loc / 'snpEff.jar').exists():
            return loc

    return None


def convert_genbank_to_gff3(gb_files, output_dir, genome_name, verbose=False):
    """
    Convert GenBank files to GFF3 and FASTA for SnpEff.

    Args:
        gb_files: List of GenBank file paths
        output_dir: Output directory for files
        genome_name: Genome name for file naming
        verbose: Print progress

    Returns:
        Tuple of (gff_path, fasta_path)
    """
    if SeqIO is None:
        print("Error: BioPython required for GenBank conversion")
        print("  Install with: pip install biopython")
        sys.exit(1)

    gff_path = output_dir / "genes.gff"
    fasta_path = output_dir / "sequences.fa"

    sequences = []
    gene_count = 0

    with open(gff_path, 'w') as out:
        out.write("##gff-version 3\n")

        for gb_file in sorted(gb_files):
            if verbose:
                print(f"  Processing: {gb_file}")

            for record in SeqIO.parse(gb_file, "genbank"):
                seqid = record.id
                seq_len = len(record.seq)
                sequences.append(record)

                out.write(f"##sequence-region {seqid} 1 {seq_len}\n")

                id_counts = {}

                for feature in record.features:
                    if feature.type == "CDS":
                        gene_count += 1
                        strand = '+' if feature.location.strand == 1 else '-'

                        gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
                        product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                        protein_id = feature.qualifiers.get('protein_id', [''])[0]

                        product = product.replace(';', '%3B').replace('=', '%3D')

                        if gene_name in id_counts:
                            id_counts[gene_name] += 1
                            unique_id = f"{gene_name}_{id_counts[gene_name]}"
                        else:
                            id_counts[gene_name] = 1
                            unique_id = gene_name

                        if hasattr(feature.location, 'parts'):
                            parts = list(feature.location.parts)
                            gene_start = min(int(p.start) + 1 for p in parts)
                            gene_end = max(int(p.end) for p in parts)

                            out.write(f"{seqid}\tNCBI\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
                            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")

                            out.write(f"{seqid}\tNCBI\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
                            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")

                            cumulative_len = 0
                            for i, part in enumerate(parts):
                                p_start = int(part.start) + 1
                                p_end = int(part.end)
                                phase = (3 - (cumulative_len % 3)) % 3
                                cumulative_len += (p_end - p_start + 1)

                                out.write(f"{seqid}\tNCBI\tCDS\t{p_start}\t{p_end}\t.\t{strand}\t{phase}\t")
                                out.write(f"ID=cds_{unique_id}_{i+1};Parent=mRNA_{unique_id};")
                                out.write(f"Name={gene_name};gene={gene_name};product={product}")
                                if protein_id:
                                    out.write(f";protein_id={protein_id}")
                                out.write("\n")
                        else:
                            start = int(feature.location.start) + 1
                            end = int(feature.location.end)

                            out.write(f"{seqid}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\t")
                            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")

                            out.write(f"{seqid}\tNCBI\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t")
                            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")

                            out.write(f"{seqid}\tNCBI\tCDS\t{start}\t{end}\t.\t{strand}\t0\t")
                            out.write(f"ID=cds_{unique_id};Parent=mRNA_{unique_id};")
                            out.write(f"Name={gene_name};gene={gene_name};product={product}")
                            if protein_id:
                                out.write(f";protein_id={protein_id}")
                            out.write("\n")

    # Write FASTA
    SeqIO.write(sequences, fasta_path, "fasta")

    if verbose:
        print(f"  Created GFF3: {gff_path} ({gene_count} genes)")
        print(f"  Created FASTA: {fasta_path} ({len(sequences)} sequences)")

    return gff_path, fasta_path


def add_genome_to_config(snpeff_dir, genome_name, description=None):
    """Add genome to snpEff.config."""
    config_path = snpeff_dir / "snpEff.config"

    if description is None:
        description = genome_name

    # Check if already present
    with open(config_path, 'r') as f:
        if f"{genome_name}.genome" in f.read():
            return False

    # Add to config
    with open(config_path, 'a') as f:
        f.write(f"\n# VICAST custom genome\n")
        f.write(f"{genome_name}.genome : {description}\n")

    return True


def build_snpeff_database(snpeff_dir, genome_name, verbose=False):
    """Build SnpEff database."""
    jar_path = snpeff_dir / "snpEff.jar"

    cmd = ["java", "-jar", str(jar_path), "build", "-gff3", "-v", genome_name]

    if verbose:
        print(f"  Running: {' '.join(cmd)}")

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    # Check for success indicators
    if "Done" in result.stderr or result.returncode == 0:
        return True

    print("  Build output:")
    for line in result.stderr.split('\n'):
        if any(x in line for x in ['Error', 'Warning', 'Done', 'Protein', 'CDS']):
            print(f"    {line}")

    return result.returncode == 0


def validate_files(gff_path, fasta_path):
    """Validate GFF3 and FASTA files using VICAST validation."""
    try:
        from vicast.validation import validate_gff_for_snpeff
        is_valid, errors, warnings = validate_gff_for_snpeff(str(gff_path), str(fasta_path))

        for w in warnings:
            print(f"  WARN: {w}")
        for e in errors:
            print(f"  ERROR: {e}")

        return is_valid
    except ImportError:
        return True  # Skip validation if module not available


def main():
    parser = argparse.ArgumentParser(
        description="Set up a SnpEff genome database for VICAST",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # From GenBank files
  %(prog)s dengue_1 --genbank NC_001477.gb

  # From GFF3 and FASTA
  %(prog)s my_virus --gff genes.gff3 --fasta genome.fasta

  # Multiple GenBank files (segmented virus)
  %(prog)s influenza_A --genbank segment*.gb --description "Influenza A H1N1"
"""
    )

    parser.add_argument(
        'genome_name',
        help='Name for the genome database (e.g., dengue_1, influenza_A)'
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '--genbank',
        nargs='+',
        metavar='FILE',
        help='Input GenBank file(s)'
    )
    input_group.add_argument(
        '--gff',
        metavar='FILE',
        help='Input GFF3 annotation file (requires --fasta)'
    )

    parser.add_argument(
        '--fasta',
        metavar='FILE',
        help='Input FASTA reference file (required with --gff)'
    )
    parser.add_argument(
        '--description',
        help='Human-readable genome description'
    )
    parser.add_argument(
        '--snpeff-dir',
        metavar='DIR',
        help='SnpEff installation directory (default: auto-detect)'
    )
    parser.add_argument(
        '--validate',
        action='store_true',
        help='Validate files before building'
    )
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Print progress messages'
    )

    args = parser.parse_args()

    # Validate arguments
    if args.gff and not args.fasta:
        parser.error("--fasta is required when using --gff")

    # Find SnpEff
    if args.snpeff_dir:
        snpeff_dir = Path(args.snpeff_dir)
    else:
        snpeff_dir = find_snpeff_dir()

    if not snpeff_dir or not (snpeff_dir / "snpEff.jar").exists():
        print("Error: SnpEff not found")
        print("  Please install SnpEff or specify --snpeff-dir")
        print("  Or run: tools/setup_snpeff.sh --install-only")
        sys.exit(1)

    if args.verbose:
        print(f"Using SnpEff at: {snpeff_dir}")

    # Create data directory
    data_dir = snpeff_dir / "data" / args.genome_name
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Setting up genome: {args.genome_name}")

    # Get or create GFF3 and FASTA
    if args.genbank:
        # Check input files
        for gb_file in args.genbank:
            if not Path(gb_file).exists():
                print(f"Error: File not found: {gb_file}")
                sys.exit(1)

        print("Converting GenBank to GFF3/FASTA...")
        gff_path, fasta_path = convert_genbank_to_gff3(
            args.genbank,
            data_dir,
            args.genome_name,
            verbose=args.verbose
        )
    else:
        # Check input files
        if not Path(args.gff).exists():
            print(f"Error: GFF3 file not found: {args.gff}")
            sys.exit(1)
        if not Path(args.fasta).exists():
            print(f"Error: FASTA file not found: {args.fasta}")
            sys.exit(1)

        print("Copying reference files...")
        gff_path = data_dir / "genes.gff"
        fasta_path = data_dir / "sequences.fa"
        shutil.copy(args.gff, gff_path)
        shutil.copy(args.fasta, fasta_path)
        if args.verbose:
            print(f"  Copied GFF3: {gff_path}")
            print(f"  Copied FASTA: {fasta_path}")

    # Validate if requested
    if args.validate:
        print("Validating files...")
        if not validate_files(gff_path, fasta_path):
            print("Validation FAILED. Fix errors before proceeding.")
            sys.exit(1)
        print("  Validation PASSED")

    # Add to config
    print("Updating snpEff.config...")
    added = add_genome_to_config(
        snpeff_dir,
        args.genome_name,
        args.description or args.genome_name
    )
    if added:
        if args.verbose:
            print(f"  Added {args.genome_name} to config")
    else:
        if args.verbose:
            print(f"  {args.genome_name} already in config")

    # Build database
    print("Building SnpEff database...")
    if build_snpeff_database(snpeff_dir, args.genome_name, verbose=args.verbose):
        print("  Build completed")
    else:
        print("  Build may have warnings (check output above)")

    # Success
    print("")
    print("=" * 50)
    print(f"Genome '{args.genome_name}' is ready!")
    print("=" * 50)
    print("")
    print("To annotate variants:")
    print(f"  java -jar {snpeff_dir}/snpEff.jar {args.genome_name} input.vcf > output.vcf")
    print("")


if __name__ == "__main__":
    main()
