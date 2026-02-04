#!/usr/bin/env python3
"""
Add evolutionary conservation scores to variant TSV files.

This script annotates variant TSV files (e.g., from SnpEff) with conservation
scores derived from multiple sequence alignments (MSAs). It supports:
- ColabFold/AlphaFold A3M format
- Standard MSA formats (Stockholm, Clustal, FASTA)
- Gene-specific MSA files
- Batch processing of multiple files

Usage:
    # Basic usage with single MSA
    python add_conservation_scores.py variants.snpEff.ann.tsv protein.a3m -o output.tsv

    # ColabFold directory with multiple MSAs
    python add_conservation_scores.py variants.tsv --colabfold-dir ./msas/

    # Gene-specific MSA mapping
    python add_conservation_scores.py variants.tsv --gene-msas "NS1:ns1.a3m,E:e.a3m"

    # Batch processing
    python add_conservation_scores.py *.tsv --msa protein.a3m --output-dir ./annotated/

Output columns added:
    CONSERVATION_SCORE    Normalized 0-1 (1 = highly conserved)
    SHANNON_ENTROPY       Lower = more conserved
    PERCENT_IDENTITY      % matching consensus residue
    CONSENSUS_AA          Most common amino acid at position
    MSA_DEPTH             Number of sequences in alignment
    GAP_FRACTION          Fraction of gaps at position
    CONSERVATION_CATEGORY Category: highly_conserved/conserved/moderately_conserved/variable
"""

import argparse
import sys
from pathlib import Path

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pandas as pd

from vicast.conservation import (
    annotate_tsv_file,
    annotate_variants_with_conservation,
    annotate_with_gene_msas,
    batch_annotate,
    get_conservation_summary,
    get_gene_msa_mapping,
    load_and_prepare_msa,
)


def find_colabfold_msa(colabfold_dir: Path, gene: str = None) -> Path:
    """Find MSA file in a ColabFold output directory.

    ColabFold typically creates files like:
    - protein_name.a3m
    - protein_name_env/protein_name.a3m

    Args:
        colabfold_dir: Path to ColabFold output directory
        gene: Optional gene name to look for

    Returns:
        Path to the MSA file

    Raises:
        FileNotFoundError: If no MSA found
    """
    # Try common patterns
    patterns = [
        "*.a3m",
        "*/*.a3m",
        "**/*.a3m",
    ]

    for pattern in patterns:
        matches = list(colabfold_dir.glob(pattern))
        if matches:
            if gene:
                # Try to find one matching the gene name
                for m in matches:
                    if gene.lower() in m.stem.lower():
                        return m
            # Return first match
            return matches[0]

    raise FileNotFoundError(f"No A3M file found in {colabfold_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Add conservation scores to variant TSV files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file with single MSA
  python add_conservation_scores.py variants.tsv protein.a3m -o output.tsv

  # Batch processing
  python add_conservation_scores.py *.tsv --msa protein.a3m --output-dir ./annotated/

  # Gene-specific MSAs
  python add_conservation_scores.py variants.tsv --gene-msas "NS1:ns1.a3m,E:e.a3m"

  # ColabFold directory
  python add_conservation_scores.py variants.tsv --colabfold-dir ./colabfold_output/
        """,
    )

    # Input files
    parser.add_argument(
        "input_files",
        nargs="+",
        help="Input TSV file(s) with variant annotations",
    )

    # MSA source options (use --msa for clarity, or the last positional arg if it's an MSA file)
    parser.add_argument(
        "--msa",
        dest="msa_file",
        help="MSA file path (A3M, Stockholm, Clustal, or FASTA)",
    )
    parser.add_argument(
        "--colabfold-dir",
        type=Path,
        help="ColabFold output directory containing MSA files",
    )
    parser.add_argument(
        "--gene-msas",
        help='Gene-specific MSA mapping as "GENE1:path1.a3m,GENE2:path2.a3m"',
    )

    # Output options
    parser.add_argument(
        "-o", "--output",
        help="Output file path (for single input file)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory for batch processing",
    )
    parser.add_argument(
        "--suffix",
        default=".conserv",
        help="Suffix to add to output filenames (default: .conserv)",
    )

    # Column names
    parser.add_argument(
        "--hgvsp-column",
        default="HGVSp",
        help="Name of column containing HGVSp notation (default: HGVSp)",
    )
    parser.add_argument(
        "--gene-column",
        default="GENE_NAME",
        help="Name of column containing gene names (default: GENE_NAME)",
    )

    # Other options
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Print progress messages",
    )
    parser.add_argument(
        "--summary",
        action="store_true",
        help="Print MSA conservation summary statistics",
    )

    args = parser.parse_args()

    # Resolve MSA source - check if last input file looks like an MSA
    msa_file = args.msa_file
    input_files = list(args.input_files)  # Make a copy we can modify

    # If last file looks like an MSA file and no --msa specified, use it as MSA
    if not msa_file and not args.colabfold_dir and not args.gene_msas:
        if len(input_files) >= 2:
            last_file = input_files[-1]
            if any(last_file.lower().endswith(ext) for ext in ['.a3m', '.sto', '.stockholm', '.aln', '.clustal']):
                msa_file = last_file
                input_files = input_files[:-1]

    if not any([msa_file, args.colabfold_dir, args.gene_msas]):
        parser.error("Must specify an MSA source: --msa FILE, --colabfold-dir, or --gene-msas")

    try:
        # Handle different MSA source modes
        if args.gene_msas:
            # Gene-specific MSA mode
            gene_mapping = get_gene_msa_mapping(args.gene_msas)

            if args.verbose:
                print(f"Using gene-specific MSAs:")
                for gene, path in gene_mapping.items():
                    print(f"  {gene}: {path}")

            for input_file in input_files:
                input_path = Path(input_file)
                if args.verbose:
                    print(f"\nProcessing: {input_path}")

                df = pd.read_csv(input_path, sep="\t")
                annotated = annotate_with_gene_msas(
                    df,
                    gene_mapping,
                    gene_column=args.gene_column,
                    hgvsp_column=args.hgvsp_column,
                    verbose=args.verbose,
                )

                # Determine output path
                if args.output and len(input_files) == 1:
                    output_path = args.output
                elif args.output_dir:
                    args.output_dir.mkdir(parents=True, exist_ok=True)
                    output_path = args.output_dir / f"{input_path.stem}{args.suffix}.tsv"
                else:
                    output_path = input_path.parent / f"{input_path.stem}{args.suffix}.tsv"

                annotated.to_csv(output_path, sep="\t", index=False)
                if args.verbose:
                    print(f"Saved to: {output_path}")

        elif args.colabfold_dir:
            # ColabFold directory mode
            msa_path = find_colabfold_msa(args.colabfold_dir)

            if args.verbose:
                print(f"Found MSA: {msa_path}")

            msa, validation = load_and_prepare_msa(msa_path, verbose=args.verbose)

            if args.summary:
                summary = get_conservation_summary(msa)
                print("\nMSA Conservation Summary:")
                for key, value in summary.items():
                    if isinstance(value, float):
                        print(f"  {key}: {value:.4f}")
                    else:
                        print(f"  {key}: {value}")

            # Process files
            for input_file in input_files:
                input_path = Path(input_file)
                if args.verbose:
                    print(f"\nProcessing: {input_path}")

                df = pd.read_csv(input_path, sep="\t")
                annotated = annotate_variants_with_conservation(
                    df, msa, hgvsp_column=args.hgvsp_column, verbose=args.verbose
                )

                # Determine output
                if args.output and len(input_files) == 1:
                    output_path = args.output
                elif args.output_dir:
                    args.output_dir.mkdir(parents=True, exist_ok=True)
                    output_path = args.output_dir / f"{input_path.stem}{args.suffix}.tsv"
                else:
                    output_path = input_path.parent / f"{input_path.stem}{args.suffix}.tsv"

                annotated.to_csv(output_path, sep="\t", index=False)
                if args.verbose:
                    print(f"Saved to: {output_path}")

        else:
            # Single MSA file mode
            msa_path = Path(msa_file)

            if len(input_files) == 1 and args.output:
                # Single file with explicit output
                annotate_tsv_file(
                    input_files[0],
                    msa_path,
                    output_tsv=args.output,
                    hgvsp_column=args.hgvsp_column,
                    verbose=args.verbose,
                )

                if args.summary:
                    msa, _ = load_and_prepare_msa(msa_path)
                    summary = get_conservation_summary(msa)
                    print("\nMSA Conservation Summary:")
                    for key, value in summary.items():
                        if isinstance(value, float):
                            print(f"  {key}: {value:.4f}")
                        else:
                            print(f"  {key}: {value}")

            else:
                # Batch mode
                results = batch_annotate(
                    input_files,
                    msa_path,
                    output_dir=args.output_dir,
                    suffix=args.suffix,
                    hgvsp_column=args.hgvsp_column,
                    verbose=args.verbose,
                )

                if args.summary:
                    msa, _ = load_and_prepare_msa(msa_path)
                    summary = get_conservation_summary(msa)
                    print("\nMSA Conservation Summary:")
                    for key, value in summary.items():
                        if isinstance(value, float):
                            print(f"  {key}: {value:.4f}")
                        else:
                            print(f"  {key}: {value}")

                # Print batch summary
                if args.verbose:
                    print(f"\nBatch processing complete:")
                    print(f"  Processed: {len(results)} files")
                    print(f"  Successful: {sum(1 for v in results.values() if v is not None)}")

        print("Done.")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {e}", file=sys.stderr)
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
