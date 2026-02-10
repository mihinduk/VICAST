"""
Variant annotation with conservation scores.

This module provides the high-level interface for adding conservation scores
to variant TSV files. It ties together MSA parsing, score calculation, and
position mapping to annotate variants.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

from .mapper import get_conservation_for_variant
from .models import ConservationResult, MSAlignment, MSAValidationResult
from .msa_parser import parse_msa, validate_msa
from .scores import calculate_conservation_scores, get_conservation_summary


def load_and_prepare_msa(
    msa_path: Union[str, Path],
    verbose: bool = False,
) -> Tuple[MSAlignment, MSAValidationResult]:
    """Load, validate, and calculate scores for an MSA file.

    This is the main entry point for preparing an MSA for conservation annotation.
    It performs:
    1. File parsing (auto-detect format)
    2. Validation
    3. Score calculation for all positions
    4. Query-to-MSA position mapping

    Args:
        msa_path: Path to the MSA file
        verbose: If True, print progress messages

    Returns:
        Tuple of (MSAlignment, MSAValidationResult)

    Raises:
        FileNotFoundError: If MSA file doesn't exist
        ValueError: If MSA is invalid for conservation scoring
    """
    if verbose:
        print(f"Loading MSA from: {msa_path}")

    # Parse MSA
    msa = parse_msa(msa_path)

    if verbose:
        print(f"  Parsed {msa.num_sequences} sequences, alignment length {msa.alignment_length}")

    # Validate
    validation = validate_msa(msa)

    if verbose:
        if validation.errors:
            print(f"  Validation errors: {validation.errors}")
        if validation.warnings:
            print(f"  Validation warnings: {validation.warnings}")

    if not validation.is_valid:
        raise ValueError(f"MSA validation failed: {validation.errors}")

    # Calculate conservation scores
    if verbose:
        print("  Calculating conservation scores...")

    calculate_conservation_scores(msa)

    if verbose:
        summary = get_conservation_summary(msa)
        print(f"  Mean conservation: {summary.get('mean_conservation', 0):.3f}")
        print(f"  Highly conserved positions: {summary.get('highly_conserved_positions', 0)}")

    return msa, validation


def annotate_variants_with_conservation(
    variants_df: pd.DataFrame,
    msa: MSAlignment,
    hgvsp_column: str = "HGVSp",
    verbose: bool = False,
) -> pd.DataFrame:
    """Add conservation score columns to a variants DataFrame.

    Args:
        variants_df: DataFrame with variant annotations
        msa: Prepared MSAlignment with calculated scores
        hgvsp_column: Name of column containing HGVSp notation
        verbose: If True, print progress messages

    Returns:
        DataFrame with added conservation columns
    """
    # Make a copy to avoid modifying original
    df = variants_df.copy()

    # Check if HGVSp column exists
    if hgvsp_column not in df.columns:
        # Try common alternatives
        alt_columns = ["hgvsp", "HGVSP", "HGVSp", "protein_change", "AA_change"]
        for alt in alt_columns:
            if alt in df.columns:
                hgvsp_column = alt
                break
        else:
            if verbose:
                print(f"Warning: HGVSp column '{hgvsp_column}' not found. Conservation columns will be empty.")
            # Add empty columns
            for col in ["CONSERVATION_SCORE", "SHANNON_ENTROPY", "PERCENT_IDENTITY",
                       "CONSENSUS_AA", "MSA_DEPTH", "GAP_FRACTION", "CONSERVATION_CATEGORY"]:
                df[col] = None
            return df

    if verbose:
        print(f"Annotating {len(df)} variants...")
        annotated_count = 0

    # Initialize new columns
    new_columns = {
        "CONSERVATION_SCORE": [],
        "SHANNON_ENTROPY": [],
        "PERCENT_IDENTITY": [],
        "CONSENSUS_AA": [],
        "MSA_DEPTH": [],
        "GAP_FRACTION": [],
        "CONSERVATION_CATEGORY": [],
    }

    # Process each variant
    for idx, row in df.iterrows():
        hgvsp = row.get(hgvsp_column, "")

        if pd.isna(hgvsp) or hgvsp in ("", ".", "-", "NA"):
            # No protein change - add None values
            result = ConservationResult(protein_position=0, found_in_msa=False)
        else:
            result = get_conservation_for_variant(str(hgvsp), msa)

        # Convert to dict and add to columns
        result_dict = result.to_dict()
        for col, values in new_columns.items():
            values.append(result_dict.get(col))

        if result.found_in_msa and verbose:
            annotated_count += 1

    # Add columns to DataFrame
    for col, values in new_columns.items():
        df[col] = values

    if verbose:
        print(f"  Successfully annotated {annotated_count}/{len(df)} variants")

    return df


def annotate_tsv_file(
    input_tsv: Union[str, Path],
    msa_path: Union[str, Path],
    output_tsv: Optional[Union[str, Path]] = None,
    hgvsp_column: str = "HGVSp",
    verbose: bool = False,
) -> pd.DataFrame:
    """Annotate a TSV file with conservation scores.

    Convenience function that handles loading both files, annotating,
    and optionally saving the output.

    Args:
        input_tsv: Path to input TSV file with variants
        msa_path: Path to MSA file
        output_tsv: Optional output path. If None, returns DataFrame without saving
        hgvsp_column: Name of column containing HGVSp notation
        verbose: If True, print progress messages

    Returns:
        Annotated DataFrame
    """
    if verbose:
        print(f"Loading variants from: {input_tsv}")

    # Load variants
    df = pd.read_csv(input_tsv, sep="\t")

    if verbose:
        print(f"  Found {len(df)} variants")

    # Load and prepare MSA
    msa, validation = load_and_prepare_msa(msa_path, verbose=verbose)

    # Annotate
    annotated_df = annotate_variants_with_conservation(
        df, msa, hgvsp_column=hgvsp_column, verbose=verbose
    )

    # Save if output path provided
    if output_tsv:
        annotated_df.to_csv(output_tsv, sep="\t", index=False)
        if verbose:
            print(f"Saved annotated variants to: {output_tsv}")

    return annotated_df


def batch_annotate(
    input_files: List[Union[str, Path]],
    msa_path: Union[str, Path],
    output_dir: Optional[Union[str, Path]] = None,
    suffix: str = ".conserv",
    hgvsp_column: str = "HGVSp",
    verbose: bool = False,
) -> Dict[str, pd.DataFrame]:
    """Annotate multiple TSV files with conservation scores.

    Args:
        input_files: List of input TSV file paths
        msa_path: Path to MSA file (used for all files)
        output_dir: Directory for output files. If None, outputs go alongside inputs
        suffix: Suffix to add before .tsv extension
        hgvsp_column: Name of column containing HGVSp notation
        verbose: If True, print progress messages

    Returns:
        Dictionary mapping input paths to annotated DataFrames
    """
    results = {}

    # Load MSA once
    if verbose:
        print("Loading MSA for batch processing...")
    msa, validation = load_and_prepare_msa(msa_path, verbose=verbose)

    for input_path in input_files:
        input_path = Path(input_path)

        if verbose:
            print(f"\nProcessing: {input_path.name}")

        try:
            # Load variants
            df = pd.read_csv(input_path, sep="\t")

            # Annotate
            annotated_df = annotate_variants_with_conservation(
                df, msa, hgvsp_column=hgvsp_column, verbose=verbose
            )

            # Determine output path
            if output_dir:
                output_dir = Path(output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                output_path = output_dir / f"{input_path.stem}{suffix}.tsv"
            else:
                output_path = input_path.parent / f"{input_path.stem}{suffix}.tsv"

            # Save
            annotated_df.to_csv(output_path, sep="\t", index=False)
            if verbose:
                print(f"  Saved to: {output_path}")

            results[str(input_path)] = annotated_df

        except Exception as e:
            if verbose:
                print(f"  Error processing {input_path}: {e}")
            results[str(input_path)] = None

    return results


def get_gene_msa_mapping(gene_msas_str: str) -> Dict[str, Path]:
    """Parse gene-to-MSA mapping from command line argument.

    Args:
        gene_msas_str: Comma-separated "GENE:path.a3m" pairs

    Returns:
        Dictionary mapping gene names to MSA file paths
    """
    mapping = {}

    for pair in gene_msas_str.split(","):
        pair = pair.strip()
        if ":" in pair:
            gene, path = pair.split(":", 1)
            mapping[gene.strip()] = Path(path.strip())

    return mapping


def annotate_with_gene_msas(
    variants_df: pd.DataFrame,
    gene_msa_mapping: Dict[str, Path],
    gene_column: str = "GENE_NAME",
    hgvsp_column: str = "HGVSp",
    verbose: bool = False,
) -> pd.DataFrame:
    """Annotate variants using gene-specific MSA files.

    For each variant, looks up the gene-specific MSA and calculates
    conservation scores using that MSA.

    Args:
        variants_df: DataFrame with variant annotations
        gene_msa_mapping: Dictionary mapping gene names to MSA paths
        gene_column: Column containing gene names
        hgvsp_column: Column containing HGVSp notation
        verbose: If True, print progress messages

    Returns:
        DataFrame with conservation columns added
    """
    df = variants_df.copy()

    # Initialize columns
    for col in ["CONSERVATION_SCORE", "SHANNON_ENTROPY", "PERCENT_IDENTITY",
               "CONSENSUS_AA", "MSA_DEPTH", "GAP_FRACTION", "CONSERVATION_CATEGORY"]:
        df[col] = None

    # Load MSAs for each gene
    msa_cache: Dict[str, MSAlignment] = {}

    if gene_column not in df.columns:
        if verbose:
            print(f"Warning: Gene column '{gene_column}' not found")
        return df

    # Process variants grouped by gene
    for gene in df[gene_column].unique():
        if pd.isna(gene):
            continue

        gene_str = str(gene)

        # Check if we have an MSA for this gene
        if gene_str not in gene_msa_mapping:
            if verbose:
                print(f"  No MSA for gene: {gene_str}")
            continue

        msa_path = gene_msa_mapping[gene_str]

        # Load MSA (with caching)
        if gene_str not in msa_cache:
            try:
                msa, _ = load_and_prepare_msa(msa_path, verbose=False)
                msa_cache[gene_str] = msa
                if verbose:
                    print(f"  Loaded MSA for {gene_str}: {msa.num_sequences} sequences")
            except Exception as e:
                if verbose:
                    print(f"  Failed to load MSA for {gene_str}: {e}")
                continue

        msa = msa_cache[gene_str]

        # Annotate variants for this gene
        gene_mask = df[gene_column] == gene
        gene_indices = df[gene_mask].index

        for idx in gene_indices:
            hgvsp = df.loc[idx, hgvsp_column] if hgvsp_column in df.columns else None

            if pd.isna(hgvsp) or hgvsp in ("", ".", "-", "NA"):
                continue

            result = get_conservation_for_variant(str(hgvsp), msa)
            result_dict = result.to_dict()

            for col, value in result_dict.items():
                df.loc[idx, col] = value

    return df
