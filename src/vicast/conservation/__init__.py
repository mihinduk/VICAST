"""
VICAST Conservation Module: MSA-based evolutionary conservation scoring.

This module provides tools for annotating viral variants with evolutionary
conservation context from multiple sequence alignments (MSAs). It supports
ColabFold/AlphaFold A3M format as well as standard alignment formats.

Key Features:
- Parse MSA files (A3M, Stockholm, Clustal, FASTA)
- Calculate Shannon entropy and percent identity
- Map HGVSp protein positions to MSA columns
- Annotate variant TSV files with conservation scores

Example Usage:
    >>> from vicast.conservation import annotate_tsv_file
    >>> annotate_tsv_file(
    ...     "variants.snpEff.ann.tsv",
    ...     "protein.a3m",
    ...     output_tsv="variants_with_conservation.tsv",
    ...     verbose=True
    ... )

    # Or for more control:
    >>> from vicast.conservation import load_and_prepare_msa, annotate_variants_with_conservation
    >>> import pandas as pd
    >>> msa, validation = load_and_prepare_msa("protein.a3m")
    >>> df = pd.read_csv("variants.tsv", sep="\\t")
    >>> annotated = annotate_variants_with_conservation(df, msa)
"""

# Data models
from .models import (
    ConservationCategory,
    ConservationResult,
    MSAFormat,
    MSAlignment,
    MSAPosition,
    MSASequence,
    MSAValidationResult,
    SequenceType,
)

# MSA parsing
from .msa_parser import (
    detect_format,
    parse_a3m,
    parse_msa,
    validate_msa,
    write_msa,
)

# Score calculations
from .scores import (
    calculate_all_positions,
    calculate_conservation_scores,
    calculate_position_stats,
    categorize_conservation,
    get_conservation_summary,
    normalized_conservation_score,
    percent_identity,
    shannon_entropy,
)

# Position mapping
from .mapper import (
    build_query_to_msa_map,
    get_conservation_for_position,
    get_conservation_for_variant,
    parse_protein_position_from_hgvsp,
    protein_pos_to_msa_col,
)

# High-level annotation
from .annotator import (
    annotate_tsv_file,
    annotate_variants_with_conservation,
    annotate_with_gene_msas,
    batch_annotate,
    get_gene_msa_mapping,
    load_and_prepare_msa,
)


__all__ = [
    # Models
    "ConservationCategory",
    "ConservationResult",
    "MSAFormat",
    "MSAlignment",
    "MSAPosition",
    "MSASequence",
    "MSAValidationResult",
    "SequenceType",
    # Parsing
    "detect_format",
    "parse_a3m",
    "parse_msa",
    "validate_msa",
    "write_msa",
    # Scores
    "calculate_all_positions",
    "calculate_conservation_scores",
    "calculate_position_stats",
    "categorize_conservation",
    "get_conservation_summary",
    "normalized_conservation_score",
    "percent_identity",
    "shannon_entropy",
    # Mapping
    "build_query_to_msa_map",
    "get_conservation_for_position",
    "get_conservation_for_variant",
    "parse_protein_position_from_hgvsp",
    "protein_pos_to_msa_col",
    # Annotation
    "annotate_tsv_file",
    "annotate_variants_with_conservation",
    "annotate_with_gene_msas",
    "batch_annotate",
    "get_gene_msa_mapping",
    "load_and_prepare_msa",
]
