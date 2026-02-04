"""
Conservation score calculations for MSA positions.

This module provides functions to calculate evolutionary conservation metrics
from multiple sequence alignments, including:
- Shannon entropy (information-theoretic measure)
- Percent identity (simple frequency measure)
- Normalized conservation scores
- Conservation category classification
"""

import math
from collections import Counter
from typing import Dict, List, Optional, Tuple

from .models import (
    ConservationCategory,
    MSAlignment,
    MSAPosition,
    SequenceType,
)


# Standard amino acid background frequencies (from UniProt)
AA_BACKGROUND_FREQ = {
    "A": 0.0826, "R": 0.0553, "N": 0.0406, "D": 0.0546, "C": 0.0137,
    "Q": 0.0393, "E": 0.0674, "G": 0.0708, "H": 0.0227, "I": 0.0593,
    "L": 0.0966, "K": 0.0583, "M": 0.0242, "F": 0.0386, "P": 0.0470,
    "S": 0.0656, "T": 0.0534, "W": 0.0109, "Y": 0.0292, "V": 0.0687,
}

# Nucleotide background (equal distribution)
NT_BACKGROUND_FREQ = {
    "A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25, "U": 0.25,
}


def shannon_entropy(
    residue_counts: Dict[str, int],
    include_gaps: bool = False,
    pseudocount: float = 0.0,
) -> float:
    """Calculate Shannon entropy for a position.

    Shannon entropy H = -sum(p_i * log2(p_i)) where p_i is the frequency
    of residue i. Lower entropy means higher conservation.

    Args:
        residue_counts: Dictionary mapping residues to their counts
        include_gaps: Whether to include gap characters in the calculation
        pseudocount: Small value added to avoid log(0); typically 0 or 1/N

    Returns:
        Shannon entropy in bits. Range is [0, log2(num_residue_types)]
        For proteins: 0 (perfectly conserved) to ~4.32 (20 AAs equally frequent)
    """
    if not residue_counts:
        return 0.0

    # Filter out gaps if requested
    if not include_gaps:
        counts = {k: v for k, v in residue_counts.items() if k not in "-.*"}
    else:
        counts = residue_counts.copy()

    if not counts:
        return 0.0

    total = sum(counts.values()) + pseudocount * len(counts)
    if total == 0:
        return 0.0

    entropy = 0.0
    for count in counts.values():
        if count > 0:
            freq = (count + pseudocount) / total
            if freq > 0:
                entropy -= freq * math.log2(freq)

    return entropy


def percent_identity(
    residue_counts: Dict[str, int],
    include_gaps: bool = False,
) -> Tuple[float, Optional[str]]:
    """Calculate percent identity (frequency of most common residue).

    Args:
        residue_counts: Dictionary mapping residues to their counts
        include_gaps: Whether to include gaps in total count

    Returns:
        Tuple of (percent_identity, consensus_residue)
        percent_identity: Fraction (0-1) matching most common residue
        consensus_residue: The most common non-gap residue
    """
    if not residue_counts:
        return 0.0, None

    # Find most common non-gap residue
    non_gap_counts = {k: v for k, v in residue_counts.items() if k not in "-.*"}

    if not non_gap_counts:
        return 0.0, None

    consensus = max(non_gap_counts.keys(), key=lambda x: non_gap_counts[x])
    consensus_count = non_gap_counts[consensus]

    # Calculate total based on include_gaps flag
    if include_gaps:
        total = sum(residue_counts.values())
    else:
        total = sum(non_gap_counts.values())

    if total == 0:
        return 0.0, consensus

    return consensus_count / total, consensus


def normalized_conservation_score(
    shannon_entropy_value: float,
    sequence_type: SequenceType = SequenceType.PROTEIN,
) -> float:
    """Convert Shannon entropy to normalized conservation score.

    Maps entropy to a 0-1 scale where 1 = highly conserved (low entropy)
    and 0 = highly variable (high entropy).

    Args:
        shannon_entropy_value: Shannon entropy value
        sequence_type: Type of sequence (affects max entropy)

    Returns:
        Normalized score from 0 (variable) to 1 (conserved)
    """
    # Maximum theoretical entropy
    if sequence_type == SequenceType.PROTEIN:
        max_entropy = math.log2(20)  # ~4.32 bits for 20 amino acids
    else:
        max_entropy = math.log2(4)  # 2 bits for 4 nucleotides

    if max_entropy == 0:
        return 1.0

    # Normalize and invert so 1 = conserved
    normalized = 1.0 - (shannon_entropy_value / max_entropy)
    return max(0.0, min(1.0, normalized))  # Clamp to [0, 1]


def categorize_conservation(
    normalized_score: float,
    thresholds: Optional[Dict[str, float]] = None,
) -> ConservationCategory:
    """Categorize conservation level based on normalized score.

    Default thresholds:
    - highly_conserved: >= 0.9
    - conserved: >= 0.7
    - moderately_conserved: >= 0.5
    - variable: < 0.5

    Args:
        normalized_score: Conservation score from 0-1
        thresholds: Optional custom thresholds dict

    Returns:
        ConservationCategory enum value
    """
    if thresholds is None:
        thresholds = {
            "highly_conserved": 0.9,
            "conserved": 0.7,
            "moderately_conserved": 0.5,
        }

    if normalized_score >= thresholds.get("highly_conserved", 0.9):
        return ConservationCategory.HIGHLY_CONSERVED
    elif normalized_score >= thresholds.get("conserved", 0.7):
        return ConservationCategory.CONSERVED
    elif normalized_score >= thresholds.get("moderately_conserved", 0.5):
        return ConservationCategory.MODERATELY_CONSERVED
    else:
        return ConservationCategory.VARIABLE


def calculate_position_stats(
    column: List[str],
    sequence_type: SequenceType = SequenceType.PROTEIN,
    query_residue: Optional[str] = None,
) -> MSAPosition:
    """Calculate all conservation statistics for a single column.

    Args:
        column: List of residues at this position from all sequences
        sequence_type: Type of sequences in alignment
        query_residue: The residue in the query sequence

    Returns:
        MSAPosition with all statistics calculated
    """
    # Convert to uppercase for counting
    column_upper = [r.upper() for r in column]

    # Count residues
    counts = Counter(column_upper)

    # Separate gap counts
    gap_chars = {"-", ".", "*"}
    gap_count = sum(counts.get(g, 0) for g in gap_chars)
    residue_counts = {k: v for k, v in counts.items() if k not in gap_chars}

    total = len(column)

    # Calculate metrics
    entropy = shannon_entropy(residue_counts)
    pct_id, consensus = percent_identity(residue_counts)
    norm_score = normalized_conservation_score(entropy, sequence_type)
    category = categorize_conservation(norm_score)

    return MSAPosition(
        column_index=-1,  # Will be set by caller
        residue_counts=residue_counts,
        gap_count=gap_count,
        total_sequences=total,
        query_residue=query_residue,
        shannon_entropy=entropy,
        percent_identity=pct_id,
        normalized_score=norm_score,
        consensus_residue=consensus,
        conservation_category=category,
    )


def calculate_all_positions(
    msa: MSAlignment,
    min_depth: int = 1,
) -> List[MSAPosition]:
    """Calculate conservation statistics for all MSA positions.

    Args:
        msa: MSAlignment object
        min_depth: Minimum number of non-gap residues required

    Returns:
        List of MSAPosition objects, one per column
    """
    if not msa.sequences:
        return []

    alignment_length = len(msa.sequences[0])
    positions: List[MSAPosition] = []

    query = msa.query_sequence

    for col_idx in range(alignment_length):
        # Get column
        column = [seq.sequence[col_idx] for seq in msa.sequences]

        # Get query residue if available
        query_residue = None
        if query and col_idx < len(query.sequence):
            query_residue = query.sequence[col_idx]

        # Calculate stats
        pos = calculate_position_stats(
            column,
            sequence_type=msa.sequence_type,
            query_residue=query_residue,
        )
        pos.column_index = col_idx

        # Check depth requirement
        if pos.depth < min_depth:
            # Mark as low confidence but still include
            pass

        positions.append(pos)

    return positions


def calculate_conservation_scores(msa: MSAlignment) -> MSAlignment:
    """Calculate conservation scores for entire MSA and update in place.

    This is the main entry point for score calculation. It:
    1. Calculates position statistics for all columns
    2. Builds the query-to-MSA position mapping
    3. Updates the MSAlignment object with results

    Args:
        msa: MSAlignment object to process

    Returns:
        The same MSAlignment object with positions populated
    """
    # Calculate all positions
    msa.positions = calculate_all_positions(msa)

    # Build query-to-MSA mapping
    _build_query_to_msa_map(msa)

    return msa


def _build_query_to_msa_map(msa: MSAlignment) -> None:
    """Build mapping from ungapped query positions to MSA columns.

    This mapping allows conversion from protein positions (1-based, ungapped)
    to MSA column indices (0-based, may have gaps).

    Args:
        msa: MSAlignment object to update
    """
    query = msa.query_sequence
    if not query:
        msa._query_to_msa_map = {}
        return

    position_map: Dict[int, int] = {}
    ungapped_pos = 0  # 0-based ungapped position

    for col_idx, residue in enumerate(query.sequence):
        if residue not in "-.*":
            ungapped_pos += 1  # Now 1-based
            position_map[ungapped_pos] = col_idx

    msa._query_to_msa_map = position_map


def get_conservation_summary(msa: MSAlignment) -> Dict:
    """Get summary statistics for the MSA conservation.

    Args:
        msa: MSAlignment with calculated positions

    Returns:
        Dictionary with summary statistics
    """
    if not msa.positions:
        return {}

    entropies = [p.shannon_entropy for p in msa.positions]
    identities = [p.percent_identity for p in msa.positions]
    scores = [p.normalized_score for p in msa.positions]

    # Count categories
    category_counts = Counter(p.conservation_category for p in msa.positions)

    return {
        "num_positions": len(msa.positions),
        "num_sequences": msa.num_sequences,
        "mean_entropy": sum(entropies) / len(entropies) if entropies else 0,
        "mean_identity": sum(identities) / len(identities) if identities else 0,
        "mean_conservation": sum(scores) / len(scores) if scores else 0,
        "highly_conserved_positions": category_counts.get(ConservationCategory.HIGHLY_CONSERVED, 0),
        "conserved_positions": category_counts.get(ConservationCategory.CONSERVED, 0),
        "moderately_conserved_positions": category_counts.get(ConservationCategory.MODERATELY_CONSERVED, 0),
        "variable_positions": category_counts.get(ConservationCategory.VARIABLE, 0),
    }
