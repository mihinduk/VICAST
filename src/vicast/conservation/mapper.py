"""
Position mapping between HGVSp notation and MSA columns.

This module provides functions to:
- Parse protein positions from HGVSp notation (e.g., "p.Ala123Gly")
- Map protein positions to MSA column indices
- Look up conservation scores for variants
"""

import re
from typing import Dict, Optional, Tuple

from .models import ConservationResult, MSAlignment, MSAPosition


# Amino acid 3-letter to 1-letter mapping
AA_3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Ter": "*", "Stop": "*", "X": "X",
}

AA_1TO3 = {v: k for k, v in AA_3TO1.items()}


def parse_protein_position_from_hgvsp(hgvsp: str) -> Tuple[Optional[int], Optional[str], Optional[str]]:
    """Extract protein position and amino acids from HGVSp notation.

    Handles various HGVSp formats:
    - p.Ala123Gly (missense)
    - p.Gly45Ter (nonsense)
    - p.Leu100del (deletion)
    - p.Ser50_Gly52del (range deletion)
    - p.Ala10fs (frameshift)
    - p.Met1? (initiation codon variant)
    - p.= (synonymous)

    Args:
        hgvsp: HGVSp notation string (e.g., "p.Ala123Gly")

    Returns:
        Tuple of (position, reference_aa, variant_aa)
        position: 1-based protein position (None if parsing failed)
        reference_aa: Reference amino acid in 1-letter code
        variant_aa: Variant amino acid in 1-letter code (None for some variants)
    """
    if not hgvsp or not isinstance(hgvsp, str):
        return None, None, None

    # Handle empty or non-standard values
    hgvsp = hgvsp.strip()
    if hgvsp in ("", ".", "-", "NA", "N/A", "p.="):
        return None, None, None

    # Remove protein prefix if present (handle both "p." and transcript prefix)
    if ":" in hgvsp:
        hgvsp = hgvsp.split(":")[-1]  # Take part after colon

    if hgvsp.startswith("p."):
        hgvsp = hgvsp[2:]

    # Pattern for standard missense: Ref123Alt
    # e.g., Ala123Gly, Met1Val, Ter456Leu
    missense_pattern = r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|[A-Z]|\*|\?)$"
    match = re.match(missense_pattern, hgvsp)
    if match:
        ref_aa_3 = match.group(1)
        position = int(match.group(2))
        alt = match.group(3)

        ref_aa = AA_3TO1.get(ref_aa_3, ref_aa_3[0] if ref_aa_3 else None)

        if alt in ("*", "?"):
            alt_aa = alt
        elif len(alt) == 1:
            alt_aa = alt
        else:
            alt_aa = AA_3TO1.get(alt, alt[0] if alt else None)

        return position, ref_aa, alt_aa

    # Pattern for frameshift: Ref123fs
    fs_pattern = r"^([A-Z][a-z]{2})(\d+)fs"
    match = re.match(fs_pattern, hgvsp)
    if match:
        ref_aa_3 = match.group(1)
        position = int(match.group(2))
        ref_aa = AA_3TO1.get(ref_aa_3, ref_aa_3[0] if ref_aa_3 else None)
        return position, ref_aa, "fs"

    # Pattern for deletion: Ref123del or Ref123_Ref456del
    del_pattern = r"^([A-Z][a-z]{2})(\d+)(?:_[A-Z][a-z]{2}\d+)?del"
    match = re.match(del_pattern, hgvsp)
    if match:
        ref_aa_3 = match.group(1)
        position = int(match.group(2))
        ref_aa = AA_3TO1.get(ref_aa_3, ref_aa_3[0] if ref_aa_3 else None)
        return position, ref_aa, "del"

    # Pattern for insertion: Ref123_Ref124insAlt
    ins_pattern = r"^([A-Z][a-z]{2})(\d+)_[A-Z][a-z]{2}\d+ins"
    match = re.match(ins_pattern, hgvsp)
    if match:
        ref_aa_3 = match.group(1)
        position = int(match.group(2))
        ref_aa = AA_3TO1.get(ref_aa_3, ref_aa_3[0] if ref_aa_3 else None)
        return position, ref_aa, "ins"

    # Pattern for duplication: Ref123dup
    dup_pattern = r"^([A-Z][a-z]{2})(\d+)(?:_[A-Z][a-z]{2}\d+)?dup"
    match = re.match(dup_pattern, hgvsp)
    if match:
        ref_aa_3 = match.group(1)
        position = int(match.group(2))
        ref_aa = AA_3TO1.get(ref_aa_3, ref_aa_3[0] if ref_aa_3 else None)
        return position, ref_aa, "dup"

    # Pattern for extension: Ter123ExtTer456 or *123Ext*456
    ext_pattern = r"^(?:Ter|\*)(\d+)"
    match = re.match(ext_pattern, hgvsp)
    if match:
        position = int(match.group(1))
        return position, "*", "ext"

    # Pattern for just position extraction (fallback)
    # Try to find any 3-letter AA code followed by number
    fallback_pattern = r"([A-Z][a-z]{2})(\d+)"
    match = re.search(fallback_pattern, hgvsp)
    if match:
        ref_aa_3 = match.group(1)
        position = int(match.group(2))
        ref_aa = AA_3TO1.get(ref_aa_3, ref_aa_3[0] if ref_aa_3 else None)
        return position, ref_aa, None

    return None, None, None


def build_query_to_msa_map(msa: MSAlignment) -> Dict[int, int]:
    """Build mapping from ungapped query positions to MSA columns.

    Args:
        msa: MSAlignment with query sequence

    Returns:
        Dict mapping 1-based protein positions to 0-based MSA columns
    """
    query = msa.query_sequence
    if not query:
        return {}

    position_map: Dict[int, int] = {}
    ungapped_pos = 0  # Will become 1-based

    for col_idx, residue in enumerate(query.sequence):
        if residue not in "-.*":
            ungapped_pos += 1  # 1-based position
            position_map[ungapped_pos] = col_idx

    return position_map


def protein_pos_to_msa_col(
    position: int,
    msa: MSAlignment,
) -> Optional[int]:
    """Convert 1-based protein position to 0-based MSA column.

    Args:
        position: 1-based protein position
        msa: MSAlignment with position mapping

    Returns:
        0-based MSA column index, or None if position not found
    """
    if not msa._query_to_msa_map:
        # Build map if not already built
        msa._query_to_msa_map = build_query_to_msa_map(msa)

    return msa._query_to_msa_map.get(position)


def get_conservation_for_variant(
    hgvsp: str,
    msa: MSAlignment,
) -> ConservationResult:
    """Get conservation scores for a variant.

    Args:
        hgvsp: HGVSp notation for the variant
        msa: MSAlignment with calculated positions

    Returns:
        ConservationResult with all conservation metrics
    """
    # Parse the HGVSp
    position, ref_aa, var_aa = parse_protein_position_from_hgvsp(hgvsp)

    if position is None:
        return ConservationResult(
            protein_position=0,
            found_in_msa=False,
            error_message=f"Could not parse protein position from: {hgvsp}",
        )

    # Get MSA column
    msa_col = protein_pos_to_msa_col(position, msa)

    if msa_col is None:
        return ConservationResult(
            protein_position=position,
            found_in_msa=False,
            reference_aa=ref_aa,
            variant_aa=var_aa,
            error_message=f"Position {position} not found in MSA (query length: {msa.query_sequence.get_ungapped_length() if msa.query_sequence else 0})",
        )

    # Get position stats
    pos_stats = msa.get_position(msa_col)

    if pos_stats is None:
        return ConservationResult(
            protein_position=position,
            msa_column=msa_col,
            found_in_msa=False,
            reference_aa=ref_aa,
            variant_aa=var_aa,
            error_message=f"No statistics for MSA column {msa_col}",
        )

    # Build result
    return ConservationResult(
        protein_position=position,
        msa_column=msa_col,
        conservation_score=pos_stats.normalized_score,
        shannon_entropy=pos_stats.shannon_entropy,
        percent_identity=pos_stats.percent_identity,
        consensus_aa=pos_stats.consensus_residue,
        msa_depth=msa.num_sequences,
        gap_fraction=pos_stats.gap_fraction,
        conservation_category=pos_stats.conservation_category.value if pos_stats.conservation_category else "",
        reference_aa=ref_aa,
        variant_aa=var_aa,
        found_in_msa=True,
    )


def get_conservation_for_position(
    protein_position: int,
    msa: MSAlignment,
) -> ConservationResult:
    """Get conservation scores for a specific protein position.

    Args:
        protein_position: 1-based protein position
        msa: MSAlignment with calculated positions

    Returns:
        ConservationResult with all conservation metrics
    """
    # Get MSA column
    msa_col = protein_pos_to_msa_col(protein_position, msa)

    if msa_col is None:
        return ConservationResult(
            protein_position=protein_position,
            found_in_msa=False,
            error_message=f"Position {protein_position} not found in MSA",
        )

    # Get position stats
    pos_stats = msa.get_position(msa_col)

    if pos_stats is None:
        return ConservationResult(
            protein_position=protein_position,
            msa_column=msa_col,
            found_in_msa=False,
            error_message=f"No statistics for MSA column {msa_col}",
        )

    # Build result
    return ConservationResult(
        protein_position=protein_position,
        msa_column=msa_col,
        conservation_score=pos_stats.normalized_score,
        shannon_entropy=pos_stats.shannon_entropy,
        percent_identity=pos_stats.percent_identity,
        consensus_aa=pos_stats.consensus_residue,
        msa_depth=msa.num_sequences,
        gap_fraction=pos_stats.gap_fraction,
        conservation_category=pos_stats.conservation_category.value if pos_stats.conservation_category else "",
        reference_aa=pos_stats.query_residue,
        found_in_msa=True,
    )
