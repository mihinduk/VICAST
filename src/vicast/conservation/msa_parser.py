"""
MSA file parsing for multiple sequence alignment formats.

This module handles parsing of MSA files from various formats including:
- A3M (ColabFold/AlphaFold format)
- Stockholm
- Clustal
- FASTA

A3M format is the primary format for ColabFold MSAs and uses lowercase letters
to denote insertions relative to the query sequence. These are removed during
parsing to maintain alignment integrity.
"""

import re
from pathlib import Path
from typing import List, Optional, Tuple, Union

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .models import (
    MSAFormat,
    MSAlignment,
    MSASequence,
    MSAValidationResult,
    SequenceType,
)


def detect_format(filepath: Union[str, Path]) -> MSAFormat:
    """Detect MSA format from file extension.

    Args:
        filepath: Path to the MSA file

    Returns:
        Detected MSAFormat enum value

    Raises:
        ValueError: If format cannot be determined
    """
    path = Path(filepath)
    suffix = path.suffix.lower()

    # Handle common double extensions
    name_lower = path.name.lower()
    if name_lower.endswith(".a3m"):
        return MSAFormat.A3M
    if name_lower.endswith(".sto") or name_lower.endswith(".stockholm"):
        return MSAFormat.STOCKHOLM
    if name_lower.endswith(".aln") or name_lower.endswith(".clustal"):
        return MSAFormat.CLUSTAL
    if suffix in (".fa", ".fasta", ".faa", ".fas"):
        return MSAFormat.FASTA

    raise ValueError(f"Cannot determine MSA format from file: {filepath}")


def parse_a3m(filepath: Union[str, Path], remove_insertions: bool = True) -> MSAlignment:
    """Parse A3M format MSA file.

    A3M is the format used by ColabFold and AlphaFold. It's similar to FASTA
    but uses lowercase letters to indicate insertions relative to the query.
    These insertions are typically removed to maintain alignment.

    Args:
        filepath: Path to the A3M file
        remove_insertions: If True, remove lowercase insertion characters

    Returns:
        Parsed MSAlignment object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file is malformed
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"A3M file not found: {filepath}")

    sequences: List[MSASequence] = []
    current_id = None
    current_desc = ""
    current_seq_parts: List[str] = []

    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n\r")

            if line.startswith(">"):
                # Save previous sequence if exists
                if current_id is not None:
                    seq_str = "".join(current_seq_parts)
                    if remove_insertions:
                        # Remove lowercase letters (insertions in A3M format)
                        seq_str = re.sub(r"[a-z]", "", seq_str)
                    sequences.append(
                        MSASequence(
                            id=current_id,
                            sequence=seq_str,
                            description=current_desc,
                            is_query=(len(sequences) == 0),
                        )
                    )

                # Parse new header
                header = line[1:].strip()
                parts = header.split(None, 1)
                current_id = parts[0] if parts else ""
                current_desc = parts[1] if len(parts) > 1 else ""
                current_seq_parts = []
            else:
                # Sequence line
                current_seq_parts.append(line)

    # Don't forget the last sequence
    if current_id is not None:
        seq_str = "".join(current_seq_parts)
        if remove_insertions:
            seq_str = re.sub(r"[a-z]", "", seq_str)
        sequences.append(
            MSASequence(
                id=current_id,
                sequence=seq_str,
                description=current_desc,
                is_query=(len(sequences) == 0),
            )
        )

    if not sequences:
        raise ValueError(f"No sequences found in A3M file: {filepath}")

    # Mark first sequence as query
    if sequences:
        sequences[0].is_query = True

    msa = MSAlignment(
        sequences=sequences,
        format=MSAFormat.A3M,
        sequence_type=_detect_sequence_type(sequences[0].sequence if sequences else ""),
        query_index=0,
    )

    return msa


def parse_fasta_msa(filepath: Union[str, Path]) -> MSAlignment:
    """Parse FASTA format MSA file.

    Args:
        filepath: Path to the FASTA alignment file

    Returns:
        Parsed MSAlignment object
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {filepath}")

    try:
        alignment = AlignIO.read(str(path), "fasta")
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA alignment: {e}")

    sequences = []
    for i, record in enumerate(alignment):
        sequences.append(
            MSASequence(
                id=record.id,
                sequence=str(record.seq),
                description=record.description,
                is_query=(i == 0),
            )
        )

    if not sequences:
        raise ValueError(f"No sequences found in FASTA file: {filepath}")

    return MSAlignment(
        sequences=sequences,
        format=MSAFormat.FASTA,
        sequence_type=_detect_sequence_type(sequences[0].sequence),
        query_index=0,
    )


def parse_stockholm(filepath: Union[str, Path]) -> MSAlignment:
    """Parse Stockholm format MSA file.

    Args:
        filepath: Path to the Stockholm file

    Returns:
        Parsed MSAlignment object
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Stockholm file not found: {filepath}")

    try:
        alignment = AlignIO.read(str(path), "stockholm")
    except Exception as e:
        raise ValueError(f"Failed to parse Stockholm alignment: {e}")

    sequences = []
    for i, record in enumerate(alignment):
        sequences.append(
            MSASequence(
                id=record.id,
                sequence=str(record.seq),
                description=record.description,
                is_query=(i == 0),
            )
        )

    if not sequences:
        raise ValueError(f"No sequences found in Stockholm file: {filepath}")

    return MSAlignment(
        sequences=sequences,
        format=MSAFormat.STOCKHOLM,
        sequence_type=_detect_sequence_type(sequences[0].sequence),
        query_index=0,
    )


def parse_clustal(filepath: Union[str, Path]) -> MSAlignment:
    """Parse Clustal format MSA file.

    Args:
        filepath: Path to the Clustal file

    Returns:
        Parsed MSAlignment object
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Clustal file not found: {filepath}")

    try:
        alignment = AlignIO.read(str(path), "clustal")
    except Exception as e:
        raise ValueError(f"Failed to parse Clustal alignment: {e}")

    sequences = []
    for i, record in enumerate(alignment):
        sequences.append(
            MSASequence(
                id=record.id,
                sequence=str(record.seq),
                description=record.description,
                is_query=(i == 0),
            )
        )

    if not sequences:
        raise ValueError(f"No sequences found in Clustal file: {filepath}")

    return MSAlignment(
        sequences=sequences,
        format=MSAFormat.CLUSTAL,
        sequence_type=_detect_sequence_type(sequences[0].sequence),
        query_index=0,
    )


def parse_msa(
    filepath: Union[str, Path],
    format: Optional[MSAFormat] = None,
    remove_insertions: bool = True,
) -> MSAlignment:
    """Parse MSA file with auto-format detection.

    Args:
        filepath: Path to the MSA file
        format: Optional explicit format. If None, auto-detect from extension
        remove_insertions: For A3M format, remove lowercase insertions

    Returns:
        Parsed MSAlignment object

    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If format is unknown or file is malformed
    """
    if format is None:
        format = detect_format(filepath)

    if format == MSAFormat.A3M:
        return parse_a3m(filepath, remove_insertions=remove_insertions)
    elif format == MSAFormat.FASTA:
        return parse_fasta_msa(filepath)
    elif format == MSAFormat.STOCKHOLM:
        return parse_stockholm(filepath)
    elif format == MSAFormat.CLUSTAL:
        return parse_clustal(filepath)
    else:
        raise ValueError(f"Unsupported MSA format: {format}")


def validate_msa(msa: MSAlignment) -> MSAValidationResult:
    """Validate an MSA for use in conservation scoring.

    Checks for:
    - Minimum number of sequences
    - Consistent alignment length
    - Valid characters
    - Presence of query sequence

    Args:
        msa: MSAlignment object to validate

    Returns:
        MSAValidationResult with validation status, errors, and warnings
    """
    errors: List[str] = []
    warnings: List[str] = []
    stats = {}

    # Check for minimum sequences
    if msa.num_sequences < 2:
        errors.append(f"MSA has only {msa.num_sequences} sequence(s); need at least 2 for conservation scoring")

    # Check alignment length consistency
    if msa.sequences:
        first_len = len(msa.sequences[0])
        inconsistent = [
            (seq.id, len(seq))
            for seq in msa.sequences
            if len(seq) != first_len
        ]
        if inconsistent:
            errors.append(
                f"Inconsistent sequence lengths: expected {first_len}, "
                f"found {inconsistent[:5]}{'...' if len(inconsistent) > 5 else ''}"
            )

    # Check for query sequence
    query = msa.query_sequence
    if query is None:
        errors.append("No query sequence found (query_index out of range)")
    elif not query.is_query:
        warnings.append("Query sequence is not marked as is_query=True")

    # Check for excessive gaps in query
    if query:
        query_gap_fraction = query.sequence.count("-") / len(query.sequence) if len(query) > 0 else 0
        if query_gap_fraction > 0.5:
            warnings.append(f"Query sequence has {query_gap_fraction:.1%} gaps - may affect position mapping")

    # Check character validity
    valid_protein = set("ACDEFGHIKLMNPQRSTVWY-.*X")
    valid_nucleotide = set("ACGTURYSWKMBDHVN-.*")

    if msa.sequences:
        sample_seq = msa.sequences[0].sequence.upper()
        invalid_chars = set(sample_seq) - valid_protein - valid_nucleotide
        if invalid_chars:
            warnings.append(f"Unusual characters found: {invalid_chars}")

    # Calculate stats
    stats["num_sequences"] = msa.num_sequences
    stats["alignment_length"] = msa.alignment_length
    stats["sequence_type"] = msa.sequence_type.value
    stats["format"] = msa.format.value
    if query:
        stats["query_length"] = query.get_ungapped_length()
        stats["query_id"] = query.id

    is_valid = len(errors) == 0

    return MSAValidationResult(
        is_valid=is_valid,
        errors=errors,
        warnings=warnings,
        stats=stats,
    )


def _detect_sequence_type(sequence: str) -> SequenceType:
    """Detect whether a sequence is protein or nucleotide.

    Args:
        sequence: Sequence string (may include gaps)

    Returns:
        SequenceType.PROTEIN or SequenceType.NUCLEOTIDE
    """
    # Remove gaps and convert to uppercase
    clean_seq = sequence.replace("-", "").replace(".", "").upper()

    if not clean_seq:
        return SequenceType.PROTEIN  # Default

    # Count nucleotide vs protein characters
    nucleotide_chars = set("ACGTU")
    protein_only_chars = set("EFIPQLDHKRMWY")

    nuc_count = sum(1 for c in clean_seq if c in nucleotide_chars)
    protein_count = sum(1 for c in clean_seq if c in protein_only_chars)

    # If we see protein-only characters, it's protein
    if protein_count > 0:
        return SequenceType.PROTEIN

    # If >80% nucleotide characters, treat as nucleotide
    if len(clean_seq) > 0 and nuc_count / len(clean_seq) > 0.8:
        return SequenceType.NUCLEOTIDE

    return SequenceType.PROTEIN


def write_msa(
    msa: MSAlignment,
    filepath: Union[str, Path],
    format: Optional[MSAFormat] = None,
) -> None:
    """Write MSA to file in specified format.

    Args:
        msa: MSAlignment to write
        filepath: Output file path
        format: Output format (defaults to FASTA)
    """
    if format is None:
        format = MSAFormat.FASTA

    path = Path(filepath)

    # Convert to BioPython alignment
    records = []
    for seq in msa.sequences:
        record = SeqRecord(
            Seq(seq.sequence),
            id=seq.id,
            description=seq.description,
        )
        records.append(record)

    format_map = {
        MSAFormat.FASTA: "fasta",
        MSAFormat.STOCKHOLM: "stockholm",
        MSAFormat.CLUSTAL: "clustal",
        MSAFormat.A3M: "fasta",  # A3M is written as FASTA
    }

    with open(path, "w") as f:
        AlignIO.write([records], f, format_map.get(format, "fasta"))
