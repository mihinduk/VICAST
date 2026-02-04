"""
GFF3 validation utilities for SnpEff compatibility.

This module provides functions to validate GFF3 annotation files
for use with SnpEff variant annotation.
"""

import os
import re
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Set


def validate_gff_for_snpeff(
    gff_path: str,
    fasta_path: Optional[str] = None,
    strict: bool = False
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate a GFF3 file for SnpEff compatibility.

    Args:
        gff_path: Path to the GFF3 file to validate
        fasta_path: Optional path to reference FASTA for coordinate validation
        strict: If True, treat warnings as errors

    Returns:
        Tuple of (is_valid, errors, warnings)
        - is_valid: True if file passes validation
        - errors: List of critical error messages
        - warnings: List of warning messages
    """
    errors: List[str] = []
    warnings: List[str] = []

    # Check if file exists
    if not os.path.exists(gff_path):
        errors.append(f"GFF file not found: {gff_path}")
        return False, errors, warnings

    # Track IDs for duplicate detection
    seen_ids: Set[str] = set()

    # Track sequence lengths from FASTA if provided
    seq_lengths: Dict[str, int] = {}
    if fasta_path and os.path.exists(fasta_path):
        seq_lengths = _parse_fasta_lengths(fasta_path)

    # Parse and validate GFF
    line_number = 0
    has_version_header = False

    try:
        with open(gff_path, 'r') as f:
            for line in f:
                line_number += 1
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Check headers
                if line.startswith('##'):
                    if line.startswith('##gff-version'):
                        has_version_header = True
                        if '3' not in line:
                            warnings.append(f"Line {line_number}: Expected GFF version 3")
                    continue

                # Skip comment lines
                if line.startswith('#'):
                    continue

                # Parse feature line
                parts = line.split('\t')

                # Validate column count
                if len(parts) < 9:
                    errors.append(
                        f"Line {line_number}: Invalid GFF line - expected 9 columns, got {len(parts)}"
                    )
                    continue

                seqid, source, feature_type, start_str, end_str, score, strand, phase, attributes = parts

                # Validate coordinates
                try:
                    start = int(start_str)
                    end = int(end_str)
                except ValueError:
                    errors.append(
                        f"Line {line_number}: Invalid coordinates - start='{start_str}', end='{end_str}'"
                    )
                    continue

                if start > end:
                    errors.append(
                        f"Line {line_number}: Start ({start}) greater than end ({end})"
                    )

                if start < 1:
                    errors.append(
                        f"Line {line_number}: Start coordinate must be >= 1, got {start}"
                    )

                # Validate against FASTA if available
                if seqid in seq_lengths:
                    if end > seq_lengths[seqid]:
                        warnings.append(
                            f"Line {line_number}: End coordinate ({end}) exceeds sequence length ({seq_lengths[seqid]})"
                        )

                # Validate strand
                if strand not in ['+', '-', '.']:
                    warnings.append(
                        f"Line {line_number}: Invalid strand '{strand}', expected +, -, or ."
                    )

                # Parse and validate attributes
                attr_dict = _parse_attributes(attributes)

                # Check for required ID attribute on CDS/mRNA/gene features
                if feature_type in ['CDS', 'mRNA', 'gene', 'exon']:
                    if 'ID' not in attr_dict:
                        errors.append(
                            f"Line {line_number}: Missing ID attribute for {feature_type} feature"
                        )
                    else:
                        feature_id = attr_dict['ID']
                        if feature_id in seen_ids:
                            errors.append(
                                f"Line {line_number}: Duplicate ID '{feature_id}'"
                            )
                        seen_ids.add(feature_id)

                # Warn about missing gene attribute
                if feature_type in ['CDS', 'mRNA'] and 'gene' not in attr_dict:
                    warnings.append(
                        f"Line {line_number}: {feature_type} feature missing 'gene' attribute"
                    )

                # Validate phase for CDS
                if feature_type == 'CDS' and phase not in ['0', '1', '2']:
                    warnings.append(
                        f"Line {line_number}: CDS phase should be 0, 1, or 2, got '{phase}'"
                    )

    except Exception as e:
        errors.append(f"Error reading GFF file: {str(e)}")

    # Check for version header
    if not has_version_header:
        warnings.append("Missing ##gff-version header")

    # Determine validity
    is_valid = len(errors) == 0
    if strict and len(warnings) > 0:
        is_valid = False

    return is_valid, errors, warnings


def _parse_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GFF3 attribute string into dictionary.

    Args:
        attr_string: Semicolon-separated attribute string

    Returns:
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    for attr in attr_string.split(';'):
        attr = attr.strip()
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key] = value
    return attributes


def _parse_fasta_lengths(fasta_path: str) -> Dict[str, int]:
    """
    Parse FASTA file to get sequence lengths.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary mapping sequence IDs to lengths
    """
    lengths = {}
    current_id = None
    current_length = 0

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    lengths[current_id] = current_length
                # Extract ID (first word after >)
                current_id = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

        # Don't forget the last sequence
        if current_id is not None:
            lengths[current_id] = current_length

    return lengths


def _parse_fasta_sequences(fasta_path: str) -> Dict[str, str]:
    """
    Parse FASTA file to get full sequences.

    Args:
        fasta_path: Path to FASTA file

    Returns:
        Dictionary mapping sequence IDs to sequences (uppercase)
    """
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq).upper()
                # Extract ID (first word after >)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq).upper()

    return sequences


def validate_vcf(vcf_path: str) -> Tuple[bool, List[str], List[str]]:
    """
    Basic VCF file validation.

    Args:
        vcf_path: Path to VCF file

    Returns:
        Tuple of (is_valid, errors, warnings)
    """
    errors = []
    warnings = []

    if not os.path.exists(vcf_path):
        errors.append(f"VCF file not found: {vcf_path}")
        return False, errors, warnings

    has_header = False
    line_number = 0

    with open(vcf_path, 'r') as f:
        for line in f:
            line_number += 1
            line = line.strip()

            if not line:
                continue

            if line.startswith('##'):
                if line.startswith('##fileformat=VCF'):
                    has_header = True
                continue

            if line.startswith('#CHROM'):
                continue

            # Data line
            parts = line.split('\t')
            if len(parts) < 8:
                errors.append(f"Line {line_number}: Invalid VCF line - expected at least 8 columns")
                continue

            # Validate POS
            try:
                pos = int(parts[1])
                if pos < 1:
                    errors.append(f"Line {line_number}: Position must be >= 1")
            except ValueError:
                errors.append(f"Line {line_number}: Invalid position '{parts[1]}'")

            # Validate QUAL
            qual = parts[5]
            if qual != '.':
                try:
                    float(qual)
                except ValueError:
                    warnings.append(f"Line {line_number}: Invalid QUAL '{qual}'")

    if not has_header:
        warnings.append("Missing ##fileformat header")

    is_valid = len(errors) == 0
    return is_valid, errors, warnings


def validate_vcf_ref_bases(
    vcf_path: str,
    fasta_path: str,
    max_errors: int = 100
) -> Tuple[bool, List[str], List[str]]:
    """
    Validate that VCF REF bases match the reference genome.

    VCF REF bases must exactly match the bases at those positions in the
    reference FASTA. This catches common errors like using placeholder REF
    values or coordinate mismatches.

    Args:
        vcf_path: Path to VCF file
        fasta_path: Path to reference FASTA file
        max_errors: Stop after this many errors (to avoid overwhelming output)

    Returns:
        Tuple of (is_valid, errors, warnings)
        - is_valid: True if all REF bases match the reference
        - errors: List of mismatch error messages
        - warnings: List of warning messages (e.g., unknown chromosome)
    """
    errors: List[str] = []
    warnings: List[str] = []

    # Check files exist
    if not os.path.exists(vcf_path):
        errors.append(f"VCF file not found: {vcf_path}")
        return False, errors, warnings

    if not os.path.exists(fasta_path):
        errors.append(f"Reference FASTA not found: {fasta_path}")
        return False, errors, warnings

    # Parse reference sequences
    try:
        sequences = _parse_fasta_sequences(fasta_path)
    except Exception as e:
        errors.append(f"Error parsing FASTA: {str(e)}")
        return False, errors, warnings

    if not sequences:
        errors.append("No sequences found in reference FASTA")
        return False, errors, warnings

    # Track statistics
    total_variants = 0
    matched_variants = 0
    unknown_chroms: Set[str] = set()

    # Parse VCF and check REF bases
    line_number = 0
    error_count = 0

    with open(vcf_path, 'r') as f:
        for line in f:
            line_number += 1
            line = line.strip()

            # Skip empty lines and headers
            if not line or line.startswith('#'):
                continue

            # Parse variant line
            parts = line.split('\t')
            if len(parts) < 5:
                continue

            chrom = parts[0]
            try:
                pos = int(parts[1])
            except ValueError:
                continue

            ref = parts[3].upper()
            total_variants += 1

            # Check if chromosome is in reference
            if chrom not in sequences:
                if chrom not in unknown_chroms:
                    unknown_chroms.add(chrom)
                    warnings.append(
                        f"Chromosome '{chrom}' not found in reference FASTA"
                    )
                continue

            # Get the actual reference bases
            seq = sequences[chrom]
            seq_len = len(seq)

            # VCF is 1-based, Python is 0-based
            start_idx = pos - 1
            end_idx = start_idx + len(ref)

            # Check bounds
            if start_idx < 0:
                errors.append(
                    f"Line {line_number}: Position {pos} is invalid (< 1)"
                )
                error_count += 1
                if error_count >= max_errors:
                    warnings.append(f"Stopped after {max_errors} errors")
                    break
                continue

            if end_idx > seq_len:
                errors.append(
                    f"Line {line_number}: REF extends beyond sequence end "
                    f"(pos={pos}, ref_len={len(ref)}, seq_len={seq_len})"
                )
                error_count += 1
                if error_count >= max_errors:
                    warnings.append(f"Stopped after {max_errors} errors")
                    break
                continue

            # Get actual bases from reference
            actual_ref = seq[start_idx:end_idx]

            # Compare
            if ref != actual_ref:
                errors.append(
                    f"Line {line_number}: REF mismatch at {chrom}:{pos} - "
                    f"VCF has '{ref}', reference has '{actual_ref}'"
                )
                error_count += 1
                if error_count >= max_errors:
                    warnings.append(f"Stopped after {max_errors} errors")
                    break
            else:
                matched_variants += 1

    # Summary warning if there were issues
    if total_variants > 0 and matched_variants < total_variants:
        checked = total_variants - len(unknown_chroms)
        if checked > 0:
            match_rate = (matched_variants / checked) * 100
            if match_rate < 100:
                warnings.append(
                    f"REF base validation: {matched_variants}/{checked} variants matched ({match_rate:.1f}%)"
                )

    is_valid = len(errors) == 0
    return is_valid, errors, warnings
