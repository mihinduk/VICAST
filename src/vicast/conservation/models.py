"""
Data models for MSA-based conservation scoring.

This module defines the core data structures used throughout the conservation
scoring pipeline, including MSA representations, position statistics, and
conservation results.
"""

from dataclasses import dataclass, field
from enum import Enum
from typing import Dict, List, Optional


class MSAFormat(Enum):
    """Supported MSA file formats."""
    A3M = "a3m"
    STOCKHOLM = "stockholm"
    CLUSTAL = "clustal"
    FASTA = "fasta"


class SequenceType(Enum):
    """Type of sequences in the alignment."""
    PROTEIN = "protein"
    NUCLEOTIDE = "nucleotide"


class ConservationCategory(Enum):
    """Conservation category classifications."""
    HIGHLY_CONSERVED = "highly_conserved"
    CONSERVED = "conserved"
    MODERATELY_CONSERVED = "moderately_conserved"
    VARIABLE = "variable"


@dataclass
class MSASequence:
    """A single sequence in the multiple sequence alignment.

    Attributes:
        id: Sequence identifier
        sequence: The aligned sequence string (including gaps)
        description: Optional description from the header
        is_query: Whether this is the query/reference sequence
    """
    id: str
    sequence: str
    description: str = ""
    is_query: bool = False

    def __len__(self) -> int:
        """Return the length of the aligned sequence (including gaps)."""
        return len(self.sequence)

    def get_ungapped_sequence(self) -> str:
        """Return the sequence without gap characters."""
        return self.sequence.replace("-", "").replace(".", "")

    def get_ungapped_length(self) -> int:
        """Return the length of the sequence without gaps."""
        return len(self.get_ungapped_sequence())


@dataclass
class MSAPosition:
    """Statistics for a single column/position in the MSA.

    Attributes:
        column_index: 0-based index of this column in the MSA
        residue_counts: Count of each residue at this position
        gap_count: Number of gaps at this position
        total_sequences: Total number of sequences
        query_residue: The residue in the query sequence (if available)
        shannon_entropy: Calculated Shannon entropy (lower = more conserved)
        percent_identity: Percentage of sequences matching consensus
        normalized_score: Conservation score normalized to 0-1 (1 = conserved)
        consensus_residue: Most common non-gap residue
        conservation_category: Categorized conservation level
    """
    column_index: int
    residue_counts: Dict[str, int] = field(default_factory=dict)
    gap_count: int = 0
    total_sequences: int = 0
    query_residue: Optional[str] = None
    shannon_entropy: float = 0.0
    percent_identity: float = 0.0
    normalized_score: float = 0.0
    consensus_residue: Optional[str] = None
    conservation_category: Optional[ConservationCategory] = None

    @property
    def gap_fraction(self) -> float:
        """Fraction of sequences with a gap at this position."""
        if self.total_sequences == 0:
            return 0.0
        return self.gap_count / self.total_sequences

    @property
    def depth(self) -> int:
        """Number of non-gap residues at this position."""
        return self.total_sequences - self.gap_count


@dataclass
class MSAlignment:
    """A multiple sequence alignment with computed position statistics.

    Attributes:
        sequences: List of aligned sequences
        positions: List of position statistics (computed after parsing)
        format: The format the MSA was parsed from
        sequence_type: Whether sequences are protein or nucleotide
        query_index: Index of the query sequence in the sequences list
        _query_to_msa_map: Mapping from ungapped query positions to MSA columns
    """
    sequences: List[MSASequence] = field(default_factory=list)
    positions: List[MSAPosition] = field(default_factory=list)
    format: MSAFormat = MSAFormat.A3M
    sequence_type: SequenceType = SequenceType.PROTEIN
    query_index: int = 0
    _query_to_msa_map: Dict[int, int] = field(default_factory=dict)

    def __len__(self) -> int:
        """Return the alignment length (number of columns)."""
        if not self.sequences:
            return 0
        return len(self.sequences[0])

    @property
    def num_sequences(self) -> int:
        """Number of sequences in the alignment."""
        return len(self.sequences)

    @property
    def alignment_length(self) -> int:
        """Length of the alignment (same as __len__)."""
        return len(self)

    @property
    def query_sequence(self) -> Optional[MSASequence]:
        """Return the query sequence."""
        if self.sequences and 0 <= self.query_index < len(self.sequences):
            return self.sequences[self.query_index]
        return None

    def get_position(self, column_index: int) -> Optional[MSAPosition]:
        """Get position statistics for a specific column.

        Args:
            column_index: 0-based column index

        Returns:
            MSAPosition for the column, or None if out of range
        """
        if 0 <= column_index < len(self.positions):
            return self.positions[column_index]
        return None

    def get_column(self, column_index: int) -> List[str]:
        """Get all residues at a specific column.

        Args:
            column_index: 0-based column index

        Returns:
            List of residues at this column from all sequences
        """
        if column_index < 0 or column_index >= len(self):
            return []
        return [seq.sequence[column_index] for seq in self.sequences]


@dataclass
class ConservationResult:
    """Conservation scores for a single variant.

    Attributes:
        protein_position: 1-based position in the protein
        msa_column: 0-based MSA column index (None if position not mapped)
        conservation_score: Normalized 0-1 score (1 = highly conserved)
        shannon_entropy: Shannon entropy (lower = more conserved)
        percent_identity: Percentage matching consensus
        consensus_aa: Most common amino acid at this position
        msa_depth: Number of sequences in the alignment
        gap_fraction: Fraction of gaps at this position
        conservation_category: Category string
        reference_aa: Amino acid in the reference/query sequence
        variant_aa: Amino acid after the variant (if available)
        found_in_msa: Whether the position was found in the MSA
        error_message: Any error message if lookup failed
    """
    protein_position: int
    msa_column: Optional[int] = None
    conservation_score: float = 0.0
    shannon_entropy: float = 0.0
    percent_identity: float = 0.0
    consensus_aa: Optional[str] = None
    msa_depth: int = 0
    gap_fraction: float = 0.0
    conservation_category: str = ""
    reference_aa: Optional[str] = None
    variant_aa: Optional[str] = None
    found_in_msa: bool = False
    error_message: Optional[str] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary for DataFrame integration."""
        return {
            "CONSERVATION_SCORE": self.conservation_score if self.found_in_msa else None,
            "SHANNON_ENTROPY": self.shannon_entropy if self.found_in_msa else None,
            "PERCENT_IDENTITY": self.percent_identity if self.found_in_msa else None,
            "CONSENSUS_AA": self.consensus_aa if self.found_in_msa else None,
            "MSA_DEPTH": self.msa_depth if self.found_in_msa else None,
            "GAP_FRACTION": self.gap_fraction if self.found_in_msa else None,
            "CONSERVATION_CATEGORY": self.conservation_category if self.found_in_msa else None,
        }


@dataclass
class MSAValidationResult:
    """Result of MSA validation.

    Attributes:
        is_valid: Whether the MSA passed validation
        errors: List of critical errors
        warnings: List of non-critical warnings
        stats: Summary statistics about the MSA
    """
    is_valid: bool
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    stats: Dict[str, any] = field(default_factory=dict)
