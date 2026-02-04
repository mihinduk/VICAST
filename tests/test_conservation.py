"""
Tests for VICAST conservation scoring module.

Tests cover:
- MSA parsing (A3M, FASTA formats)
- Shannon entropy and percent identity calculations
- HGVSp position parsing
- Position mapping
- TSV annotation
"""

import math
import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import pytest

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vicast.conservation import (
    ConservationCategory,
    ConservationResult,
    MSAFormat,
    MSAlignment,
    MSAPosition,
    MSASequence,
    SequenceType,
    annotate_variants_with_conservation,
    build_query_to_msa_map,
    calculate_conservation_scores,
    calculate_position_stats,
    categorize_conservation,
    detect_format,
    get_conservation_for_position,
    get_conservation_for_variant,
    load_and_prepare_msa,
    normalized_conservation_score,
    parse_a3m,
    parse_msa,
    parse_protein_position_from_hgvsp,
    percent_identity,
    protein_pos_to_msa_col,
    shannon_entropy,
    validate_msa,
)

# Test fixtures directory
FIXTURES_DIR = Path(__file__).parent.parent / "examples" / "data"


class TestMSAModels:
    """Tests for MSA data models."""

    def test_msa_sequence_basic(self):
        """Test MSASequence creation and properties."""
        seq = MSASequence(
            id="test_seq",
            sequence="MRCVGI--GNRDF",
            description="Test sequence",
            is_query=True,
        )

        assert seq.id == "test_seq"
        assert len(seq) == 13
        assert seq.get_ungapped_sequence() == "MRCVGIGNRDF"
        assert seq.get_ungapped_length() == 11
        assert seq.is_query is True

    def test_msa_alignment_basic(self):
        """Test MSAlignment creation and properties."""
        seqs = [
            MSASequence(id="query", sequence="ACDEF", is_query=True),
            MSASequence(id="seq2", sequence="ACDEF"),
            MSASequence(id="seq3", sequence="ACDGF"),
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)

        assert msa.num_sequences == 3
        assert msa.alignment_length == 5
        assert msa.query_sequence.id == "query"
        assert msa.get_column(0) == ["A", "A", "A"]

    def test_conservation_result_to_dict(self):
        """Test ConservationResult conversion to dict."""
        result = ConservationResult(
            protein_position=10,
            msa_column=9,
            conservation_score=0.85,
            shannon_entropy=0.5,
            percent_identity=0.9,
            consensus_aa="A",
            msa_depth=100,
            gap_fraction=0.05,
            conservation_category="conserved",
            found_in_msa=True,
        )

        d = result.to_dict()
        assert d["CONSERVATION_SCORE"] == 0.85
        assert d["SHANNON_ENTROPY"] == 0.5
        assert d["MSA_DEPTH"] == 100

    def test_conservation_result_not_found(self):
        """Test ConservationResult when not found."""
        result = ConservationResult(
            protein_position=999,
            found_in_msa=False,
            error_message="Not found",
        )

        d = result.to_dict()
        assert d["CONSERVATION_SCORE"] is None
        assert d["MSA_DEPTH"] is None


class TestMSAParsing:
    """Tests for MSA file parsing."""

    def test_detect_format_a3m(self, tmp_path):
        """Test format detection for A3M files."""
        a3m_file = tmp_path / "test.a3m"
        a3m_file.touch()
        assert detect_format(a3m_file) == MSAFormat.A3M

    def test_detect_format_fasta(self, tmp_path):
        """Test format detection for FASTA files."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.touch()
        assert detect_format(fasta_file) == MSAFormat.FASTA

        fa_file = tmp_path / "test.fa"
        fa_file.touch()
        assert detect_format(fa_file) == MSAFormat.FASTA

    def test_detect_format_unknown(self, tmp_path):
        """Test format detection raises for unknown format."""
        unknown_file = tmp_path / "test.xyz"
        unknown_file.touch()
        with pytest.raises(ValueError, match="Cannot determine"):
            detect_format(unknown_file)

    def test_parse_a3m_basic(self, tmp_path):
        """Test basic A3M parsing."""
        a3m_content = """>query_seq Query protein
MRCVGIGNRDF
>seq2 Second sequence
MRCVGIGNRDFinsertionLETTERS
>seq3 Third sequence
MRCVGIaNRDF
"""
        a3m_file = tmp_path / "test.a3m"
        a3m_file.write_text(a3m_content)

        msa = parse_a3m(a3m_file)

        assert msa.num_sequences == 3
        assert msa.sequences[0].id == "query_seq"
        assert msa.sequences[0].is_query is True
        # Lowercase should be removed
        assert "insertion" not in msa.sequences[1].sequence
        assert msa.sequences[1].sequence == "MRCVGIGNRDFLETTERS"

    def test_parse_a3m_preserves_gaps(self, tmp_path):
        """Test that A3M parsing preserves gap characters."""
        a3m_content = """>query
MRCVGI--GNRDF
>seq2
MRCVGIADGNRDF
"""
        a3m_file = tmp_path / "test.a3m"
        a3m_file.write_text(a3m_content)

        msa = parse_a3m(a3m_file)

        assert "--" in msa.sequences[0].sequence
        assert msa.alignment_length == 13

    def test_parse_msa_auto_detect(self):
        """Test auto-format detection in parse_msa."""
        test_file = FIXTURES_DIR / "test_virus_msa.a3m"
        if test_file.exists():
            msa = parse_msa(test_file)
            assert msa.num_sequences > 0
            assert msa.format == MSAFormat.A3M

    def test_parse_a3m_file_not_found(self):
        """Test error handling for missing file."""
        with pytest.raises(FileNotFoundError):
            parse_a3m("/nonexistent/file.a3m")

    def test_validate_msa_valid(self):
        """Test validation of valid MSA."""
        seqs = [
            MSASequence(id="query", sequence="ACDEF", is_query=True),
            MSASequence(id="seq2", sequence="ACDEF"),
            MSASequence(id="seq3", sequence="ACDGF"),
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)

        result = validate_msa(msa)

        assert result.is_valid is True
        assert len(result.errors) == 0
        assert result.stats["num_sequences"] == 3

    def test_validate_msa_too_few_sequences(self):
        """Test validation fails with single sequence."""
        seqs = [MSASequence(id="query", sequence="ACDEF", is_query=True)]
        msa = MSAlignment(sequences=seqs)

        result = validate_msa(msa)

        assert result.is_valid is False
        assert any("at least 2" in e for e in result.errors)

    def test_validate_msa_inconsistent_length(self):
        """Test validation catches inconsistent lengths."""
        seqs = [
            MSASequence(id="query", sequence="ACDEF", is_query=True),
            MSASequence(id="seq2", sequence="ACD"),  # Different length
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)

        result = validate_msa(msa)

        assert result.is_valid is False
        assert any("Inconsistent" in e for e in result.errors)


class TestConservationScores:
    """Tests for conservation score calculations."""

    def test_shannon_entropy_conserved(self):
        """Test entropy for perfectly conserved position."""
        counts = {"A": 100}
        entropy = shannon_entropy(counts)
        assert entropy == 0.0  # Perfect conservation

    def test_shannon_entropy_variable(self):
        """Test entropy for variable position."""
        # Equal distribution of 4 residues = 2 bits
        counts = {"A": 25, "C": 25, "D": 25, "E": 25}
        entropy = shannon_entropy(counts)
        assert abs(entropy - 2.0) < 0.001

    def test_shannon_entropy_with_gaps(self):
        """Test entropy calculation ignores gaps by default."""
        counts = {"A": 50, "-": 50}
        entropy = shannon_entropy(counts, include_gaps=False)
        assert entropy == 0.0  # Only A counted

    def test_percent_identity_full(self):
        """Test percent identity with perfect conservation."""
        counts = {"A": 100}
        pct, consensus = percent_identity(counts)
        assert pct == 1.0
        assert consensus == "A"

    def test_percent_identity_partial(self):
        """Test percent identity with partial conservation."""
        counts = {"A": 75, "G": 25}
        pct, consensus = percent_identity(counts)
        assert pct == 0.75
        assert consensus == "A"

    def test_normalized_score_conserved(self):
        """Test normalized score for conserved position."""
        score = normalized_conservation_score(0.0, SequenceType.PROTEIN)
        assert score == 1.0

    def test_normalized_score_variable(self):
        """Test normalized score for variable position."""
        max_entropy = math.log2(20)  # Max for 20 AAs
        score = normalized_conservation_score(max_entropy, SequenceType.PROTEIN)
        assert abs(score) < 0.001

    def test_categorize_conservation(self):
        """Test conservation categorization."""
        assert categorize_conservation(0.95) == ConservationCategory.HIGHLY_CONSERVED
        assert categorize_conservation(0.8) == ConservationCategory.CONSERVED
        assert categorize_conservation(0.6) == ConservationCategory.MODERATELY_CONSERVED
        assert categorize_conservation(0.3) == ConservationCategory.VARIABLE

    def test_calculate_position_stats(self):
        """Test full position statistics calculation."""
        column = ["A", "A", "A", "G", "-"]
        stats = calculate_position_stats(column, SequenceType.PROTEIN)

        assert stats.residue_counts["A"] == 3
        assert stats.residue_counts["G"] == 1
        assert stats.gap_count == 1
        assert stats.total_sequences == 5
        assert stats.consensus_residue == "A"
        assert stats.percent_identity == 0.75  # 3/4 non-gap


class TestPositionMapping:
    """Tests for HGVSp parsing and position mapping."""

    def test_parse_hgvsp_missense(self):
        """Test parsing standard missense variant."""
        pos, ref, alt = parse_protein_position_from_hgvsp("p.Ala123Gly")
        assert pos == 123
        assert ref == "A"
        assert alt == "G"

    def test_parse_hgvsp_with_transcript(self):
        """Test parsing HGVSp with transcript prefix."""
        pos, ref, alt = parse_protein_position_from_hgvsp("YP_001234:p.Met1Val")
        assert pos == 1
        assert ref == "M"
        assert alt == "V"

    def test_parse_hgvsp_nonsense(self):
        """Test parsing nonsense (stop) variant."""
        pos, ref, alt = parse_protein_position_from_hgvsp("p.Gln45Ter")
        assert pos == 45
        assert ref == "Q"
        assert alt == "*"

    def test_parse_hgvsp_frameshift(self):
        """Test parsing frameshift variant."""
        pos, ref, alt = parse_protein_position_from_hgvsp("p.Leu100fs")
        assert pos == 100
        assert ref == "L"
        assert alt == "fs"

    def test_parse_hgvsp_deletion(self):
        """Test parsing deletion variant."""
        pos, ref, alt = parse_protein_position_from_hgvsp("p.Ser50del")
        assert pos == 50
        assert ref == "S"
        assert alt == "del"

    def test_parse_hgvsp_range_deletion(self):
        """Test parsing range deletion."""
        pos, ref, alt = parse_protein_position_from_hgvsp("p.Ser50_Gly52del")
        assert pos == 50
        assert ref == "S"
        assert alt == "del"

    def test_parse_hgvsp_invalid(self):
        """Test parsing invalid HGVSp returns None."""
        pos, ref, alt = parse_protein_position_from_hgvsp("")
        assert pos is None

        pos, ref, alt = parse_protein_position_from_hgvsp(".")
        assert pos is None

        pos, ref, alt = parse_protein_position_from_hgvsp(None)
        assert pos is None

    def test_build_query_to_msa_map(self):
        """Test query-to-MSA position mapping."""
        seqs = [
            MSASequence(id="query", sequence="AC-DEF", is_query=True),
            MSASequence(id="seq2", sequence="ACGDEF"),
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)

        mapping = build_query_to_msa_map(msa)

        # 1-based protein positions to 0-based MSA columns
        assert mapping[1] == 0  # A at column 0
        assert mapping[2] == 1  # C at column 1
        assert mapping[3] == 3  # D at column 3 (gap at 2)
        assert mapping[4] == 4  # E at column 4
        assert mapping[5] == 5  # F at column 5

    def test_protein_pos_to_msa_col(self):
        """Test position conversion."""
        seqs = [
            MSASequence(id="query", sequence="AC-DEF", is_query=True),
            MSASequence(id="seq2", sequence="ACGDEF"),
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)
        msa._query_to_msa_map = build_query_to_msa_map(msa)

        assert protein_pos_to_msa_col(1, msa) == 0
        assert protein_pos_to_msa_col(3, msa) == 3
        assert protein_pos_to_msa_col(99, msa) is None  # Out of range


class TestConservationLookup:
    """Tests for conservation score lookup."""

    @pytest.fixture
    def prepared_msa(self):
        """Create a prepared MSA with scores."""
        seqs = [
            MSASequence(id="query", sequence="ACDEF", is_query=True),
            MSASequence(id="seq2", sequence="ACDEF"),
            MSASequence(id="seq3", sequence="ACDGF"),
            MSASequence(id="seq4", sequence="ACDGF"),
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)
        return calculate_conservation_scores(msa)

    def test_get_conservation_for_variant(self, prepared_msa):
        """Test getting conservation for a variant."""
        # Position 4 should have E in 2 seqs, G in 2 seqs
        result = get_conservation_for_variant("p.Glu4Gly", prepared_msa)

        assert result.found_in_msa is True
        assert result.protein_position == 4
        assert result.msa_column == 3
        assert result.msa_depth == 4
        assert result.percent_identity == 0.5  # 2/4 for most common

    def test_get_conservation_for_conserved_position(self, prepared_msa):
        """Test conservation for fully conserved position."""
        result = get_conservation_for_variant("p.Ala1Val", prepared_msa)

        assert result.found_in_msa is True
        assert result.conservation_score == 1.0
        assert result.percent_identity == 1.0
        assert result.consensus_aa == "A"

    def test_get_conservation_invalid_hgvsp(self, prepared_msa):
        """Test handling of invalid HGVSp."""
        result = get_conservation_for_variant("invalid", prepared_msa)

        assert result.found_in_msa is False
        assert result.error_message is not None

    def test_get_conservation_position_out_of_range(self, prepared_msa):
        """Test handling of position beyond MSA."""
        result = get_conservation_for_variant("p.Ala999Gly", prepared_msa)

        assert result.found_in_msa is False
        assert "not found" in result.error_message.lower()


class TestAnnotation:
    """Tests for TSV annotation functions."""

    @pytest.fixture
    def sample_variants_df(self):
        """Create sample variants DataFrame."""
        return pd.DataFrame({
            "CHROM": ["NC_001", "NC_001", "NC_001"],
            "POS": [100, 200, 300],
            "REF": ["A", "G", "C"],
            "ALT": ["G", "T", "A"],
            "HGVSp": ["p.Ala1Val", "p.Glu4Gly", "p.Phe5Leu"],
            "GENE_NAME": ["NS1", "NS1", "NS1"],
        })

    @pytest.fixture
    def prepared_msa(self):
        """Create a prepared MSA."""
        seqs = [
            MSASequence(id="query", sequence="ACDEF", is_query=True),
            MSASequence(id="seq2", sequence="ACDEF"),
            MSASequence(id="seq3", sequence="ACDGF"),
        ]
        msa = MSAlignment(sequences=seqs, query_index=0)
        return calculate_conservation_scores(msa)

    def test_annotate_variants_adds_columns(self, sample_variants_df, prepared_msa):
        """Test that annotation adds expected columns."""
        annotated = annotate_variants_with_conservation(
            sample_variants_df, prepared_msa
        )

        expected_cols = [
            "CONSERVATION_SCORE",
            "SHANNON_ENTROPY",
            "PERCENT_IDENTITY",
            "CONSENSUS_AA",
            "MSA_DEPTH",
            "GAP_FRACTION",
            "CONSERVATION_CATEGORY",
        ]

        for col in expected_cols:
            assert col in annotated.columns

    def test_annotate_variants_preserves_original(self, sample_variants_df, prepared_msa):
        """Test that annotation preserves original columns."""
        annotated = annotate_variants_with_conservation(
            sample_variants_df, prepared_msa
        )

        assert "CHROM" in annotated.columns
        assert "POS" in annotated.columns
        assert "HGVSp" in annotated.columns

    def test_annotate_handles_missing_hgvsp(self, prepared_msa):
        """Test annotation handles missing HGVSp values."""
        df = pd.DataFrame({
            "CHROM": ["NC_001", "NC_001"],
            "POS": [100, 200],
            "HGVSp": ["p.Ala1Val", None],
        })

        annotated = annotate_variants_with_conservation(df, prepared_msa)

        # First row should have scores
        assert annotated.iloc[0]["CONSERVATION_SCORE"] is not None
        assert not pd.isna(annotated.iloc[0]["CONSERVATION_SCORE"])
        # Second row should be None/NaN
        assert pd.isna(annotated.iloc[1]["CONSERVATION_SCORE"])


class TestIntegration:
    """Integration tests using real test files."""

    def test_load_and_prepare_msa(self):
        """Test loading and preparing a real MSA file."""
        test_file = FIXTURES_DIR / "test_virus_msa.a3m"
        if not test_file.exists():
            pytest.skip("Test MSA file not found")

        msa, validation = load_and_prepare_msa(test_file)

        assert validation.is_valid
        assert msa.num_sequences > 0
        assert len(msa.positions) > 0
        assert len(msa._query_to_msa_map) > 0

    def test_full_annotation_pipeline(self, tmp_path):
        """Test full annotation pipeline with temp files."""
        # Create test MSA
        a3m_content = """>query_protein NS1
MRCVGIGNRDF
>seq2 Related protein
MRCVGIGNRDF
>seq3 Related protein with variant
MRCVGMGNRDF
"""
        msa_file = tmp_path / "test.a3m"
        msa_file.write_text(a3m_content)

        # Create test variants
        variants_content = """CHROM\tPOS\tREF\tALT\tHGVSp\tGENE_NAME
NC_001\t100\tA\tG\tp.Ile6Met\tNS1
NC_001\t200\tG\tT\tp.Arg10Trp\tNS1
"""
        variants_file = tmp_path / "variants.tsv"
        variants_file.write_text(variants_content)

        # Run annotation
        msa, validation = load_and_prepare_msa(msa_file)
        df = pd.read_csv(variants_file, sep="\t")
        annotated = annotate_variants_with_conservation(df, msa)

        # Check results
        assert "CONSERVATION_SCORE" in annotated.columns
        assert len(annotated) == 2

        # Position 6 (I) - 2/3 have I, 1/3 has M
        row1 = annotated.iloc[0]
        assert row1["CONSERVATION_SCORE"] is not None

    def test_cli_script_exists(self):
        """Verify CLI script exists."""
        cli_script = Path(__file__).parent.parent / "vicast-analyze" / "add_conservation_scores.py"
        assert cli_script.exists()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
