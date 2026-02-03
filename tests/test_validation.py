"""Tests for the validation module."""

import os
import tempfile
import pytest
from vicast.validation import (
    validate_vcf_ref_bases,
    validate_vcf,
    validate_gff_for_snpeff,
    _parse_fasta_sequences,
)


class TestParseFastaSequences:
    """Tests for _parse_fasta_sequences helper."""

    def test_parse_single_sequence(self, tmp_path):
        """Test parsing a single sequence."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACGT\nTGCA\n")

        result = _parse_fasta_sequences(str(fasta_file))

        assert result == {"seq1": "ACGTTGCA"}

    def test_parse_multiple_sequences(self, tmp_path):
        """Test parsing multiple sequences."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nACGT\n>seq2\nTTTT\nAAAA\n")

        result = _parse_fasta_sequences(str(fasta_file))

        assert result == {"seq1": "ACGT", "seq2": "TTTTAAAA"}

    def test_parse_sequence_with_description(self, tmp_path):
        """Test that only the ID (first word) is used."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1 description here\nACGT\n")

        result = _parse_fasta_sequences(str(fasta_file))

        assert "seq1" in result
        assert "seq1 description here" not in result

    def test_uppercase_conversion(self, tmp_path):
        """Test that sequences are converted to uppercase."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nacgt\n")

        result = _parse_fasta_sequences(str(fasta_file))

        assert result["seq1"] == "ACGT"


class TestValidateVcfRefBases:
    """Tests for validate_vcf_ref_bases function."""

    def test_matching_ref_bases(self, tmp_path):
        """Test validation passes when REF bases match."""
        # Create reference FASTA
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGTACGTACGT\n")

        # Create VCF with matching REF
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tA\tG\t100\tPASS\t.\n"
            "chr1\t5\t.\tA\tT\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        assert is_valid
        assert len(errors) == 0

    def test_mismatching_ref_bases(self, tmp_path):
        """Test validation fails when REF bases don't match."""
        # Create reference FASTA
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGTACGTACGT\n")

        # Create VCF with wrong REF (position 1 is A, not T)
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tT\tG\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        assert not is_valid
        assert len(errors) == 1
        assert "REF mismatch" in errors[0]
        assert "VCF has 'T'" in errors[0]
        assert "reference has 'A'" in errors[0]

    def test_multi_base_ref(self, tmp_path):
        """Test validation of multi-base REF alleles."""
        # Create reference FASTA
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGTACGTACGT\n")

        # Create VCF with multi-base REF (positions 2-4 should be CGT)
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t2\t.\tCGT\tC\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        assert is_valid
        assert len(errors) == 0

    def test_unknown_chromosome_warning(self, tmp_path):
        """Test warning when chromosome not in reference."""
        # Create reference FASTA with only chr1
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGTACGTACGT\n")

        # Create VCF with chr2
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr2\t1\t.\tA\tG\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        # Should be valid (no errors) but with warning
        assert is_valid
        assert len(errors) == 0
        assert any("chr2" in w and "not found" in w for w in warnings)

    def test_position_beyond_sequence(self, tmp_path):
        """Test error when variant position exceeds sequence length."""
        # Create short reference
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGT\n")  # Only 4 bases

        # Create VCF with position beyond end
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t10\t.\tA\tG\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        assert not is_valid
        assert len(errors) == 1
        assert "extends beyond" in errors[0]

    def test_missing_vcf_file(self, tmp_path):
        """Test error when VCF file doesn't exist."""
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGT\n")

        is_valid, errors, warnings = validate_vcf_ref_bases(
            "/nonexistent/file.vcf", str(fasta_file)
        )

        assert not is_valid
        assert any("not found" in e for e in errors)

    def test_missing_fasta_file(self, tmp_path):
        """Test error when FASTA file doesn't exist."""
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text("##fileformat=VCFv4.2\n")

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), "/nonexistent/ref.fasta"
        )

        assert not is_valid
        assert any("not found" in e for e in errors)

    def test_case_insensitive_comparison(self, tmp_path):
        """Test that comparison is case-insensitive."""
        # Create reference with lowercase
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nacgt\n")

        # Create VCF with uppercase REF
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tA\tG\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        assert is_valid

    def test_multiple_chromosomes(self, tmp_path):
        """Test validation across multiple chromosomes."""
        # Create multi-sequence reference
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nACGT\n>chr2\nTTTT\n>chr3\nGGGG\n")

        # Create VCF with variants on multiple chroms
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1\t.\tA\tG\t100\tPASS\t.\n"
            "chr2\t2\t.\tT\tA\t100\tPASS\t.\n"
            "chr3\t3\t.\tG\tC\t100\tPASS\t.\n"
        )

        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file)
        )

        assert is_valid
        assert len(errors) == 0

    def test_max_errors_limit(self, tmp_path):
        """Test that validation stops after max_errors."""
        # Create reference
        fasta_file = tmp_path / "ref.fasta"
        fasta_file.write_text(">chr1\nAAAAAAAAAAAA\n")

        # Create VCF with many wrong REFs
        vcf_lines = ["##fileformat=VCFv4.2\n",
                     "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"]
        for i in range(1, 11):
            vcf_lines.append(f"chr1\t{i}\t.\tT\tG\t100\tPASS\t.\n")  # All wrong

        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text("".join(vcf_lines))

        # Set max_errors to 3
        is_valid, errors, warnings = validate_vcf_ref_bases(
            str(vcf_file), str(fasta_file), max_errors=3
        )

        assert not is_valid
        assert len(errors) == 3
        assert any("Stopped after" in w for w in warnings)


class TestValidateVcf:
    """Tests for basic validate_vcf function."""

    def test_valid_vcf(self, tmp_path):
        """Test validation of a valid VCF."""
        vcf_file = tmp_path / "test.vcf"
        vcf_file.write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t100\tPASS\tDP=50\n"
        )

        is_valid, errors, warnings = validate_vcf(str(vcf_file))

        assert is_valid
        assert len(errors) == 0

    def test_missing_file(self):
        """Test validation of nonexistent file."""
        is_valid, errors, warnings = validate_vcf("/nonexistent/file.vcf")

        assert not is_valid
        assert any("not found" in e for e in errors)
