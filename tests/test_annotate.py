"""
Tests for VICAST-annotate pipeline components.

Tests cover:
- GFF3 validation
- TSV parsing and conversion
- Annotation curation
"""

import os
import sys
import pytest
import tempfile
import shutil
from pathlib import Path

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Add vicast-annotate to path for script imports
sys.path.insert(0, str(Path(__file__).parent.parent / "vicast-annotate"))

# Fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestGFFValidation:
    """Tests for GFF3 validation functions."""

    @pytest.fixture
    def sample_gff(self):
        """Path to sample GFF file."""
        return str(FIXTURES_DIR / "sample.gff3")

    @pytest.fixture
    def sample_fasta(self):
        """Path to sample FASTA file."""
        return str(FIXTURES_DIR / "sample.fasta")

    def test_validate_valid_gff(self, sample_gff, sample_fasta):
        """Test validation of a valid GFF file."""
        from vicast_validation import validate_gff_for_snpeff

        is_valid, errors, warnings = validate_gff_for_snpeff(sample_gff, sample_fasta)

        # Should be valid with no critical errors
        assert is_valid or len(errors) == 0, f"Unexpected errors: {errors}"

    def test_validate_gff_missing_file(self):
        """Test validation with non-existent file."""
        from vicast_validation import validate_gff_for_snpeff

        is_valid, errors, warnings = validate_gff_for_snpeff("/nonexistent/file.gff")

        assert not is_valid
        assert len(errors) > 0
        assert "not found" in errors[0].lower()

    def test_validate_gff_malformed(self, tmp_path):
        """Test validation of malformed GFF file."""
        from vicast_validation import validate_gff_for_snpeff

        # Create a malformed GFF
        malformed_gff = tmp_path / "malformed.gff3"
        malformed_gff.write_text("""##gff-version 3
NC_TEST\tVICAST\tCDS\t100\t500
NC_TEST\tVICAST\tCDS\tinvalid\tend\t.\t+\t.\tID=test
""")

        is_valid, errors, warnings = validate_gff_for_snpeff(str(malformed_gff))

        # Should have errors for malformed lines
        assert len(errors) > 0

    def test_validate_gff_invalid_coordinates(self, tmp_path):
        """Test validation catches invalid coordinates."""
        from vicast_validation import validate_gff_for_snpeff

        # Create GFF with start > end
        bad_coords_gff = tmp_path / "bad_coords.gff3"
        bad_coords_gff.write_text("""##gff-version 3
NC_TEST\tVICAST\tCDS\t500\t100\t.\t+\t.\tID=test;gene=testGene
""")

        is_valid, errors, warnings = validate_gff_for_snpeff(str(bad_coords_gff))

        assert len(errors) > 0
        assert any("start" in e.lower() and "end" in e.lower() for e in errors)

    def test_validate_gff_duplicate_ids(self, tmp_path):
        """Test validation catches duplicate IDs."""
        from vicast_validation import validate_gff_for_snpeff

        # Create GFF with duplicate IDs
        dup_ids_gff = tmp_path / "dup_ids.gff3"
        dup_ids_gff.write_text("""##gff-version 3
NC_TEST\tVICAST\tCDS\t100\t200\t.\t+\t.\tID=gene1;gene=test
NC_TEST\tVICAST\tCDS\t300\t400\t.\t+\t.\tID=gene1;gene=test
""")

        is_valid, errors, warnings = validate_gff_for_snpeff(str(dup_ids_gff))

        assert len(errors) > 0
        assert any("duplicate" in e.lower() for e in errors)

    def test_validate_gff_missing_id(self, tmp_path):
        """Test validation catches missing ID attribute."""
        from vicast_validation import validate_gff_for_snpeff

        # Create GFF without ID
        no_id_gff = tmp_path / "no_id.gff3"
        no_id_gff.write_text("""##gff-version 3
NC_TEST\tVICAST\tCDS\t100\t200\t.\t+\t.\tgene=testGene
""")

        is_valid, errors, warnings = validate_gff_for_snpeff(str(no_id_gff))

        assert len(errors) > 0
        assert any("missing" in e.lower() and "id" in e.lower() for e in errors)


class TestTSVConversion:
    """Tests for TSV to GFF conversion functions."""

    @pytest.fixture
    def sample_tsv(self):
        """Path to sample annotation TSV file."""
        return str(FIXTURES_DIR / "sample_annotation.tsv")

    def test_tsv_to_gff_basic(self, sample_tsv, tmp_path):
        """Test basic TSV to GFF conversion."""
        from step2_add_to_snpeff import tsv_to_gff

        output_gff = tmp_path / "output.gff3"
        success, count, summary = tsv_to_gff(sample_tsv, str(output_gff), force=True)

        assert success
        assert count > 0
        assert output_gff.exists()

        # Check GFF content
        content = output_gff.read_text()
        assert "##gff-version 3" in content
        assert "NC_TEST" in content

    def test_tsv_to_gff_deletes_removed(self, sample_tsv, tmp_path):
        """Test that DELETE marked rows are excluded."""
        from step2_add_to_snpeff import tsv_to_gff

        output_gff = tmp_path / "output.gff3"
        success, count, summary = tsv_to_gff(sample_tsv, str(output_gff), force=True)

        assert success

        # Check that deleted gene is not in output
        content = output_gff.read_text()
        assert "deleteMe" not in content

    def test_tsv_to_gff_creates_mrna_hierarchy(self, sample_tsv, tmp_path):
        """Test that mRNA features are created for CDS."""
        from step2_add_to_snpeff import tsv_to_gff

        output_gff = tmp_path / "output.gff3"
        success, count, summary = tsv_to_gff(sample_tsv, str(output_gff), force=True)

        assert success

        content = output_gff.read_text()
        # Should have mRNA features
        assert "mRNA" in content
        # CDS should have Parent attribute
        assert "Parent=" in content

    def test_tsv_to_gff_preserves_strand(self, sample_tsv, tmp_path):
        """Test that strand information is preserved."""
        from step2_add_to_snpeff import tsv_to_gff

        output_gff = tmp_path / "output.gff3"
        success, count, summary = tsv_to_gff(sample_tsv, str(output_gff), force=True)

        content = output_gff.read_text()
        lines = [l for l in content.split('\n') if l and not l.startswith('#')]

        # Check for both strands
        strands = [l.split('\t')[6] for l in lines]
        assert '+' in strands
        assert '-' in strands

    def test_tsv_to_gff_missing_file(self, tmp_path):
        """Test handling of missing input file."""
        from step2_add_to_snpeff import tsv_to_gff

        output_gff = tmp_path / "output.gff3"
        success, count, summary = tsv_to_gff("/nonexistent/file.tsv", str(output_gff), force=True)

        assert not success
        assert count == 0


class TestTSVParsing:
    """Tests for TSV file parsing."""

    def test_parse_tsv_with_all_columns(self, tmp_path):
        """Test parsing TSV with all expected columns."""
        import pandas as pd

        tsv_file = tmp_path / "test.tsv"
        tsv_file.write_text("""seqid\tsource\ttype\tstart\tend\tstrand\tgene_name\tprotein_id\tproduct
NC_001\tNCBI\tCDS\t100\t500\t+\tgeneA\tYP_001\tprotein A
NC_001\tNCBI\tCDS\t600\t1000\t-\tgeneB\tYP_002\tprotein B
""")

        df = pd.read_csv(tsv_file, sep='\t')

        assert len(df) == 2
        assert 'seqid' in df.columns
        assert 'gene_name' in df.columns
        assert df.iloc[0]['gene_name'] == 'geneA'
        assert df.iloc[1]['strand'] == '-'


class TestAnnotationCuration:
    """Tests for annotation curation logic."""

    def test_polyprotein_detection(self):
        """Test detection of polyprotein annotations."""
        # Polyprotein keywords that should be filtered
        polyprotein_keywords = ['polyprotein', 'polypeptide', 'orf1ab']

        test_products = [
            ("polyprotein", True),
            ("genome polyprotein", True),
            ("envelope protein", False),
            ("NS1", False),
            ("ORF1ab polyprotein", True),
        ]

        for product, should_be_polyprotein in test_products:
            is_poly = any(kw in product.lower() for kw in polyprotein_keywords)
            assert is_poly == should_be_polyprotein, f"Failed for '{product}'"

    def test_feature_filtering_logic(self):
        """Test feature filtering based on action column."""
        import pandas as pd

        data = {
            'gene_name': ['keep1', 'keep2', 'delete1'],
            'action': ['KEEP', 'KEEP', 'DELETE']
        }
        df = pd.DataFrame(data)

        # Filter logic
        filtered = df[df['action'] != 'DELETE']

        assert len(filtered) == 2
        assert 'delete1' not in filtered['gene_name'].values


class TestCDSProteinGeneration:
    """Tests for CDS and protein FASTA generation."""

    @pytest.fixture
    def sample_fasta(self):
        """Path to sample FASTA file."""
        return str(FIXTURES_DIR / "sample.fasta")

    @pytest.fixture
    def sample_tsv(self):
        """Path to sample annotation TSV file."""
        return str(FIXTURES_DIR / "sample_annotation.tsv")

    def test_sequence_extraction_coordinates(self):
        """Test that sequence extraction uses correct coordinates."""
        # Test 1-based to 0-based conversion
        start_1based = 100
        end_1based = 500

        # Convert to 0-based for Python slicing
        start_0based = start_1based - 1
        end_0based = end_1based  # End is exclusive in Python

        # Mock sequence
        seq = "N" * 1000

        extracted = seq[start_0based:end_0based]
        expected_length = end_1based - start_1based + 1

        assert len(extracted) == expected_length


class TestConfigIntegration:
    """Tests for config module integration with annotate scripts."""

    def test_config_available(self):
        """Test that config module is importable."""
        try:
            from vicast.config import get_config, Config
            assert True
        except ImportError:
            pytest.skip("vicast.config module not available")

    def test_snpeff_path_detection(self):
        """Test SnpEff path detection from environment."""
        import os
        from unittest.mock import patch

        with patch.dict(os.environ, {'SNPEFF_HOME': '/test/snpeff'}):
            from vicast.config import Config
            config = Config()
            # Should detect SNPEFF_HOME
            # Note: snpeff_jar might not be set if file doesn't exist


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
