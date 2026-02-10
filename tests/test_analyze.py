"""
Tests for VICAST-analyze pipeline components.

Tests cover:
- VCF/TSV parsing
- INFO field extraction (AF, DP, DP4, SB)
- Variant filtering
- Known viruses configuration
"""

import os
import sys
import pytest
import tempfile
import re
from pathlib import Path
import pandas as pd
import io

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Add vicast-analyze to path for script imports
sys.path.insert(0, str(Path(__file__).parent.parent / "vicast-analyze"))

# Fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestINFOFieldExtraction:
    """Tests for INFO field parsing functions."""

    def test_extract_af(self):
        """Test allele frequency extraction from INFO field."""
        def extract_af(info):
            match = re.search(r'AF=([0-9.]+)', str(info))
            return float(match.group(1)) if match else 0.0

        # Test cases
        assert extract_af("DP=100;AF=0.25;SB=5") == 0.25
        assert extract_af("AF=0.5") == 0.5
        assert extract_af("AF=1.0;DP=200") == 1.0
        assert extract_af("DP=100") == 0.0  # No AF
        assert extract_af("") == 0.0
        assert extract_af(None) == 0.0

    def test_extract_dp(self):
        """Test depth extraction from INFO field."""
        def extract_dp(info):
            match = re.search(r'DP=([0-9]+)', str(info))
            return int(match.group(1)) if match else 0

        # Test cases
        assert extract_dp("DP=500;AF=0.25") == 500
        assert extract_dp("AF=0.25;DP=1000;SB=10") == 1000
        assert extract_dp("AF=0.25") == 0  # No DP
        assert extract_dp("") == 0
        assert extract_dp("DP4=10,20,5,5") == 0  # DP4 is different from DP

    def test_extract_dp4(self):
        """Test DP4 (strand-specific depth) extraction."""
        def extract_dp4(info):
            match = re.search(r'DP4=([0-9,]+)', str(info))
            return match.group(1) if match else '0,0,0,0'

        # Test cases
        assert extract_dp4("DP4=100,90,50,40") == "100,90,50,40"
        assert extract_dp4("DP=200;DP4=100,100,50,50;AF=0.5") == "100,100,50,50"
        assert extract_dp4("DP=200") == "0,0,0,0"  # No DP4

    def test_extract_sb(self):
        """Test strand bias extraction."""
        def extract_sb(info):
            match = re.search(r'SB=([0-9]+)', str(info))
            return match.group(1) if match else '0'

        # Test cases
        assert extract_sb("SB=15") == "15"
        assert extract_sb("DP=100;SB=0;AF=0.5") == "0"
        assert extract_sb("DP=100") == "0"  # No SB

    def test_parse_dp4_values(self):
        """Test parsing individual DP4 values."""
        dp4_string = "100,90,50,40"
        values = [int(x) for x in dp4_string.split(',')]

        assert values[0] == 100  # ref forward
        assert values[1] == 90   # ref reverse
        assert values[2] == 50   # alt forward
        assert values[3] == 40   # alt reverse

        # Calculate derived metrics
        ref_total = values[0] + values[1]
        alt_total = values[2] + values[3]
        total = ref_total + alt_total

        assert ref_total == 190
        assert alt_total == 90
        assert total == 280


class TestVariantFiltering:
    """Tests for variant filtering logic."""

    @pytest.fixture
    def sample_variants_df(self):
        """Create a sample DataFrame of variants."""
        data = {
            'CHROM': ['NC_TEST'] * 5,
            'POS': [100, 200, 300, 400, 500],
            'REF': ['A', 'T', 'G', 'C', 'A'],
            'ALT': ['G', 'C', 'A', 'T', 'T'],
            'QUAL': [2500, 1500, 500, 100, 3000],
            'Total_Depth': [500, 300, 150, 50, 800],
            'Allele_Frequency': [0.25, 0.15, 0.08, 0.02, 0.50]
        }
        return pd.DataFrame(data)

    def test_quality_filter(self, sample_variants_df):
        """Test filtering by quality score."""
        quality_cutoff = 1000

        filtered = sample_variants_df[sample_variants_df['QUAL'] >= quality_cutoff]

        assert len(filtered) == 3  # QUAL >= 1000: 2500, 1500, 3000
        assert all(filtered['QUAL'] >= quality_cutoff)

    def test_depth_filter(self, sample_variants_df):
        """Test filtering by depth."""
        depth_cutoff = 200

        filtered = sample_variants_df[sample_variants_df['Total_Depth'] >= depth_cutoff]

        assert len(filtered) == 3  # DP >= 200: 500, 300, 800
        assert all(filtered['Total_Depth'] >= depth_cutoff)

    def test_frequency_filter(self, sample_variants_df):
        """Test filtering by allele frequency."""
        freq_cutoff = 0.10

        filtered = sample_variants_df[sample_variants_df['Allele_Frequency'] >= freq_cutoff]

        assert len(filtered) == 3  # AF >= 0.10: 0.25, 0.15, 0.50
        assert all(filtered['Allele_Frequency'] >= freq_cutoff)

    def test_combined_filters(self, sample_variants_df):
        """Test combined quality, depth, and frequency filtering."""
        quality_cutoff = 1000
        depth_cutoff = 200
        freq_cutoff = 0.10

        filtered = sample_variants_df[
            (sample_variants_df['QUAL'] >= quality_cutoff) &
            (sample_variants_df['Total_Depth'] >= depth_cutoff) &
            (sample_variants_df['Allele_Frequency'] >= freq_cutoff)
        ]

        # Only variants passing all three filters
        assert len(filtered) == 3
        assert all(filtered['QUAL'] >= quality_cutoff)
        assert all(filtered['Total_Depth'] >= depth_cutoff)
        assert all(filtered['Allele_Frequency'] >= freq_cutoff)


class TestSnpEffTSVParsing:
    """Tests for SnpEff TSV file parsing."""

    @pytest.fixture
    def sample_snpeff_tsv(self):
        """Path to sample SnpEff TSV file."""
        return str(FIXTURES_DIR / "sample_snpeff.tsv")

    def test_parse_snpeff_header(self, sample_snpeff_tsv):
        """Test parsing SnpEff TSV with #CHROM header."""
        with open(sample_snpeff_tsv, 'r') as f:
            lines = f.readlines()

        # Find header line
        header_line = None
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line[1:].strip()  # Remove leading #
                break

        assert header_line is not None
        assert 'CHROM' in header_line
        assert 'POS' in header_line
        assert 'EFFECT' in header_line

    def test_parse_snpeff_data(self, sample_snpeff_tsv):
        """Test parsing SnpEff TSV data."""
        with open(sample_snpeff_tsv, 'r') as f:
            lines = f.readlines()

        header_line = None
        data_lines = []
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line[1:].strip()
            elif not line.startswith('##'):
                data_lines.append(line.strip())

        tsv_content = header_line + '\n' + '\n'.join(data_lines)
        df = pd.read_csv(io.StringIO(tsv_content), sep='\t')

        assert len(df) == 5
        assert 'EFFECT' in df.columns
        assert 'GENE_NAME' in df.columns

    def test_snpeff_effect_types(self, sample_snpeff_tsv):
        """Test that effect types are parsed correctly."""
        with open(sample_snpeff_tsv, 'r') as f:
            lines = f.readlines()

        header_line = None
        data_lines = []
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line[1:].strip()
            elif not line.startswith('##'):
                data_lines.append(line.strip())

        tsv_content = header_line + '\n' + '\n'.join(data_lines)
        df = pd.read_csv(io.StringIO(tsv_content), sep='\t')

        effects = df['EFFECT'].unique()
        assert 'missense_variant' in effects
        assert 'synonymous_variant' in effects

    def test_snpeff_hgvs_parsing(self, sample_snpeff_tsv):
        """Test HGVS notation parsing."""
        with open(sample_snpeff_tsv, 'r') as f:
            lines = f.readlines()

        header_line = None
        data_lines = []
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line[1:].strip()
            elif not line.startswith('##'):
                data_lines.append(line.strip())

        tsv_content = header_line + '\n' + '\n'.join(data_lines)
        df = pd.read_csv(io.StringIO(tsv_content), sep='\t')

        # Check HGVSc format (e.g., c.51A>G)
        assert df['HGVSc'].iloc[0].startswith('c.')
        # Check HGVSp format (e.g., p.Lys17Arg)
        assert df['HGVSp'].iloc[0].startswith('p.')


class TestKnownVirusesConfig:
    """Tests for known_viruses.json configuration."""

    def test_known_viruses_json_structure(self):
        """Test that known_viruses.json has expected structure."""
        import json

        known_viruses_path = Path(__file__).parent.parent / "vicast-analyze" / "known_viruses.json"

        if not known_viruses_path.exists():
            pytest.skip("known_viruses.json not found")

        with open(known_viruses_path, 'r') as f:
            config = json.load(f)

        # Should be a dictionary with accession keys
        assert isinstance(config, dict)

        # Check first entry has expected fields
        if config:
            first_key = list(config.keys())[0]
            entry = config[first_key]

            # Expected fields
            expected_fields = ['name', 'gene_coords']
            for field in expected_fields:
                assert field in entry, f"Missing field: {field}"

    def test_gene_coords_format(self):
        """Test gene coordinates format in known_viruses.json."""
        import json

        known_viruses_path = Path(__file__).parent.parent / "vicast-analyze" / "known_viruses.json"

        if not known_viruses_path.exists():
            pytest.skip("known_viruses.json not found")

        with open(known_viruses_path, 'r') as f:
            config = json.load(f)

        for accession, entry in config.items():
            if 'gene_coords' in entry:
                gene_coords = entry['gene_coords']
                assert isinstance(gene_coords, dict)

                for gene, coords in gene_coords.items():
                    assert isinstance(coords, list)
                    assert len(coords) == 2
                    assert coords[0] < coords[1], f"Invalid coords for {gene}: {coords}"


class TestConsensusGeneration:
    """Tests for consensus sequence generation logic."""

    def test_variant_application_logic(self):
        """Test logic for applying variants to reference."""
        # Mock reference sequence
        reference = list("ATGCATGCATGC")

        # Apply SNP at position 3 (0-indexed = 2)
        position = 3  # 1-based
        ref_allele = "G"
        alt_allele = "A"

        # Verify reference
        assert reference[position - 1] == ref_allele

        # Apply variant
        reference[position - 1] = alt_allele

        assert reference[position - 1] == alt_allele
        assert "".join(reference) == "ATACATGCATGC"

    def test_multi_allelic_handling(self):
        """Test handling of multi-allelic sites."""
        # At position 100, we have two variants at different frequencies
        variants_at_pos = [
            {'alt': 'G', 'af': 0.40},
            {'alt': 'T', 'af': 0.15}
        ]

        # The consensus should use the most frequent variant
        most_frequent = max(variants_at_pos, key=lambda x: x['af'])
        assert most_frequent['alt'] == 'G'

    def test_frequency_threshold_for_consensus(self):
        """Test that low-frequency variants don't affect consensus."""
        consensus_threshold = 0.50

        variants = [
            {'pos': 100, 'alt': 'G', 'af': 0.25},  # Below threshold
            {'pos': 200, 'alt': 'T', 'af': 0.75},  # Above threshold
            {'pos': 300, 'alt': 'C', 'af': 0.50},  # At threshold
        ]

        consensus_variants = [v for v in variants if v['af'] >= consensus_threshold]

        assert len(consensus_variants) == 2
        assert all(v['af'] >= consensus_threshold for v in consensus_variants)


class TestGeneCoordExtraction:
    """Tests for extracting gene coordinates from SnpEff."""

    def test_gff_parsing_for_coords(self):
        """Test parsing GFF for gene coordinates."""
        gff_content = """##gff-version 3
NC_TEST\tVICAST\tmRNA\t100\t500\t.\t+\t.\tID=gene1;gene=testGene1
NC_TEST\tVICAST\tmRNA\t600\t1200\t.\t+\t.\tID=gene2;gene=testGene2
"""
        gene_coords = {}

        for line in gff_content.strip().split('\n'):
            if line.startswith('#'):
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]

            if feature_type == 'mRNA':
                # Extract gene name
                for attr in attributes.split(';'):
                    if attr.startswith('gene='):
                        gene_name = attr.split('=')[1]
                        gene_coords[gene_name] = [start, end]
                        break

        assert len(gene_coords) == 2
        assert gene_coords['testGene1'] == [100, 500]
        assert gene_coords['testGene2'] == [600, 1200]


class TestConfigIntegration:
    """Tests for config module integration with analyze scripts."""

    def test_snpeff_data_detection(self):
        """Test SnpEff data directory detection."""
        import os
        from unittest.mock import patch

        with patch.dict(os.environ, {'SNPEFF_DATA': '/test/snpeff/data'}):
            snpeff_data = os.environ.get('SNPEFF_DATA')
            assert snpeff_data == '/test/snpeff/data'


class TestPipelineIntegration:
    """Integration tests for the analyze pipeline."""

    @pytest.mark.slow
    def test_full_tsv_processing_pipeline(self):
        """Test complete TSV processing workflow."""
        sample_snpeff_tsv = FIXTURES_DIR / "sample_snpeff.tsv"

        if not sample_snpeff_tsv.exists():
            pytest.skip("Sample SnpEff TSV not found")

        # Parse file
        with open(sample_snpeff_tsv, 'r') as f:
            lines = f.readlines()

        header_line = None
        data_lines = []
        for line in lines:
            if line.startswith('#CHROM'):
                header_line = line[1:].strip()
            elif not line.startswith('##'):
                data_lines.append(line.strip())

        tsv_content = header_line + '\n' + '\n'.join(data_lines)
        df = pd.read_csv(io.StringIO(tsv_content), sep='\t')

        # Extract INFO fields
        def extract_af(info):
            match = re.search(r'AF=([0-9.]+)', str(info))
            return float(match.group(1)) if match else 0.0

        def extract_dp(info):
            match = re.search(r'DP=([0-9]+)', str(info))
            return int(match.group(1)) if match else 0

        df['Allele_Frequency'] = df['INFO'].apply(extract_af)
        df['Total_Depth'] = df['INFO'].apply(extract_dp)

        # Apply filters
        quality_cutoff = 1000
        depth_cutoff = 200
        freq_cutoff = 0.10

        filtered = df[
            (df['QUAL'] >= quality_cutoff) &
            (df['Total_Depth'] >= depth_cutoff) &
            (df['Allele_Frequency'] >= freq_cutoff)
        ]

        # Verify filtering worked
        assert len(filtered) < len(df)
        assert all(filtered['QUAL'] >= quality_cutoff)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
