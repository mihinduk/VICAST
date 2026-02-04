"""Tests for VICAST configuration module."""

import os
import pytest
from pathlib import Path
from unittest.mock import patch

# Add src to path for imports
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from vicast.config import Config, get_config, reset_config


class TestConfig:
    """Test suite for Config class."""

    def setup_method(self):
        """Reset config before each test."""
        reset_config()

    def teardown_method(self):
        """Clean up after each test."""
        reset_config()

    def test_config_creation(self):
        """Test basic config creation."""
        config = Config()
        assert config is not None
        assert config.threads == 4  # default value

    def test_config_singleton(self):
        """Test that get_config returns singleton."""
        config1 = get_config()
        config2 = get_config()
        assert config1 is config2

    def test_reset_config(self):
        """Test config reset."""
        config1 = get_config()
        reset_config()
        config2 = get_config()
        assert config1 is not config2

    @patch.dict(os.environ, {"SNPEFF_JAR": "/path/to/snpEff.jar"})
    def test_env_snpeff_jar(self):
        """Test loading SNPEFF_JAR from environment."""
        config = Config()
        assert config.snpeff_jar == Path("/path/to/snpEff.jar")

    @patch.dict(os.environ, {"SNPEFF_DATA": "/path/to/data"})
    def test_env_snpeff_data(self):
        """Test loading SNPEFF_DATA from environment."""
        config = Config()
        assert config.snpeff_data == Path("/path/to/data")

    @patch.dict(os.environ, {"VICAST_THREADS": "8"})
    def test_env_threads(self):
        """Test loading VICAST_THREADS from environment."""
        config = Config()
        assert config.threads == 8

    @patch.dict(os.environ, {"VICAST_THREADS": "invalid"})
    def test_env_threads_invalid(self):
        """Test invalid VICAST_THREADS value."""
        config = Config()
        assert config.threads == 4  # should fall back to default

    @patch.dict(os.environ, {"NCBI_EMAIL": "test@example.com"})
    def test_env_ncbi_email(self):
        """Test loading NCBI_EMAIL from environment."""
        config = Config()
        assert config.ncbi_email == "test@example.com"

    def test_to_dict(self):
        """Test exporting config as dictionary."""
        config = Config()
        d = config.to_dict()
        assert isinstance(d, dict)
        assert "threads" in d
        assert "ncbi_email" in d

    def test_to_shell_exports(self):
        """Test generating shell export statements."""
        config = Config()
        config.threads = 8
        exports = config.to_shell_exports()
        assert "export VICAST_THREADS=" in exports

    def test_validate_missing_snpeff(self):
        """Test validation fails without SnpEff."""
        config = Config()
        config.snpeff_jar = None
        config.snpeff_data = None
        is_valid, errors = config.validate(require_snpeff=True)
        assert not is_valid
        assert len(errors) > 0

    def test_validate_skip_snpeff(self):
        """Test validation can skip SnpEff requirement."""
        config = Config()
        config.snpeff_jar = None
        config.snpeff_data = None
        # Mock java being available
        config.java_path = Path("/usr/bin/java")
        with patch.object(Path, "exists", return_value=True):
            is_valid, errors = config.validate(require_snpeff=False)
            # Should only fail if java doesn't exist
            assert any("Java" in e for e in errors) or is_valid

    def test_get_snpeff_genome_dir(self):
        """Test getting genome-specific SnpEff directory."""
        config = Config()
        config.snpeff_data = Path("/path/to/data")
        genome_dir = config.get_snpeff_genome_dir("NC_001477")
        assert genome_dir == Path("/path/to/data/NC_001477")

    def test_get_snpeff_genome_dir_not_configured(self):
        """Test error when snpeff_data not configured."""
        config = Config()
        config.snpeff_data = None
        with pytest.raises(ValueError):
            config.get_snpeff_genome_dir("NC_001477")


class TestConfigAutoDetection:
    """Test suite for auto-detection features."""

    def setup_method(self):
        reset_config()

    def teardown_method(self):
        reset_config()

    @patch.dict(os.environ, {"SNPEFF_HOME": "/opt/snpEff"}, clear=False)
    def test_derive_jar_from_home(self):
        """Test deriving snpeff_jar from SNPEFF_HOME."""
        with patch.object(Path, "exists", return_value=True):
            config = Config()
            assert config.snpeff_home == Path("/opt/snpEff")
            # JAR should be derived
            assert config.snpeff_jar == Path("/opt/snpEff/snpEff.jar")

    @patch.dict(os.environ, {}, clear=False)
    def test_scratch_defaults_to_tmpdir(self):
        """Test scratch directory defaults to temp."""
        config = Config()
        assert config.scratch_dir is not None
        # Should be a Path object pointing to system temp
        assert isinstance(config.scratch_dir, Path)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
