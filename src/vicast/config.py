"""
VICAST Configuration Module

Centralized configuration management for VICAST pipelines.
Supports environment variables, config files, and sensible defaults.

Configuration Priority (highest to lowest):
1. Explicit function arguments
2. Environment variables
3. Config file (~/.vicast/config.yml or VICAST_CONFIG)
4. Auto-detected defaults

Environment Variables:
    VICAST_CONFIG       - Path to config file (default: ~/.vicast/config.yml)
    SNPEFF_JAR          - Path to snpEff.jar
    SNPEFF_DATA         - Path to snpEff data directory
    SNPEFF_HOME         - Path to snpEff installation directory
    JAVA_HOME           - Path to Java installation
    VICAST_SCRATCH      - Scratch/temp directory for large files
    VICAST_THREADS      - Default number of threads
    NCBI_EMAIL          - Email for NCBI Entrez queries (required by NCBI)
"""

import os
import shutil
import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

# Singleton config instance
_config_instance: Optional["Config"] = None


@dataclass
class Config:
    """
    VICAST configuration container.

    Attributes:
        snpeff_jar: Path to snpEff.jar executable
        snpeff_data: Path to snpEff data directory
        snpeff_home: Path to snpEff installation root
        java_home: Path to Java installation
        java_path: Path to java executable
        scratch_dir: Scratch directory for temporary/large files
        threads: Default number of threads for parallel operations
        ncbi_email: Email for NCBI Entrez queries
        blast_db: Path to BLAST database directory
        conda_env: Name of conda environment to use
    """

    # SnpEff paths
    snpeff_jar: Optional[Path] = None
    snpeff_data: Optional[Path] = None
    snpeff_home: Optional[Path] = None

    # Java paths
    java_home: Optional[Path] = None
    java_path: Optional[Path] = None

    # Working directories
    scratch_dir: Optional[Path] = None
    work_dir: Optional[Path] = None

    # Resources
    threads: int = 4
    memory_gb: int = 8

    # External services
    ncbi_email: str = "vicast_user@example.com"

    # Optional tool paths
    blast_db: Optional[Path] = None
    conda_env: Optional[str] = None

    # Internal state
    _config_file: Optional[Path] = field(default=None, repr=False)
    _initialized: bool = field(default=False, repr=False)

    def __post_init__(self):
        """Initialize configuration from environment and auto-detection."""
        if not self._initialized:
            self._load_from_environment()
            self._auto_detect_paths()
            self._initialized = True

    def _load_from_environment(self) -> None:
        """Load configuration from environment variables."""

        # SnpEff configuration
        if os.environ.get("SNPEFF_JAR"):
            self.snpeff_jar = Path(os.environ["SNPEFF_JAR"])
        if os.environ.get("SNPEFF_DATA"):
            self.snpeff_data = Path(os.environ["SNPEFF_DATA"])
        if os.environ.get("SNPEFF_HOME"):
            self.snpeff_home = Path(os.environ["SNPEFF_HOME"])

        # Java configuration
        if os.environ.get("JAVA_HOME"):
            self.java_home = Path(os.environ["JAVA_HOME"])
            java_bin = self.java_home / "bin" / "java"
            if java_bin.exists():
                self.java_path = java_bin

        # Working directories
        if os.environ.get("VICAST_SCRATCH"):
            self.scratch_dir = Path(os.environ["VICAST_SCRATCH"])
        elif os.environ.get("SCRATCH_DIR"):
            self.scratch_dir = Path(os.environ["SCRATCH_DIR"])
        elif os.environ.get("TMPDIR"):
            self.scratch_dir = Path(os.environ["TMPDIR"])

        # Resources
        if os.environ.get("VICAST_THREADS"):
            try:
                self.threads = int(os.environ["VICAST_THREADS"])
            except ValueError:
                pass

        # NCBI
        if os.environ.get("NCBI_EMAIL"):
            self.ncbi_email = os.environ["NCBI_EMAIL"]
        elif os.environ.get("ENTREZ_EMAIL"):
            self.ncbi_email = os.environ["ENTREZ_EMAIL"]

        # BLAST database
        if os.environ.get("BLASTDB"):
            self.blast_db = Path(os.environ["BLASTDB"])

        # Conda environment
        if os.environ.get("VICAST_CONDA_ENV"):
            self.conda_env = os.environ["VICAST_CONDA_ENV"]
        elif os.environ.get("CONDA_DEFAULT_ENV"):
            self.conda_env = os.environ["CONDA_DEFAULT_ENV"]

    def _auto_detect_paths(self) -> None:
        """Auto-detect paths for tools if not explicitly configured."""

        # Auto-detect Java
        if not self.java_path:
            java_cmd = shutil.which("java")
            if java_cmd:
                self.java_path = Path(java_cmd)

        # Auto-detect SnpEff from snpeff_home
        if self.snpeff_home and not self.snpeff_jar:
            candidate = self.snpeff_home / "snpEff.jar"
            if candidate.exists():
                self.snpeff_jar = candidate

        if self.snpeff_home and not self.snpeff_data:
            candidate = self.snpeff_home / "data"
            if candidate.exists():
                self.snpeff_data = candidate

        # Derive snpeff_home from snpeff_jar if needed
        if self.snpeff_jar and not self.snpeff_home:
            self.snpeff_home = self.snpeff_jar.parent
            if not self.snpeff_data:
                candidate = self.snpeff_home / "data"
                if candidate.exists():
                    self.snpeff_data = candidate

        # Default scratch to temp directory
        if not self.scratch_dir:
            import tempfile
            self.scratch_dir = Path(tempfile.gettempdir())

        # Default work_dir to current directory
        if not self.work_dir:
            self.work_dir = Path.cwd()

    def validate(self, require_snpeff: bool = True) -> tuple[bool, list[str]]:
        """
        Validate configuration.

        Args:
            require_snpeff: Whether SnpEff paths are required

        Returns:
            Tuple of (is_valid, list of error messages)
        """
        errors = []

        if require_snpeff:
            if not self.snpeff_jar:
                errors.append("SNPEFF_JAR not configured. Set SNPEFF_JAR environment variable.")
            elif not self.snpeff_jar.exists():
                errors.append(f"SNPEFF_JAR not found: {self.snpeff_jar}")

            if not self.snpeff_data:
                errors.append("SNPEFF_DATA not configured. Set SNPEFF_DATA environment variable.")
            elif not self.snpeff_data.exists():
                errors.append(f"SNPEFF_DATA directory not found: {self.snpeff_data}")

        if not self.java_path:
            errors.append("Java not found. Set JAVA_HOME or ensure java is in PATH.")
        elif not self.java_path.exists():
            errors.append(f"Java executable not found: {self.java_path}")

        return len(errors) == 0, errors

    def get_snpeff_genome_dir(self, genome_id: str) -> Path:
        """Get the SnpEff data directory for a specific genome."""
        if not self.snpeff_data:
            raise ValueError("SNPEFF_DATA not configured")
        return self.snpeff_data / genome_id

    def genome_exists_in_snpeff(self, genome_id: str) -> bool:
        """Check if a genome exists in the SnpEff database."""
        genome_dir = self.get_snpeff_genome_dir(genome_id)
        predictor = genome_dir / "snpEffectPredictor.bin"
        return predictor.exists()

    def to_dict(self) -> Dict[str, Any]:
        """Export configuration as dictionary."""
        return {
            "snpeff_jar": str(self.snpeff_jar) if self.snpeff_jar else None,
            "snpeff_data": str(self.snpeff_data) if self.snpeff_data else None,
            "snpeff_home": str(self.snpeff_home) if self.snpeff_home else None,
            "java_home": str(self.java_home) if self.java_home else None,
            "java_path": str(self.java_path) if self.java_path else None,
            "scratch_dir": str(self.scratch_dir) if self.scratch_dir else None,
            "work_dir": str(self.work_dir) if self.work_dir else None,
            "threads": self.threads,
            "memory_gb": self.memory_gb,
            "ncbi_email": self.ncbi_email,
            "blast_db": str(self.blast_db) if self.blast_db else None,
            "conda_env": self.conda_env,
        }

    def to_shell_exports(self) -> str:
        """Generate shell export statements for this configuration."""
        lines = ["# VICAST Configuration Exports"]

        if self.snpeff_jar:
            lines.append(f'export SNPEFF_JAR="{self.snpeff_jar}"')
        if self.snpeff_data:
            lines.append(f'export SNPEFF_DATA="{self.snpeff_data}"')
        if self.snpeff_home:
            lines.append(f'export SNPEFF_HOME="{self.snpeff_home}"')
        if self.java_home:
            lines.append(f'export JAVA_HOME="{self.java_home}"')
        if self.scratch_dir:
            lines.append(f'export VICAST_SCRATCH="{self.scratch_dir}"')
        lines.append(f'export VICAST_THREADS="{self.threads}"')
        if self.ncbi_email != "vicast_user@example.com":
            lines.append(f'export NCBI_EMAIL="{self.ncbi_email}"')

        return "\n".join(lines)

    def print_status(self) -> None:
        """Print configuration status to stdout."""
        print("VICAST Configuration Status")
        print("=" * 50)

        def status_icon(path: Optional[Path]) -> str:
            if path is None:
                return "[ ] Not configured"
            elif path.exists():
                return f"[✓] {path}"
            else:
                return f"[✗] {path} (NOT FOUND)"

        print(f"SnpEff JAR:  {status_icon(self.snpeff_jar)}")
        print(f"SnpEff Data: {status_icon(self.snpeff_data)}")
        print(f"Java:        {status_icon(self.java_path)}")
        print(f"Scratch:     {status_icon(self.scratch_dir)}")
        print(f"Threads:     {self.threads}")
        print(f"NCBI Email:  {self.ncbi_email}")
        print("=" * 50)


def get_config() -> Config:
    """
    Get the global configuration instance.

    Returns:
        Config: The singleton configuration instance
    """
    global _config_instance
    if _config_instance is None:
        _config_instance = Config()
    return _config_instance


def reset_config() -> None:
    """Reset the global configuration instance (useful for testing)."""
    global _config_instance
    _config_instance = None


def print_setup_instructions() -> None:
    """Print setup instructions for users."""
    print("""
VICAST Configuration Setup
===========================

VICAST requires the following to be configured:

1. SNPEFF INSTALLATION (Required for annotation)

   Option A - Environment Variables (Recommended):

     export SNPEFF_HOME="/path/to/snpEff"
     export SNPEFF_JAR="${SNPEFF_HOME}/snpEff.jar"
     export SNPEFF_DATA="${SNPEFF_HOME}/data"

   Option B - If SnpEff not installed:

     # Download SnpEff
     wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
     unzip snpEff_latest_core.zip
     export SNPEFF_HOME="$(pwd)/snpEff"

2. JAVA 21+ (Required for SnpEff 5.2+)

   Most systems:
     export JAVA_HOME="/path/to/jdk-21"

   Or install via conda:
     conda install -c conda-forge openjdk=21

3. NCBI EMAIL (Required for genome downloads)

   export NCBI_EMAIL="your.email@institution.edu"

4. SCRATCH DIRECTORY (Recommended for HPC)

   export VICAST_SCRATCH="/scratch/$USER/vicast"

After setting environment variables, verify with:

   python -c "from vicast.config import get_config; get_config().print_status()"

For persistent configuration, add exports to ~/.bashrc or ~/.zshrc
""")


if __name__ == "__main__":
    # When run directly, print status
    config = get_config()
    config.print_status()

    is_valid, errors = config.validate()
    if not is_valid:
        print("\nConfiguration Errors:")
        for error in errors:
            print(f"  - {error}")
        print()
        print_setup_instructions()
