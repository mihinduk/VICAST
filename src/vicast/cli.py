"""
VICAST Command-Line Interface

Entry points for vicast-annotate, vicast-analyze, and vicast-validate commands.
"""

import sys
from vicast.logging import setup_logging, get_logger


def annotate_main():
    """Entry point for vicast-annotate command."""
    logger = setup_logging("vicast.annotate")

    logger.info("VICAST-Annotate")
    logger.info("===============")
    logger.info("")
    logger.info("For now, please use the Python scripts directly:")
    logger.info("")
    logger.info("  # Check SnpEff database")
    logger.info("  python vicast-annotate/step0_check_snpeff.py <genome_id>")
    logger.info("")
    logger.info("  # Parse genome (Pathway 2)")
    logger.info("  python vicast-annotate/step1_parse_viral_genome.py <genome_id>")
    logger.info("")
    logger.info("  # Add to SnpEff")
    logger.info("  python vicast-annotate/step2_add_to_snpeff.py <genome_id> <tsv_file>")
    logger.info("")
    logger.info("See vicast-annotate/README.md for full documentation.")
    sys.exit(0)


def analyze_main():
    """Entry point for vicast-analyze command."""
    logger = setup_logging("vicast.analyze")

    logger.info("VICAST-Analyze")
    logger.info("==============")
    logger.info("")
    logger.info("For now, please use the shell scripts directly:")
    logger.info("")
    logger.info("  # QC only (Steps 1-6)")
    logger.info("  ./vicast-analyze/run_vicast_analyze_qc_only.sh <R1> <R2> <accession> [threads]")
    logger.info("")
    logger.info("  # Annotation only (Steps 7-9)")
    logger.info("  ./vicast-analyze/run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>")
    logger.info("")
    logger.info("  # Full pipeline")
    logger.info("  ./vicast-analyze/run_vicast_analyze_full.sh <R1> <R2> <accession> [threads]")
    logger.info("")
    logger.info("See vicast-analyze/README.md for full documentation.")
    sys.exit(0)


def validate_main():
    """Entry point for vicast-validate command."""
    from vicast.config import get_config

    logger = setup_logging("vicast.validate")

    logger.info("VICAST Configuration Validation")
    logger.info("================================")

    config = get_config()
    config.print_status()

    is_valid, errors = config.validate()

    if errors:
        logger.error("\nConfiguration Errors:")
        for error in errors:
            logger.error(f"  ✗ {error}")
        sys.exit(1)
    else:
        logger.info("\n✓ Configuration is valid")
        sys.exit(0)


if __name__ == "__main__":
    validate_main()
