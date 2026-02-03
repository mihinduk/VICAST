"""
VICAST Command-Line Interface

Entry points for vicast-annotate, vicast-analyze, and vicast-validate commands.
"""

import sys


def annotate_main():
    """Entry point for vicast-annotate command."""
    print("VICAST-Annotate")
    print("===============")
    print()
    print("For now, please use the Python scripts directly:")
    print()
    print("  # Check SnpEff database")
    print("  python vicast-annotate/step0_check_snpeff.py <genome_id>")
    print()
    print("  # Parse genome (Pathway 2)")
    print("  python vicast-annotate/step1_parse_viral_genome.py <genome_id>")
    print()
    print("  # Add to SnpEff")
    print("  python vicast-annotate/step2_add_to_snpeff.py <genome_id> <tsv_file>")
    print()
    print("See vicast-annotate/README.md for full documentation.")
    sys.exit(0)


def analyze_main():
    """Entry point for vicast-analyze command."""
    print("VICAST-Analyze")
    print("==============")
    print()
    print("For now, please use the shell scripts directly:")
    print()
    print("  # QC only (Steps 1-6)")
    print("  ./vicast-analyze/run_vicast_analyze_qc_only.sh <R1> <R2> <accession> [threads]")
    print()
    print("  # Annotation only (Steps 7-9)")
    print("  ./vicast-analyze/run_vicast_analyze_annotate_only.sh <R1> <R2> <accession>")
    print()
    print("  # Full pipeline")
    print("  ./vicast-analyze/run_vicast_analyze_full.sh <R1> <R2> <accession> [threads]")
    print()
    print("See vicast-analyze/README.md for full documentation.")
    sys.exit(0)


def validate_main():
    """Entry point for vicast-validate command."""
    from vicast.config import get_config

    print("VICAST Configuration Validation")
    print("================================")
    print()

    config = get_config()
    config.print_status()

    is_valid, errors = config.validate()

    if errors:
        print("\nConfiguration Errors:")
        for error in errors:
            print(f"  ✗ {error}")
        sys.exit(1)
    else:
        print("\n✓ Configuration is valid")
        sys.exit(0)


if __name__ == "__main__":
    validate_main()
