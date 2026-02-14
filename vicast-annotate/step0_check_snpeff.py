#!/usr/bin/env python3
"""
STEP 0: Check if a viral genome is already available in SnpEff.

This is Pathway 1 of VICAST-annotate. If the genome is already in SnpEff's
database, no annotation work is needed - you can proceed directly to analysis.

Usage:
    # Check a single accession
    step0_check_snpeff.py NC_001477

    # Check multiple accessions
    step0_check_snpeff.py NC_001477 NC_001474.2 NC_045512.2

    # Search by keyword
    step0_check_snpeff.py --search dengue

    # List all viral genomes in SnpEff
    step0_check_snpeff.py --list-viral

    # Show details for a specific genome
    step0_check_snpeff.py NC_001477 --details
"""

import sys
import os
import argparse
import subprocess
import re
from pathlib import Path

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Import vicast configuration if available
try:
    from vicast.config import get_config
    CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False


def find_snpeff_jar():
    """Find the SnpEff JAR file."""
    # Check SNPEFF_JAR environment variable
    jar = os.environ.get('SNPEFF_JAR')
    if jar and os.path.isfile(jar):
        return jar

    # Check SNPEFF_HOME
    home = os.environ.get('SNPEFF_HOME')
    if home:
        jar = os.path.join(home, 'snpEff.jar')
        if os.path.isfile(jar):
            return jar

    # Check vicast config
    if CONFIG_AVAILABLE:
        config = get_config()
        jar = config.get('snpeff', {}).get('jar')
        if jar and os.path.isfile(jar):
            return jar

    # Check common locations
    locations = [
        Path.home() / 'snpEff' / 'snpEff.jar',
        Path('/opt/snpEff') / 'snpEff.jar',
        Path('/opt/conda/share/snpeff-5.4.0a-0') / 'snpEff.jar',
        Path(__file__).parent.parent / 'tools' / 'snpEff' / 'snpEff.jar',
    ]

    for loc in locations:
        if loc.exists():
            return str(loc)

    return None


def get_snpeff_databases(snpeff_jar):
    """Get list of all available SnpEff databases."""
    try:
        result = subprocess.run(
            ['java', '-jar', snpeff_jar, 'databases'],
            capture_output=True, text=True, timeout=120
        )
        return result.stdout
    except subprocess.TimeoutExpired:
        print("Error: SnpEff database query timed out", file=sys.stderr)
        return None
    except FileNotFoundError:
        print("Error: Java not found. Please install Java 11+", file=sys.stderr)
        return None


def check_accession(accession, db_text):
    """
    Check if an accession is in the SnpEff database list.

    Returns:
        list of matching lines (empty if not found)
    """
    matches = []
    for line in db_text.splitlines():
        # SnpEff databases output is tab-separated
        # Format: genome_name\torganism\tother_info
        if accession.lower() in line.lower():
            matches.append(line.strip())
    return matches


def search_databases(keyword, db_text):
    """Search SnpEff databases by keyword."""
    matches = []
    keyword_lower = keyword.lower()
    for line in db_text.splitlines():
        if keyword_lower in line.lower():
            matches.append(line.strip())
    return matches


def check_custom_database(accession, snpeff_data=None):
    """Check if accession exists as a custom-built database."""
    data_dirs = []

    # Check SNPEFF_DATA environment variable
    if snpeff_data:
        data_dirs.append(Path(snpeff_data))
    if os.environ.get('SNPEFF_DATA'):
        data_dirs.append(Path(os.environ['SNPEFF_DATA']))
    if os.environ.get('SNPEFF_DATA_CUSTOM'):
        data_dirs.append(Path(os.environ['SNPEFF_DATA_CUSTOM']))
    if os.environ.get('SNPEFF_DATA_BUILTIN'):
        data_dirs.append(Path(os.environ['SNPEFF_DATA_BUILTIN']))

    for data_dir in data_dirs:
        genome_dir = data_dir / accession
        if genome_dir.is_dir():
            has_predictor = (genome_dir / 'snpEffectPredictor.bin').exists()
            has_sequences = (genome_dir / 'sequences.fa').exists()
            has_genes = (genome_dir / 'genes.gff').exists()
            return {
                'found': True,
                'path': str(genome_dir),
                'has_predictor': has_predictor,
                'has_sequences': has_sequences,
                'has_genes': has_genes,
                'complete': has_predictor and has_sequences and has_genes
            }

    return {'found': False}


def get_database_details(accession, snpeff_jar):
    """Get detailed information about a SnpEff database using dump."""
    try:
        result = subprocess.run(
            ['java', '-jar', snpeff_jar, 'dump', accession],
            capture_output=True, text=True, timeout=60
        )
        if result.returncode == 0:
            return result.stdout
        return None
    except (subprocess.TimeoutExpired, FileNotFoundError):
        return None


def print_result(accession, builtin_matches, custom_info, details=False, snpeff_jar=None):
    """Print results for a single accession check."""
    found_anywhere = False

    # Check custom databases first (user-built take priority)
    if custom_info.get('found'):
        found_anywhere = True
        status = "complete" if custom_info['complete'] else "incomplete"
        print(f"  \u2705 {accession} found in CUSTOM database ({status})")
        print(f"     Path: {custom_info['path']}")
        if not custom_info['complete']:
            missing = []
            if not custom_info['has_predictor']:
                missing.append('snpEffectPredictor.bin')
            if not custom_info['has_sequences']:
                missing.append('sequences.fa')
            if not custom_info['has_genes']:
                missing.append('genes.gff')
            print(f"     Missing: {', '.join(missing)}")
            print(f"     Rebuild with: step2_add_to_snpeff.py {accession} <tsv_file>")

    # Check built-in databases
    if builtin_matches:
        found_anywhere = True
        print(f"  \u2705 {accession} found in SnpEff built-in database")
        for match in builtin_matches[:5]:  # Show max 5 matches
            print(f"     {match}")
        if len(builtin_matches) > 5:
            print(f"     ... and {len(builtin_matches) - 5} more matches")

    if not found_anywhere:
        print(f"  \u274c {accession} NOT found in SnpEff")
        print(f"     Recommended: Use Pathway 2 (standard annotation)")
        print(f"     Run: step1_parse_viral_genome.py {accession}")

    # Show details if requested
    if details and found_anywhere and snpeff_jar:
        db_name = accession
        if custom_info.get('found') and custom_info.get('complete'):
            detail_text = get_database_details(db_name, snpeff_jar)
            if detail_text:
                # Show first 30 lines of dump output
                lines = detail_text.splitlines()[:30]
                print(f"\n  Database details ({db_name}):")
                for line in lines:
                    print(f"     {line}")

    return found_anywhere


def list_viral_databases(db_text):
    """List all viral genomes in SnpEff database."""
    viral_keywords = [
        'virus', 'viral', 'phage', 'viridae', 'virales',
        'dengue', 'zika', 'ebola', 'influenza', 'hiv',
        'sars', 'corona', 'hepatitis', 'herpes', 'measles',
        'rabies', 'polio', 'rotavirus', 'norovirus', 'adenovirus',
    ]

    matches = set()
    for line in db_text.splitlines():
        line_lower = line.lower()
        for keyword in viral_keywords:
            if keyword in line_lower:
                matches.add(line.strip())
                break

    if matches:
        print(f"\nFound {len(matches)} viral genome(s) in SnpEff:\n")
        for match in sorted(matches):
            print(f"  {match}")
    else:
        print("\nNo viral genomes found in SnpEff built-in database.")
        print("Use VICAST pathways 2-4 to build custom databases.")

    return matches


def main():
    parser = argparse.ArgumentParser(
        description='Check if viral genomes are available in SnpEff database',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s NC_001477                    # Check Dengue-1
  %(prog)s NC_001477 NC_001474.2        # Check multiple
  %(prog)s --search dengue              # Search by keyword
  %(prog)s --list-viral                 # List all viral genomes
  %(prog)s NC_001477 --details          # Show database details

Pathways:
  If genome IS found     -> Pathway 1: Ready to use!
  If genome is NOT found -> Use Pathway 2, 3, or 4 to build database
        """
    )

    parser.add_argument(
        'accessions', nargs='*',
        help='Genome accession(s) to check (e.g., NC_001477, NC_001474.2)'
    )
    parser.add_argument(
        '--search', '-s',
        help='Search databases by keyword (e.g., "dengue", "influenza")'
    )
    parser.add_argument(
        '--list-viral', action='store_true',
        help='List all viral genomes in SnpEff'
    )
    parser.add_argument(
        '--details', '-d', action='store_true',
        help='Show detailed database information'
    )
    parser.add_argument(
        '--snpeff-jar',
        help='Path to snpEff.jar (default: auto-detect)'
    )
    parser.add_argument(
        '--snpeff-data',
        help='Path to SnpEff data directory (default: auto-detect)'
    )

    args = parser.parse_args()

    if not args.accessions and not args.search and not args.list_viral:
        parser.print_help()
        sys.exit(1)

    # Find SnpEff
    snpeff_jar = args.snpeff_jar or find_snpeff_jar()
    if not snpeff_jar:
        print("Error: Could not find snpEff.jar", file=sys.stderr)
        print("Set SNPEFF_JAR or SNPEFF_HOME environment variable,", file=sys.stderr)
        print("or use --snpeff-jar option", file=sys.stderr)
        sys.exit(1)

    print(f"Using SnpEff: {snpeff_jar}")
    print(f"Querying SnpEff database list...\n")

    # Get database list
    db_text = get_snpeff_databases(snpeff_jar)
    if db_text is None:
        sys.exit(1)

    # List viral databases
    if args.list_viral:
        list_viral_databases(db_text)
        return

    # Search by keyword
    if args.search:
        matches = search_databases(args.search, db_text)
        if matches:
            print(f"Found {len(matches)} match(es) for '{args.search}':\n")
            for match in matches:
                print(f"  {match}")
        else:
            print(f"No matches found for '{args.search}'")
            print("Try a broader search term, or use Pathway 2-4 to build a custom database")
        return

    # Check specific accessions
    if args.accessions:
        print(f"Checking {len(args.accessions)} accession(s):\n")
        all_found = True

        for accession in args.accessions:
            # Check built-in SnpEff databases
            builtin_matches = check_accession(accession, db_text)

            # Check custom databases
            custom_info = check_custom_database(accession, args.snpeff_data)

            # Print result
            found = print_result(
                accession, builtin_matches, custom_info,
                details=args.details, snpeff_jar=snpeff_jar
            )
            if not found:
                all_found = False
            print()

        # Summary
        if all_found:
            print("=" * 60)
            print("Pathway 1: All genomes are available in SnpEff!")
            print("You can proceed directly to variant annotation.")
            print("=" * 60)
        else:
            print("=" * 60)
            print("Some genomes need annotation. Recommended next steps:")
            print("  Pathway 2: step1_parse_viral_genome.py <accession>")
            print("  Pathway 3: step1_blastx_annotate.py <fasta>")
            print("  Pathway 4: For segmented viruses")
            print("=" * 60)

        sys.exit(0 if all_found else 1)


if __name__ == '__main__':
    main()
