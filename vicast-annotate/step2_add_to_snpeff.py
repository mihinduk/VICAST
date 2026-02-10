#!/usr/bin/env python3
"""
STEP 2: Add edited viral genome to snpEff database (with validation).
This script takes the edited TSV file, validates it, and adds the genome to snpEff.
"""

import sys
import os
import argparse
import pandas as pd
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

# Import vicast configuration
try:
    from vicast.config import get_config, print_setup_instructions
    CONFIG_AVAILABLE = True
except ImportError:
    CONFIG_AVAILABLE = False

# Import validation functions if available
try:
    from vicast_validation import validate_gff_for_snpeff, generate_curation_report
    VALIDATION_AVAILABLE = True
except ImportError:
    VALIDATION_AVAILABLE = False
    print("Note: vicast_validation module not found. Validation features disabled.")
    print("      To enable, ensure vicast_validation.py is in the same directory.")


def tsv_to_gff(tsv_file, output_gff, force=False):
    """
    Convert edited TSV back to GFF3 format.

    Args:
        tsv_file: Path to TSV file
        output_gff: Path to output GFF3 file
        force: If True, continue without prompting on warnings

    Returns:
        tuple: (success, feature_count, feature_summary)
    """
    try:
        df = pd.read_csv(tsv_file, sep='\t')
    except Exception as e:
        print(f"Error reading TSV file: {e}")
        return False, 0, {}

    # Filter out rows marked for deletion
    if 'action' in df.columns:
        original_count = len(df)
        df = df[df['action'] != 'DELETE']
        deleted_count = original_count - len(df)
        if deleted_count > 0:
            print(f"  Removed {deleted_count} features marked for deletion")

    # Validate required columns
    required_cols = ['seqid', 'source', 'type', 'start', 'end', 'strand']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: TSV file missing required columns: {missing_cols}")
        return False, 0, {}

    # Backup existing GFF if it exists
    if os.path.exists(output_gff):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        backup_file = f"{output_gff}.backup.{timestamp}"
        shutil.copy(output_gff, backup_file)
        print(f"  Backed up existing GFF to: {backup_file}")

    with open(output_gff, 'w') as f:
        # Write header
        f.write("##gff-version 3\n")

        # Get sequence ID and length
        if not df.empty:
            seqid = df.iloc[0]['seqid']
            max_end = df['end'].max()
            f.write(f"##sequence-region {seqid} 1 {max_end}\n")

        features_written = 0
        for _, row in df.iterrows():
            # Skip if action is DELETE
            if 'action' in df.columns and row.get('action') == 'DELETE':
                continue

            # Build attributes
            attributes = []

            # Generate unique ID from gene_name, protein_id, or sequential counter
            gene_name_val = row.get('gene_name', '') or row.get('gene', '')
            protein_id_val = row.get('protein_id', '')

            if pd.notna(gene_name_val) and gene_name_val:
                # Use gene_name as ID (clean and readable)
                unique_id = gene_name_val
            elif pd.notna(protein_id_val) and protein_id_val:
                # Use protein_id if no gene_name
                unique_id = protein_id_val
            else:
                # Generate sequential ID only if needed
                unique_id = f"{row['type']}_{features_written + 1}"

            # For CDS features, create proper gene -> mRNA -> CDS hierarchy
            if row['type'] == 'CDS':
                # Write mRNA feature first (parent of CDS)
                mrna_id = unique_id
                mrna_attributes = [f"ID={mrna_id}"]

                if pd.notna(gene_name_val) and gene_name_val:
                    mrna_attributes.append(f"gene={gene_name_val}")
                if pd.notna(row.get('product', '')) and row['product']:
                    mrna_attributes.append(f"product={row['product']}")

                mrna_line = f"{row['seqid']}\t{row['source']}\tmRNA\t{row['start']}\t{row['end']}\t.\t{row['strand']}\t.\t{';'.join(mrna_attributes)}\n"
                f.write(mrna_line)
                features_written += 1

                # Now write CDS with Parent pointing to mRNA
                cds_id = f"{unique_id}.cds"
                attributes.append(f"ID={cds_id}")
                attributes.append(f"Parent={mrna_id}")
            else:
                # Non-CDS features (UTR, gene, etc.)
                attributes.append(f"ID={unique_id}")

            # Add common attributes
            if pd.notna(gene_name_val) and gene_name_val and row['type'] != 'CDS':
                attributes.append(f"gene={gene_name_val}")

            if pd.notna(row.get('product', '')) and row['product']:
                attributes.append(f"product={row['product']}")
            if pd.notna(protein_id_val) and protein_id_val:
                attributes.append(f"protein_id={protein_id_val}")

            # Add user notes if present
            if pd.notna(row.get('notes', '')) and row['notes']:
                attributes.append(f"Note={row['notes']}")

            # Write GFF line
            gff_line = f"{row['seqid']}\t{row['source']}\t{row['type']}\t{row['start']}\t{row['end']}\t.\t{row['strand']}\t.\t{';'.join(attributes)}\n"
            f.write(gff_line)
            features_written += 1

    print(f"Created GFF3 file from edited TSV: {output_gff}")
    print(f"  Total features written: {features_written}")

    # Check if no features were written
    if features_written == 0:
        print("\nWarning: No features were written to GFF!")
        if not force:
            response = input("Continue anyway? (y/N): ")
            if response.lower() != 'y':
                print("Aborted. Please check your TSV file.")
                return False, 0, {}
        else:
            print("  --force specified, continuing anyway...")

    # Show summary of features
    feature_counts = df['type'].value_counts()
    print("\nFeature summary:")
    for feat_type, count in feature_counts.items():
        print(f"  {feat_type}: {count}")

    return True, features_written, feature_counts.to_dict()


def generate_cds_protein_fasta(genome_fasta, tsv_file, genome_id, snpeff_data_dir):
    """
    Generate CDS and protein FASTA files for SnpEff validation.

    Args:
        genome_fasta: Path to genome FASTA file
        tsv_file: Path to TSV annotation file
        genome_id: Genome identifier
        snpeff_data_dir: SnpEff data directory

    Returns:
        tuple: (cds_file, protein_file) paths
    """
    # Load genome sequence
    genome_record = SeqIO.read(genome_fasta, "fasta")
    genome_seq = genome_record.seq

    print(f"  Extracting CDS sequences from genome ({len(genome_seq)} bp)")

    # Load TSV annotations
    df = pd.read_csv(tsv_file, sep='\t')

    # Filter for CDS features that are marked KEEP
    cds_df = df[(df['type'] == 'CDS') & (df.get('action', 'KEEP') != 'DELETE')]

    print(f"  Found {len(cds_df)} CDS features")

    cds_sequences = []
    protein_sequences = []

    feature_counter = 0

    for idx, row in cds_df.iterrows():
        start = int(row['start']) - 1  # Convert to 0-based
        end = int(row['end'])
        strand = row['strand']
        gene_name = row.get('gene_name', '') or row.get('gene', '')
        protein_id = row.get('protein_id', '')
        product = row.get('product', '')

        # Extract CDS sequence
        cds_seq = genome_seq[start:end]

        # Reverse complement if on minus strand
        if strand == '-':
            cds_seq = cds_seq.reverse_complement()

        # Generate sequence ID - must match mRNA ID from GFF3 (the parent of CDS)
        gene_name_str = str(gene_name) if pd.notna(gene_name) and gene_name else ''
        protein_id_str = str(protein_id) if pd.notna(protein_id) and protein_id else ''

        if gene_name_str:
            seq_id = gene_name_str
        elif protein_id_str:
            seq_id = protein_id_str
        else:
            feature_counter += 1
            seq_id = f"CDS_{feature_counter}"

        # Description
        product_str = str(product) if pd.notna(product) and product else ''
        description = product_str if product_str else "hypothetical protein"

        # Create CDS record
        cds_record = SeqIO.SeqRecord(
            cds_seq,
            id=seq_id,
            description=description
        )
        cds_sequences.append(cds_record)

        # Translate to protein
        try:
            protein_seq = cds_seq.translate(to_stop=True)

            if not str(protein_seq).startswith('M') and len(cds_seq) >= 3:
                protein_seq = Seq('M') + cds_seq[3:].translate(to_stop=True)

            protein_record = SeqIO.SeqRecord(
                protein_seq,
                id=seq_id,
                description=description
            )
            protein_sequences.append(protein_record)

        except Exception as e:
            print(f"    Warning: Could not translate {seq_id}: {e}")
            protein_record = SeqIO.SeqRecord(
                Seq(''),
                id=seq_id,
                description=f"{description} [translation failed]"
            )
            protein_sequences.append(protein_record)

    # Write to genome-specific directory in SnpEff data
    genome_dir = os.path.join(snpeff_data_dir, genome_id)
    os.makedirs(genome_dir, exist_ok=True)

    cds_file = os.path.join(genome_dir, 'cds.fa')
    protein_file = os.path.join(genome_dir, 'protein.fa')

    SeqIO.write(cds_sequences, cds_file, "fasta")
    SeqIO.write(protein_sequences, protein_file, "fasta")

    print(f"  Generated CDS FASTA: {cds_file} ({len(cds_sequences)} sequences)")
    print(f"  Generated protein FASTA: {protein_file} ({len(protein_sequences)} sequences)")

    return cds_file, protein_file


def add_genome_to_snpeff(genome_id, fasta_file, gff_file, snpeff_data_dir=None,
                         validate=True, update=False, force=False):
    """
    Add a genome to snpEff database with optional validation.

    Args:
        genome_id: Genome identifier (e.g., NC_009942.1)
        fasta_file: Path to genome FASTA file
        gff_file: Path to GFF3 annotation file
        snpeff_data_dir: Path to snpEff data directory (optional)
        validate: Whether to validate GFF before adding (default: True)
        update: Whether to overwrite existing genome entry (default: False)
        force: If True, continue without prompting on errors (default: False)

    Returns:
        bool: Success status
    """

    # Validate GFF if requested and available
    if validate and VALIDATION_AVAILABLE:
        print("\n" + "-"*40)
        print("Validating GFF3 file...")
        print("-"*40)

        is_valid, errors, warnings = validate_gff_for_snpeff(gff_file, fasta_file)

        if errors:
            print("\n[ERROR] Validation Errors Found:")
            for error in errors[:10]:
                print(f"  - {error}")
            if len(errors) > 10:
                print(f"  ... and {len(errors) - 10} more errors")

            if not force:
                response = input("\nContinue despite errors? (y/N): ")
                if response.lower() != 'y':
                    print("Aborted. Please fix errors and try again.")
                    return False
            else:
                print("\n  --force specified, continuing despite errors...")

        if warnings:
            print("\n[WARNING] Validation Warnings:")
            for warning in warnings[:10]:
                print(f"  - {warning}")
            if len(warnings) > 10:
                print(f"  ... and {len(warnings) - 10} more warnings")

        if is_valid:
            print("\n[SUCCESS] GFF validation passed!")
        elif not errors:
            print("\n[SUCCESS] GFF validation passed with warnings")

    # Get configuration
    if CONFIG_AVAILABLE:
        config = get_config()
        snpeff_home = str(config.snpeff_home) if config.snpeff_home else os.environ.get('SNPEFF_HOME')
        java_home = str(config.java_home) if config.java_home else os.environ.get('JAVA_HOME')
        snpeff_jar = str(config.snpeff_jar) if config.snpeff_jar else None
    else:
        snpeff_home = os.environ.get('SNPEFF_HOME')
        java_home = os.environ.get('JAVA_HOME')
        snpeff_jar = os.environ.get('SNPEFF_JAR')

    # Pre-flight environment check
    if not snpeff_home:
        print("\n" + "="*60)
        print("ERROR: SnpEff environment not configured")
        print("="*60)
        print("\nThe SNPEFF_HOME environment variable is not set.")
        print("\nTo configure VICAST, set the following environment variables:")
        print("\n  export SNPEFF_HOME=/path/to/snpEff")
        print("  export SNPEFF_JAR=${SNPEFF_HOME}/snpEff.jar")
        print("  export SNPEFF_DATA=${SNPEFF_HOME}/data")
        print("\nOr copy and configure the template:")
        print("  cp vicast_config.template.sh ~/.vicast/config.sh")
        print("  source ~/.vicast/config.sh")
        print("\nSee README.md for detailed installation instructions.")
        print("="*60)
        return False

    if not java_home:
        print("\n" + "="*60)
        print("WARNING: JAVA_HOME not set")
        print("="*60)
        print("\nJAVA_HOME is not set. Attempting to use system Java...")
        print("If the build fails, please set JAVA_HOME:")
        print("  export JAVA_HOME=/path/to/jdk-21")
        print("="*60 + "\n")

    if not snpeff_data_dir:
        snpeff_data_dir = os.path.join(snpeff_home, 'data')

    # Create genome directory
    genome_dir = os.path.join(snpeff_data_dir, genome_id)
    os.makedirs(genome_dir, exist_ok=True)

    # If updating, remove old snpEff database file to force clean rebuild
    if update:
        snpeff_bin = os.path.join(genome_dir, 'snpEffectPredictor.bin')
        if os.path.exists(snpeff_bin):
            print(f"  Removing old database file: {snpeff_bin}")
            os.remove(snpeff_bin)

    # Copy files to snpEff structure
    sequences_file = os.path.join(genome_dir, 'sequences.fa')
    genes_file = os.path.join(genome_dir, 'genes.gff')

    print(f"\nCopying files to snpEff data directory...")
    shutil.copy(fasta_file, sequences_file)
    shutil.copy(gff_file, genes_file)

    print(f"  Genome directory: {genome_dir}")
    print(f"  Sequences: {sequences_file}")
    print(f"  Annotations: {genes_file}")

    # Update snpEff.config
    config_file = os.path.join(os.path.dirname(snpeff_data_dir), 'snpEff.config')

    # Check if genome already in config
    genome_in_config = False
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            config_content = f.read()
            if f"{genome_id}.genome" in config_content:
                genome_in_config = True

    if not genome_in_config:
        print(f"\nAdding {genome_id} to snpEff configuration...")
        with open(config_file, 'a') as f:
            f.write(f"\n# {genome_id} - Added by VICAST-annotate\n")
            f.write(f"{genome_id}.genome : {genome_id}\n")
        print(f"  Added {genome_id} to snpEff.config")
    else:
        print(f"\n  {genome_id} already in snpEff.config")

    # Build the database
    print(f"\nBuilding snpEff database for {genome_id}...")
    print("This may take a minute...")

    try:
        # Determine snpEff JAR path
        if not snpeff_jar:
            snpeff_jar = os.path.join(snpeff_home, 'snpEff.jar')

        # Build base command
        if java_home and os.path.exists(snpeff_jar):
            java_cmd = os.path.join(java_home, 'bin', 'java')
            if not os.path.exists(java_cmd):
                java_cmd = 'java'  # Fall back to PATH
            cmd = [java_cmd, '-jar', snpeff_jar, 'build', '-gff3', '-v']
        else:
            cmd = ['snpeff', 'build', '-gff3', '-v']

        # Add validation skip flags if validation disabled
        if not validate:
            cmd.extend(['-noCheckProtein', '-noCheckCds'])

        cmd.append(genome_id)

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"\n[SUCCESS] Successfully built snpEff database for {genome_id}")
            return True
        else:
            print(f"\n[ERROR] Error building database:")
            print(result.stderr)

            stderr_lower = result.stderr.lower()

            if "outofmemoryerror" in stderr_lower:
                print("\nTip: Try increasing memory allocation:")
                print("  export _JAVA_OPTIONS='-Xmx4g'")

            if any(keyword in result.stderr for keyword in [
                "START codon errors",
                "STOP codon warnings",
                "CDS sequences comparison failed",
                "Database check failed",
                "Protein check",
                "codon"
            ]):
                print("\n" + "="*60)
                print("VALIDATION ERROR DETECTED")
                print("="*60)
                print("\nThis error typically occurs with:")
                print("  - Polyprotein genomes (mat_peptide features)")
                print("  - Genomes with non-standard START/STOP codons")
                print("  - CDS features that don't have proper codon boundaries")
                print("\nRECOMMENDED SOLUTION:")
                print("  Rerun with --no-validate flag to skip CDS/protein validation:")
                print(f"  python3 step2_add_to_snpeff.py {genome_id} <tsv_file> --no-validate")
                if update:
                    print(f"  python3 step2_add_to_snpeff.py {genome_id} <tsv_file> --update --no-validate")
                print("\nThe --no-validate flag is safe for viral genomes where:")
                print("  - Proteins are cleaved from polyproteins")
                print("  - Individual CDS features don't have START/STOP codons")
                print("  - Annotations come from mature peptide (mat_peptide) features")
                print("="*60)

            return False

    except Exception as e:
        print(f"\n[ERROR] Error running snpEff: {e}")

        error_str = str(e).lower()
        if "no such file" in error_str or "command not found" in error_str:
            print("\n" + "="*60)
            print("ENVIRONMENT ERROR: snpEff command not found")
            print("="*60)
            print("\nThe snpEff environment is not properly configured.")
            print("\nPlease set the following environment variables:")
            print("  export SNPEFF_HOME=/path/to/snpEff")
            print("  export SNPEFF_JAR=${SNPEFF_HOME}/snpEff.jar")
            print("  export JAVA_HOME=/path/to/jdk-21")
            print("\nSee README.md for detailed setup instructions.")
            print("="*60)
        else:
            print("\nMake sure snpEff environment is configured correctly.")
            print("See README.md for setup instructions.")
        return False


def main():
    parser = argparse.ArgumentParser(
        description='STEP 2: Add edited viral genome to snpEff database (with validation)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Prerequisites:
  1. You must have run step1_parse_viral_genome.py first
  2. You must have reviewed and edited the TSV file
  3. SnpEff must be installed and configured:
     export SNPEFF_HOME=/path/to/snpEff
     export SNPEFF_JAR=${SNPEFF_HOME}/snpEff.jar
     export SNPEFF_DATA=${SNPEFF_HOME}/data

Examples:
  # Add new genome using edited TSV
  python3 step2_add_to_snpeff.py NC_009942.1 NC_009942.1_no_polyprotein.tsv

  # Update existing genome (e.g., after adding missing proteins)
  python3 step2_add_to_snpeff.py NC_038433.1 NC_038433.1_no_polyprotein.tsv --update

  # Skip validation (for polyproteins with expected codon warnings)
  python3 step2_add_to_snpeff.py NC_027998.1 my_edited.tsv --no-validate

  # Non-interactive mode (for scripting/automation)
  python3 step2_add_to_snpeff.py NC_009942.1 my_edited.tsv --force

After successful addition:
  You can use the genome for annotation:
    snpeff NC_009942.1 your_variants.vcf > annotated.vcf
        """
    )

    parser.add_argument('genome_id',
                       help='Genome ID (e.g., NC_009942.1)')
    parser.add_argument('tsv_file',
                       help='Path to edited TSV file')
    parser.add_argument('--fasta',
                       help='Path to FASTA file (default: genome_id.fasta)')
    parser.add_argument('--gb',
                       help='Path to original GenBank file (for report)')
    parser.add_argument('--gff',
                       help='Output GFF file name (default: genome_id_final.gff3)')
    parser.add_argument('--data-dir',
                       help='snpEff data directory (default: auto-detect from SNPEFF_HOME)')
    parser.add_argument('--no-validate', action='store_true',
                       help='Skip GFF validation')
    parser.add_argument('--report', action='store_true',
                       help='Generate curation report')
    parser.add_argument('--update', action='store_true',
                       help='Allow updating/overwriting existing genome in SnpEff database')
    parser.add_argument('--force', '-f', action='store_true',
                       help='Non-interactive mode: continue without prompting on warnings/errors')

    args = parser.parse_args()

    # Check TSV file exists
    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV file not found: {args.tsv_file}")
        sys.exit(1)

    # Determine FASTA file
    fasta_file = args.fasta if args.fasta else f"{args.genome_id}.fasta"
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        print(f"\nMake sure you have run: python3 step1_parse_viral_genome.py {args.genome_id}")
        sys.exit(1)

    # Determine output GFF name
    gff_file = args.gff if args.gff else f"{args.genome_id}_final.gff3"

    print("="*60)
    print("STEP 2: Add Edited Viral Genome to snpEff")
    print("="*60)

    print(f"\nGenome ID: {args.genome_id}")
    print(f"Input TSV: {args.tsv_file}")
    print(f"FASTA file: {fasta_file}")
    print(f"Output GFF: {gff_file}")
    print(f"Validation: {'Enabled' if not args.no_validate else 'Disabled'}")
    print(f"Force mode: {'Yes' if args.force else 'No'}")

    # Get configuration for snpEff data directory
    snpeff_data_dir = args.data_dir
    if not snpeff_data_dir:
        if CONFIG_AVAILABLE:
            config = get_config()
            if config.snpeff_data:
                snpeff_data_dir = str(config.snpeff_data)

        if not snpeff_data_dir:
            snpeff_home = os.environ.get('SNPEFF_HOME')
            if snpeff_home:
                snpeff_data_dir = os.path.join(snpeff_home, 'data')
            else:
                snpeff_data_dir = 'data'

    genome_dir = os.path.join(snpeff_data_dir, args.genome_id)

    if os.path.exists(genome_dir):
        if not args.update:
            print("\n" + "="*60)
            print("ERROR: Genome already exists in SnpEff database")
            print("="*60)
            print(f"\nGenome directory exists: {genome_dir}")
            print("\nTo prevent accidental overwrites, you must explicitly use --update flag")
            print("to update an existing genome.")
            print("\nOptions:")
            print(f"  1. Use --update to overwrite existing genome:")
            print(f"     python3 step2_add_to_snpeff.py {args.genome_id} {args.tsv_file} --update")
            print(f"\n  2. Remove existing genome first:")
            print(f"     rm -rf {genome_dir}")
            print(f"\n  3. Use a different genome ID")
            print("\n" + "="*60)
            sys.exit(1)
        else:
            print("\n" + "="*60)
            print("WARNING: Updating existing genome in SnpEff")
            print("="*60)
            print(f"\nExisting genome directory: {genome_dir}")
            print("This will overwrite existing files.")
            print("\nProceeding with update...")
            print("="*60)

    # Convert TSV to GFF
    print("\n" + "-"*40)
    print("Converting edited TSV to GFF3...")
    print("-"*40)
    success, feature_count, feature_summary = tsv_to_gff(args.tsv_file, gff_file, force=args.force)

    if not success:
        print("Failed to convert TSV to GFF")
        sys.exit(1)

    # Generate curation report if requested
    if args.report and VALIDATION_AVAILABLE:
        print("\n" + "-"*40)
        print("Generating curation report...")
        print("-"*40)

        gb_file = args.gb if args.gb else f"{args.genome_id}.gb"
        if not os.path.exists(gb_file):
            print(f"Warning: GenBank file not found for report: {gb_file}")
            print("Skipping detailed report generation")
        else:
            report_file = f"{args.genome_id}_curation_report.txt"
            generate_curation_report(gb_file, args.tsv_file, gff_file, report_file)

    # Generate CDS and protein FASTA files
    print("\n" + "-"*40)
    print("Generating CDS and protein FASTA files...")
    print("-"*40)

    try:
        cds_file, protein_file = generate_cds_protein_fasta(
            fasta_file, args.tsv_file, args.genome_id, snpeff_data_dir
        )
    except Exception as e:
        print(f"  Warning: Could not generate CDS/protein files: {e}")
        print("  Continuing without validation files...")

    # Add to snpEff
    print("\n" + "-"*40)
    print("Adding genome to snpEff database...")
    print("-"*40)
    success = add_genome_to_snpeff(args.genome_id, fasta_file, gff_file,
                                   snpeff_data_dir, validate=(not args.no_validate),
                                   update=args.update, force=args.force)

    if success:
        print("\n" + "="*60)
        print("SUCCESS!")
        print("="*60)
        print(f"\n[SUCCESS] Genome {args.genome_id} has been added to snpEff")
        print(f"\nYou can now use it for variant annotation:")
        print(f"  snpeff {args.genome_id} your_variants.vcf > annotated.vcf")
        print("\nTo test the database:")
        print(f"  snpeff dump {args.genome_id} | head")
        print("\n" + "-"*60)
        print("NEXT STEPS: Run variant calling pipeline")
        print("-"*60)
        print("\nTo analyze sequencing data with this database:")
        print(f"  ./vicast-analyze/run_vicast_analyze_full.sh R1.fastq.gz R2.fastq.gz {args.genome_id} 4")
        print("\nSee vicast-analyze/README.md for full documentation")
        print("="*60)
    else:
        print("\n" + "="*60)
        print("FAILED")
        print("="*60)
        print(f"\n[ERROR] Failed to add {args.genome_id} to snpEff")
        print("\nTroubleshooting:")
        print("1. Check that SNPEFF_HOME is set correctly")
        print("2. Check file permissions in snpEff data directory")
        print("3. Check the error messages above")
        print("\nRun: python -c \"from vicast.config import get_config; get_config().print_status()\"")
        print("to verify your configuration.")


if __name__ == '__main__':
    main()
