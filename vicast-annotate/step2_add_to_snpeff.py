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

# Import validation functions if available
try:
    from vicast_validation import validate_gff_for_snpeff, generate_curation_report
    VALIDATION_AVAILABLE = True
except ImportError:
    VALIDATION_AVAILABLE = False
    print("Note: vicast_validation module not found. Validation features disabled.")
    print("      To enable, ensure vicast_validation.py is in the same directory.")

def tsv_to_gff(tsv_file, output_gff):
    """
    Convert edited TSV back to GFF3 format.
    
    Args:
        tsv_file: Path to TSV file
        output_gff: Path to output GFF3 file
        
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
            if pd.notna(row.get('ID', '')) and row['ID']:
                attributes.append(f"ID={row['ID']}")
            else:
                # Generate ID if missing
                attributes.append(f"ID={row['type']}_{features_written + 1}")
            
            # Use gene_name column if available, otherwise fall back to gene column
            gene_value = row.get('gene_name', '') or row.get('gene', '')
            if pd.notna(gene_value) and gene_value:
                attributes.append(f"gene={gene_value}")
            
            if pd.notna(row.get('product', '')) and row['product']:
                attributes.append(f"product={row['product']}")
            if pd.notna(row.get('protein_id', '')) and row['protein_id']:
                attributes.append(f"protein_id={row['protein_id']}")
            
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
        response = input("Continue anyway? (y/N): ")
        if response.lower() != 'y':
            print("Aborted. Please check your TSV file.")
            return False, 0, {}
    
    # Show summary of features
    feature_counts = df['type'].value_counts()
    print("\nFeature summary:")
    for feat_type, count in feature_counts.items():
        print(f"  {feat_type}: {count}")
    
    return True, features_written, feature_counts.to_dict()

def add_genome_to_snpeff(genome_id, fasta_file, gff_file, snpeff_data_dir=None, validate=True):
    """
    Add a genome to snpEff database with optional validation.
    
    Args:
        genome_id: Genome identifier (e.g., NC_009942.1)
        fasta_file: Path to genome FASTA file
        gff_file: Path to GFF3 annotation file
        snpeff_data_dir: Path to snpEff data directory (optional)
        validate: Whether to validate GFF before adding (default: True)
    
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
            for error in errors[:10]:  # Show first 10 errors
                print(f"  - {error}")
            if len(errors) > 10:
                print(f"  ... and {len(errors) - 10} more errors")
            
            response = input("\nContinue despite errors? (y/N): ")
            if response.lower() != 'y':
                print("Aborted. Please fix errors and try again.")
                return False
        
        if warnings:
            print("\n[WARNING] Validation Warnings:")
            for warning in warnings[:10]:  # Show first 10 warnings
                print(f"  - {warning}")
            if len(warnings) > 10:
                print(f"  ... and {len(warnings) - 10} more warnings")
        
        if is_valid:
            print("\n[SUCCESS] GFF validation passed!")
        elif not errors:  # Only warnings
            print("\n[SUCCESS] GFF validation passed with warnings")
    
    # Determine snpEff data directory
    if not snpeff_data_dir:
        snpeff_home = os.environ.get('SNPEFF_HOME')
        if not snpeff_home:
            print("ERROR: SNPEFF_HOME not set. Please source the snpEff configuration:")
            print("  export SCRATCH_DIR=/path/to/your/scratch")
            print("  source /ref/sahlab/software/snpeff_configs/snpeff_current.sh")
            return False
        snpeff_data_dir = os.path.join(snpeff_home, 'data')
    
    # Create genome directory
    genome_dir = os.path.join(snpeff_data_dir, genome_id)
    os.makedirs(genome_dir, exist_ok=True)
    
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
        # Add genome to config
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
        # Check if snpeff command exists
        snpeff_jar = os.path.join(snpeff_home, 'snpEff.jar')
        java_home = os.environ.get('JAVA_HOME')
        
        if java_home and os.path.exists(snpeff_jar):
            # Use direct Java command
            java_cmd = os.path.join(java_home, 'bin', 'java')
            cmd = [java_cmd, '-jar', snpeff_jar, 'build', '-gff3', '-v', genome_id]
        else:
            # Try using snpeff command
            cmd = ['snpeff', 'build', '-gff3', '-v', genome_id]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"\n[SUCCESS] Successfully built snpEff database for {genome_id}")
            return True
        else:
            print(f"\n[ERROR] Error building database:")
            print(result.stderr)
            if "OutOfMemoryError" in result.stderr:
                print("\nTip: Try increasing memory allocation:")
                print("  export _JAVA_OPTIONS='-Xmx4g'")
            return False
            
    except Exception as e:
        print(f"\n[ERROR] Error running snpEff: {e}")
        print("\nMake sure snpEff environment is loaded:")
        print("  export SCRATCH_DIR=/path/to/your/scratch")
        print("  source /ref/sahlab/software/snpeff_configs/snpeff_current.sh")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='STEP 2: Add edited viral genome to snpEff database (with validation)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Prerequisites:
  1. You must have run step1_parse_viral_genome.py first
  2. You must have reviewed and edited the TSV file
  3. snpEff environment must be loaded:
     export SCRATCH_DIR=/path/to/your/scratch
     source /ref/sahlab/software/snpeff_configs/snpeff_current.sh

Examples:
  # Add genome using edited TSV
  python3 step2_add_to_snpeff.py NC_009942.1 NC_009942.1_no_polyprotein.tsv
  
  # Use custom GFF name
  python3 step2_add_to_snpeff.py NC_009942.1 my_edited.tsv --gff my_genome.gff3
  
  # Skip validation (not recommended)
  python3 step2_add_to_snpeff.py NC_009942.1 my_edited.tsv --no-validate
  
  # Generate curation report
  python3 step2_add_to_snpeff.py NC_009942.1 my_edited.tsv --report

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
    
    args = parser.parse_args()
    
    # Check TSV file exists
    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV file not found: {args.tsv_file}")
        sys.exit(1)
    
    # Determine FASTA file
    fasta_file = args.fasta if args.fasta else f"{args.genome_id}.fasta"
    if not os.path.exists(fasta_file):
        print(f"Error: FASTA file not found: {fasta_file}")
        print(f"\nMake sure you have run: download_ncbi_genome {args.genome_id}")
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
    
    # Convert TSV to GFF
    print("\n" + "-"*40)
    print("Converting edited TSV to GFF3...")
    print("-"*40)
    success, feature_count, feature_summary = tsv_to_gff(args.tsv_file, gff_file)
    
    if not success:
        print("Failed to convert TSV to GFF")
        sys.exit(1)
    
    # Generate curation report if requested
    if args.report and VALIDATION_AVAILABLE:
        print("\n" + "-"*40)
        print("Generating curation report...")
        print("-"*40)
        
        # Determine GenBank file
        gb_file = args.gb if args.gb else f"{args.genome_id}.gb"
        if not os.path.exists(gb_file):
            print(f"Warning: GenBank file not found for report: {gb_file}")
            print("Skipping detailed report generation")
        else:
            report_file = f"{args.genome_id}_curation_report.txt"
            generate_curation_report(gb_file, args.tsv_file, gff_file, report_file)
    
    # Add to snpEff
    print("\n" + "-"*40)
    print("Adding genome to snpEff database...")
    print("-"*40)
    success = add_genome_to_snpeff(args.genome_id, fasta_file, gff_file, 
                                   args.data_dir, validate=(not args.no_validate))
    
    if success:
        print("\n" + "="*60)
        print("SUCCESS!")
        print("="*60)
        print(f"\n[SUCCESS] Genome {args.genome_id} has been added to snpEff")
        print(f"\nYou can now use it for variant annotation:")
        print(f"  snpeff {args.genome_id} your_variants.vcf > annotated.vcf")
        print("\nTo test the database:")
        print(f"  snpeff dump {args.genome_id} | head")
    else:
        print("\n" + "="*60)
        print("FAILED")
        print("="*60)
        print(f"\n[ERROR] Failed to add {args.genome_id} to snpEff")
        print("\nTroubleshooting:")
        print("1. Check that snpEff environment is loaded")
        print("2. Check file permissions in snpEff data directory")
        print("3. Check the error messages above")

if __name__ == '__main__':
    main()