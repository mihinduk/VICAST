#!/usr/bin/env python3
"""
VICAST-annotate: Master controller for viral genome annotation pipeline.
Automatically determines and executes the appropriate annotation pathway.
"""

import sys
import os
import argparse
import subprocess

def run_step0(genome_id):
    """Run step0_check_snpeff.py to determine pathway."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    step0_script = os.path.join(script_dir, 'step0_check_snpeff.py')
    
    try:
        result = subprocess.run(
            ['python3', step0_script, genome_id],
            capture_output=True,
            text=True
        )
        # step0 returns pathway number as exit code
        return result.returncode, result.stdout
    except Exception as e:
        print(f"Error running step0: {e}")
        return None, None

def run_step1_standard(genome_id):
    """Run standard annotation pipeline (Pathway 2)."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    step1_script = os.path.join(script_dir, 'step1_parse_viral_genome.py')
    
    try:
        result = subprocess.run(
            ['python3', step1_script, genome_id],
            check=True
        )
        return f"{genome_id}_no_polyprotein.tsv"
    except subprocess.CalledProcessError as e:
        print(f"Error running step1: {e}")
        return None

def run_step1_vadr(genome_id):
    """Run VADR-enhanced annotation pipeline (Pathway 3)."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    step1_script = os.path.join(script_dir, 'step1_parse_viral_genome.py')
    
    try:
        result = subprocess.run(
            ['python3', step1_script, genome_id, '--use-vadr'],
            check=True
        )
        return f"{genome_id}_vadr_curated.tsv"
    except subprocess.CalledProcessError as e:
        print(f"Error running step1 with VADR: {e}")
        return None

def run_step1_blastx(fasta_file, blast_db='nr'):
    """Run BLASTx annotation pipeline (Pathway 4)."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    step1_script = os.path.join(script_dir, 'step1_blastx_annotate.py')
    
    try:
        result = subprocess.run(
            ['python3', step1_script, fasta_file, '--blast-db', blast_db],
            check=True
        )
        genome_id = os.path.basename(fasta_file).replace('.fasta', '').replace('.fa', '')
        return f"{genome_id}_blastx.tsv"
    except subprocess.CalledProcessError as e:
        print(f"Error running BLASTx annotation: {e}")
        return None

def run_step2(genome_id, tsv_file):
    """Run step2_add_to_snpeff.py."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    step2_script = os.path.join(script_dir, 'step2_add_to_snpeff.py')
    
    try:
        result = subprocess.run(
            ['python3', step2_script, genome_id, tsv_file],
            check=True
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error running step2: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='VICAST-annotate: Master controller for viral genome annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
VICAST-annotate automatically determines the best annotation pathway:

Pathway 1: Already in SnpEff
  → No action needed

Pathway 2: Well-annotated (NCBI)
  → Standard pipeline with manual curation

Pathway 3: VADR annotation
  → Enhanced validation with VADR

Pathway 4: BLASTx annotation
  → Homology-based annotation for poorly annotated genomes

Examples:
  # Automatic pathway selection
  python3 vicast_annotate.py NC_001477
  
  # Force specific pathway
  python3 vicast_annotate.py NC_001477 --pathway 3
  
  # Skip manual curation (auto-proceed to step2)
  python3 vicast_annotate.py NC_001477 --auto
  
  # Use custom BLAST database
  python3 vicast_annotate.py NC_XXXXXX --pathway 4 --blast-db /path/to/viral_db
        """
    )
    
    parser.add_argument('genome_id',
                       help='Genome ID (NCBI RefSeq accession) or FASTA file')
    parser.add_argument('--pathway', type=int, choices=[1, 2, 3, 4],
                       help='Force specific pathway (1-4)')
    parser.add_argument('--auto', action='store_true',
                       help='Auto-proceed without manual curation checkpoint')
    parser.add_argument('--blast-db', default='nr',
                       help='BLAST database for pathway 4 (default: nr)')
    parser.add_argument('--skip-step0', action='store_true',
                       help='Skip pathway detection (use --pathway)')
    
    args = parser.parse_args()
    
    # Determine pathway
    if args.pathway:
        pathway_num = args.pathway
        print(f"Using forced pathway {pathway_num}")
    elif args.skip_step0:
        print("Error: --skip-step0 requires --pathway")
        sys.exit(1)
    else:
        # Run step0 to determine pathway
        pathway_num, step0_output = run_step0(args.genome_id)
        print(step0_output)
        
        if pathway_num is None:
            print("Error: Could not determine annotation pathway")
            sys.exit(1)
    
    # Execute appropriate pathway
    if pathway_num == 1:
        print("\n" + "="*60)
        print("Pathway 1: Genome already in SnpEff")
        print("="*60)
        print(f"\nNo action needed! Use the genome directly:")
        print(f"  snpeff {args.genome_id} your_variants.vcf > annotated.vcf")
        sys.exit(0)
    
    elif pathway_num == 2:
        print("\n" + "="*60)
        print("Pathway 2: Well-annotated (NCBI)")
        print("="*60)
        
        # Run step1 standard
        tsv_file = run_step1_standard(args.genome_id)
        
        if not tsv_file:
            sys.exit(1)
        
        # Manual curation checkpoint
        if not args.auto:
            print("\n" + "="*60)
            print("MANUAL CURATION CHECKPOINT")
            print("="*60)
            print(f"\nPlease review and edit: {tsv_file}")
            input("\nPress ENTER when ready to proceed to step2...")
        
        # Run step2
        success = run_step2(args.genome_id, tsv_file)
        sys.exit(0 if success else 1)
    
    elif pathway_num == 3:
        print("\n" + "="*60)
        print("Pathway 3: VADR-enhanced annotation")
        print("="*60)
        
        # Run step1 with VADR
        tsv_file = run_step1_vadr(args.genome_id)
        
        if not tsv_file:
            sys.exit(1)
        
        # Manual curation checkpoint
        if not args.auto:
            print("\n" + "="*60)
            print("MANUAL CURATION CHECKPOINT")
            print("="*60)
            print(f"\nPlease review and edit: {tsv_file}")
            print("Also review VADR validation results")
            input("\nPress ENTER when ready to proceed to step2...")
        
        # Run step2
        success = run_step2(args.genome_id, tsv_file)
        sys.exit(0 if success else 1)
    
    elif pathway_num == 4:
        print("\n" + "="*60)
        print("Pathway 4: BLASTx-based annotation")
        print("="*60)
        
        # Determine FASTA file
        if os.path.exists(args.genome_id) and (args.genome_id.endswith('.fasta') or args.genome_id.endswith('.fa')):
            fasta_file = args.genome_id
            genome_id = os.path.basename(fasta_file).replace('.fasta', '').replace('.fa', '')
        else:
            fasta_file = f"{args.genome_id}.fasta"
            genome_id = args.genome_id
            if not os.path.exists(fasta_file):
                print(f"Error: FASTA file not found: {fasta_file}")
                sys.exit(1)
        
        # Run step1 BLASTx
        tsv_file = run_step1_blastx(fasta_file, args.blast_db)
        
        if not tsv_file:
            sys.exit(1)
        
        # Manual curation checkpoint
        if not args.auto:
            print("\n" + "="*60)
            print("MANUAL CURATION CHECKPOINT")
            print("="*60)
            print(f"\nPlease review and edit: {tsv_file}")
            print("BLASTx results may need significant curation")
            input("\nPress ENTER when ready to proceed to step2...")
        
        # Run step2
        success = run_step2(genome_id, tsv_file)
        sys.exit(0 if success else 1)
    
    else:
        print(f"Error: Invalid pathway number: {pathway_num}")
        sys.exit(1)

if __name__ == '__main__':
    main()
