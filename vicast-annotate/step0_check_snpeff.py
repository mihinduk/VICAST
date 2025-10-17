#!/usr/bin/env python3
"""
STEP 0: Check if genome is already in SnpEff database.
This script checks whether a genome needs to be added or is already available.
"""

import sys
import os
import argparse
import subprocess
from pathlib import Path

def check_snpeff_database(genome_id, snpeff_data_dir=None, snpeff_jar=None):
    """
    Check if a genome is already in the SnpEff database.
    
    Args:
        genome_id: Genome identifier (e.g., NC_009942.1)
        snpeff_data_dir: Path to SnpEff data directory (optional)
        snpeff_jar: Path to SnpEff JAR file (optional)
        
    Returns:
        tuple: (is_available, database_path, database_type)
            - is_available: True if genome is in SnpEff
            - database_path: Path to database directory if found
            - database_type: 'built-in' or 'custom' or None
    """
    
    # Get SnpEff paths from environment if not provided
    if not snpeff_data_dir:
        snpeff_home = os.environ.get('SNPEFF_HOME')
        snpeff_data = os.environ.get('SNPEFF_DATA')
        
        if snpeff_data:
            snpeff_data_dir = snpeff_data
        elif snpeff_home:
            snpeff_data_dir = os.path.join(snpeff_home, 'data')
        else:
            print("Warning: SNPEFF_HOME or SNPEFF_DATA not set")
            return False, None, None
    
    if not snpeff_jar:
        snpeff_home = os.environ.get('SNPEFF_HOME')
        if snpeff_home:
            snpeff_jar = os.path.join(snpeff_home, 'snpEff.jar')
    
    # Check 1: Look for genome directory in data folder
    genome_dir = os.path.join(snpeff_data_dir, genome_id)
    if os.path.exists(genome_dir):
        # Check if it has required files
        sequences_file = os.path.join(genome_dir, 'sequences.fa')
        genes_file = os.path.join(genome_dir, 'genes.gff')
        db_file = os.path.join(genome_dir, 'snpEffectPredictor.bin')
        
        if os.path.exists(sequences_file) or os.path.exists(genes_file) or os.path.exists(db_file):
            return True, genome_dir, 'custom'
    
    # Check 2: Query SnpEff databases list
    if snpeff_jar and os.path.exists(snpeff_jar):
        try:
            java_home = os.environ.get('JAVA_HOME')
            if java_home:
                java_cmd = os.path.join(java_home, 'bin', 'java')
                cmd = [java_cmd, '-jar', snpeff_jar, 'databases']
            else:
                cmd = ['snpeff', 'databases']
            
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            if result.returncode == 0:
                # Search for genome ID in output
                if genome_id in result.stdout:
                    return True, None, 'built-in'
        except subprocess.TimeoutExpired:
            print("Warning: SnpEff databases query timed out")
        except Exception as e:
            print(f"Warning: Could not query SnpEff databases: {e}")
    
    return False, None, None

def check_ncbi_annotation_quality(genome_id):
    """
    Check if NCBI has quality annotations for this genome.
    Downloads and inspects the GenBank file.
    
    Args:
        genome_id: NCBI RefSeq accession
        
    Returns:
        tuple: (has_annotations, quality_score, details)
            - has_annotations: True if CDS features exist
            - quality_score: 'good', 'poor', or 'none'
            - details: Dictionary with annotation statistics
    """
    gb_file = f"{genome_id}.gb"
    
    # Check if GenBank file exists, if not try to download
    if not os.path.exists(gb_file):
        print(f"Downloading GenBank file for {genome_id}...")
        try:
            from Bio import Entrez
            Entrez.email = "your_email@example.com"  # Set a default
            handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb", retmode="text")
            with open(gb_file, 'w') as f:
                f.write(handle.read())
            handle.close()
            print(f"  Downloaded: {gb_file}")
        except Exception as e:
            print(f"  Error downloading: {e}")
            return False, 'none', {}
    
    # Parse GenBank file
    try:
        from Bio import SeqIO
        details = {
            'cds_count': 0,
            'cds_with_product': 0,
            'cds_with_gene': 0,
            'has_polyprotein': False,
            'gene_count': 0,
        }
        
        for record in SeqIO.parse(gb_file, "genbank"):
            for feature in record.features:
                if feature.type == "CDS":
                    details['cds_count'] += 1
                    if 'product' in feature.qualifiers:
                        details['cds_with_product'] += 1
                        product = feature.qualifiers['product'][0].lower()
                        if 'polyprotein' in product:
                            details['has_polyprotein'] = True
                    if 'gene' in feature.qualifiers:
                        details['cds_with_gene'] += 1
                elif feature.type == "gene":
                    details['gene_count'] += 1
        
        # Determine quality
        if details['cds_count'] == 0:
            return False, 'none', details
        
        # Good quality: most CDS have products and genes
        annotation_rate = details['cds_with_product'] / details['cds_count'] if details['cds_count'] > 0 else 0
        
        if annotation_rate > 0.7 and details['cds_with_gene'] > 0:
            quality = 'good'
        elif annotation_rate > 0.3:
            quality = 'poor'
        else:
            quality = 'none'
        
        return True, quality, details
        
    except Exception as e:
        print(f"Error parsing GenBank file: {e}")
        return False, 'none', {}

def recommend_pathway(genome_id):
    """
    Recommend which annotation pathway to use for this genome.
    
    Args:
        genome_id: Genome identifier
        
    Returns:
        tuple: (pathway_number, pathway_name, reason, details)
    """
    print("="*60)
    print("STEP 0: Checking Genome Annotation Status")
    print("="*60)
    print(f"\nGenome ID: {genome_id}\n")
    
    # Check 1: Already in SnpEff?
    print("Checking SnpEff database...")
    is_available, db_path, db_type = check_snpeff_database(genome_id)
    
    if is_available:
        print(f"  ✓ Found in SnpEff ({db_type})")
        if db_path:
            print(f"    Database: {db_path}")
        return 1, "Already in SnpEff", "Genome is already available in SnpEff database", {
            'database_path': db_path,
            'database_type': db_type
        }
    else:
        print(f"  ✗ Not found in SnpEff")
    
    # Check 2: NCBI annotation quality
    print("\nChecking NCBI annotation quality...")
    has_annotations, quality, details = check_ncbi_annotation_quality(genome_id)
    
    print(f"  Annotation quality: {quality}")
    print(f"  CDS features: {details.get('cds_count', 0)}")
    print(f"  CDS with products: {details.get('cds_with_product', 0)}")
    print(f"  CDS with gene names: {details.get('cds_with_gene', 0)}")
    if details.get('has_polyprotein'):
        print(f"  ⚠ Contains polyprotein (will be handled)")
    
    # Recommend pathway
    if quality == 'good':
        return 2, "Well-annotated (NCBI)", \
               "NCBI has good quality annotations - use standard pipeline with manual curation", \
               details
    elif quality == 'poor':
        return 3, "VADR annotation", \
               "NCBI annotations are incomplete - VADR can improve them", \
               details
    else:
        return 4, "BLASTx annotation", \
               "No quality annotations available - use homology-based annotation", \
               details

def main():
    parser = argparse.ArgumentParser(
        description='STEP 0: Check genome annotation status and recommend pathway',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script determines which annotation pathway to use:

Pathway 1: Already in SnpEff
  → No action needed, use existing database

Pathway 2: Well-annotated (NCBI has good annotations)
  → Use step1_parse_viral_genome.py → step2_add_to_snpeff.py

Pathway 3: Can be annotated with VADR
  → Use step1_parse_viral_genome.py --use-vadr → step2_add_to_snpeff.py

Pathway 4: Requires BLASTx-based annotation
  → Use step1_blastx_annotate.py → step2_add_to_snpeff.py

Examples:
  python3 step0_check_snpeff.py NC_001477
  python3 step0_check_snpeff.py NC_009942.1
        """
    )
    
    parser.add_argument('genome_id', 
                       help='Genome ID (NCBI RefSeq accession, e.g., NC_009942.1)')
    parser.add_argument('--auto', action='store_true',
                       help='Automatically proceed with recommended pathway')
    
    args = parser.parse_args()
    
    # Get recommendation
    pathway_num, pathway_name, reason, details = recommend_pathway(args.genome_id)
    
    # Display recommendation
    print("\n" + "="*60)
    print("RECOMMENDATION")
    print("="*60)
    print(f"\nPathway {pathway_num}: {pathway_name}")
    print(f"\nReason: {reason}")
    
    # Provide next steps
    print("\n" + "-"*60)
    print("NEXT STEPS:")
    print("-"*60)
    
    if pathway_num == 1:
        print("\nGenome is already available in SnpEff!")
        print(f"You can use it directly:")
        print(f"  snpeff {args.genome_id} your_variants.vcf > annotated.vcf")
        
    elif pathway_num == 2:
        print("\nUse standard annotation pipeline:")
        print(f"  python3 step1_parse_viral_genome.py {args.genome_id}")
        print(f"  # Review and edit the TSV file")
        print(f"  python3 step2_add_to_snpeff.py {args.genome_id} {args.genome_id}_no_polyprotein.tsv")
        
    elif pathway_num == 3:
        print("\nUse VADR-enhanced annotation pipeline:")
        print(f"  python3 step1_parse_viral_genome.py {args.genome_id} --use-vadr")
        print(f"  # Review VADR validation results and edit TSV")
        print(f"  python3 step2_add_to_snpeff.py {args.genome_id} {args.genome_id}_vadr_curated.tsv")
        
    elif pathway_num == 4:
        print("\nUse BLASTx-based annotation pipeline:")
        print(f"  python3 step1_blastx_annotate.py {args.genome_id}")
        print(f"  # Review and edit the TSV file")
        print(f"  python3 step2_add_to_snpeff.py {args.genome_id} {args.genome_id}_blastx.tsv")
    
    print("\n" + "="*60)
    
    # Return pathway number as exit code for scripting
    sys.exit(pathway_num)

if __name__ == '__main__':
    main()
