#!/usr/bin/env python3
"""
Add a genome from SnpEff database to known_viruses.json configuration.
Extracts gene coordinates from SnpEff GFF file.
"""

import sys
import json
import argparse
from pathlib import Path

def parse_gff(gff_file):
    """Parse GFF file and extract gene coordinates."""
    gene_coords = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            
            # Extract gene name from attributes
            gene_name = None
            for attr in attributes.split(';'):
                if attr.startswith('gene='):
                    gene_name = attr.split('=')[1]
                    break
                elif attr.startswith('ID=') and feature_type == 'mRNA':
                    # For mRNA features, use ID
                    gene_name = attr.split('=')[1]
                    break
            
            if gene_name and feature_type in ['mRNA', 'gene']:
                gene_coords[gene_name] = [start, end]
    
    return gene_coords

def get_genome_length(fasta_file):
    """Get genome length from FASTA file."""
    length = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                length += len(line.strip())
    return length

def add_genome_to_config(accession, name, family, snpeff_data_dir, config_file):
    """Add genome configuration to known_viruses.json."""
    
    # Locate files in SnpEff data directory
    genome_dir = Path(snpeff_data_dir) / accession
    gff_file = genome_dir / 'genes.gff'
    fasta_file = genome_dir / 'sequences.fa'
    
    if not gff_file.exists():
        print(f"Error: GFF file not found: {gff_file}")
        return False
    
    if not fasta_file.exists():
        print(f"Error: FASTA file not found: {fasta_file}")
        return False
    
    # Parse gene coordinates
    print(f"Reading gene coordinates from: {gff_file}")
    gene_coords = parse_gff(gff_file)
    
    if not gene_coords:
        print("Error: No gene coordinates found in GFF file")
        return False
    
    print(f"Found {len(gene_coords)} genes")
    
    # Get genome length
    genome_length = get_genome_length(fasta_file)
    print(f"Genome length: {genome_length} bp")
    
    # Read existing config
    config_path = Path(config_file)
    if config_path.exists():
        with open(config_path, 'r') as f:
            config = json.load(f)
    else:
        config = {}
    
    # Check if genome already exists
    if accession in config:
        print(f"\nWarning: {accession} already exists in config")
        response = input("Overwrite? (y/N): ")
        if response.lower() != 'y':
            print("Aborted")
            return False
    
    # Default colors for common viral genes
    default_colors = {
        'ancC': '#2166ac',
        'C': '#4575b4',
        'prM': '#74add1',
        'pr': '#20b2aa',
        'M': '#74add1',
        'E': '#abd9e9',
        'Env': '#abd9e9',
        'NS1': '#fdae61',
        'NS2A': '#f46d43',
        'NS2a': '#f46d43',
        'NS2B': '#d73027',
        'NS2b': '#d73027',
        'NS3': '#a50026',
        'NS4A': '#762a83',
        'NS4a': '#762a83',
        '2K': '#e6f5d0',
        'NS4B': '#9970ab',
        'NS4b': '#9970ab',
        'NS5': '#c2a5cf',
        'nsP1': '#fdae61',
        'nsP2': '#f46d43',
        'nsP3': '#d73027',
        'nsP4': '#a50026',
        'E1': '#abd9e9',
        'E2': '#7fb3d3',
        'E3': '#5ba5c9'
    }
    
    # Generate colors for genes
    colors = {}
    for gene in gene_coords.keys():
        if gene in default_colors:
            colors[gene] = default_colors[gene]
        else:
            # Default gray for unknown genes
            colors[gene] = '#808080'
    
    # Classify genes (basic heuristic)
    structural_genes = []
    nonstructural_genes = []
    
    for gene in gene_coords.keys():
        gene_lower = gene.lower()
        if any(s in gene_lower for s in ['c', 'prm', 'pr', 'm', 'e', 'env', 'ancc']):
            structural_genes.append(gene)
        elif any(s in gene_lower for s in ['ns', 'nsp']):
            nonstructural_genes.append(gene)
        else:
            # Default to nonstructural if unclear
            nonstructural_genes.append(gene)
    
    # Build config entry
    config[accession] = {
        'name': name,
        'family': family,
        'genome_length': genome_length,
        'gene_coords': gene_coords,
        'colors': colors,
        'structural_genes': structural_genes,
        'nonstructural_genes': nonstructural_genes
    }
    
    # Write updated config
    print(f"\nWriting updated config to: {config_path}")
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✓ Successfully added {accession} to config")
    print(f"  Name: {name}")
    print(f"  Family: {family}")
    print(f"  Genes: {len(gene_coords)}")
    print(f"  Structural: {', '.join(structural_genes)}")
    print(f"  Nonstructural: {', '.join(nonstructural_genes)}")
    
    return True

def main():
    parser = argparse.ArgumentParser(
        description='Add genome from SnpEff database to known_viruses.json',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add a Dengue virus genome
  python3 add_genome_to_config.py ON109597.1 "Dengue virus" flavivirus
  
  # Specify custom SnpEff directory
  python3 add_genome_to_config.py NC_001477.1 "Dengue virus type 1" flavivirus \
    --snpeff-data /path/to/snpEff/data
"""
    )
    
    parser.add_argument('accession', help='Genome accession (e.g., ON109597.1)')
    parser.add_argument('name', help='Virus name (e.g., "Dengue virus")')
    parser.add_argument('family', help='Virus family (e.g., flavivirus, alphavirus)')
    parser.add_argument('--snpeff-data', default='/ref/sahlab/software/snpEff/data',
                       help='SnpEff data directory (default: %(default)s)')
    parser.add_argument('--config', default=None,
                       help='Config file path (default: known_viruses.json in script dir)')
    
    args = parser.parse_args()
    
    # Default config file location
    if args.config is None:
        script_dir = Path(__file__).parent
        args.config = script_dir / 'known_viruses.json'
    
    success = add_genome_to_config(
        args.accession,
        args.name,
        args.family,
        args.snpeff_data,
        args.config
    )
    
    sys.exit(0 if success else 1)

if __name__ == '__main__':
    main()
