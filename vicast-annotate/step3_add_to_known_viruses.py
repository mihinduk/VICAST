#!/usr/bin/env python3
"""
Step 3: Add curated virus to known_viruses.json for vicast-analyze

After building SnpEff database with step2, this script extracts gene coordinates
from the curated TSV and adds them to known_viruses.json so vicast-analyze can
generate individual protein sequences from consensus genomes.

Usage:
    python3 step3_add_to_known_viruses.py <accession> <tsv_file> [--name "Virus Name"] [--family flavivirus]

Example:
    python3 step3_add_to_known_viruses.py NC_001474.2 NC_001474.2_no_polyprotein.tsv \
        --name "Dengue virus type 2" --family flavivirus
"""

import sys
import os
import json
import argparse
from collections import OrderedDict

# Default color palette for viral genes (matches existing configs)
DEFAULT_COLORS = {
    # Structural proteins (cool colors - blues/greens)
    'ancC': '#2166ac',
    'C': '#4575b4',
    'prM': '#74add1',
    'pr': '#20b2aa',
    'M': '#74add1',
    'E': '#abd9e9',
    'Env': '#abd9e9',
    # Nonstructural proteins (warm colors - oranges/reds/purples)
    'NS1': '#fdae61',
    'NS2a': '#f46d43',
    'NS2A': '#f46d43',
    'NS2b': '#d73027',
    'NS2B': '#d73027',
    'NS3': '#a50026',
    'NS4a': '#762a83',
    'NS4A': '#762a83',
    '2K': '#e6f5d0',
    'NS4b': '#9970ab',
    'NS4B': '#9970ab',
    'NS5': '#c2a5cf',
    # UTRs (gray)
    '5_UTR': '#f0f0f0',
    '3_UTR': '#f0f0f0',
    # Alphavirus proteins
    'nsP1': '#fdae61',
    'nsP2': '#f46d43',
    'nsP3': '#d73027',
    'nsP4': '#a50026',
    'E1': '#abd9e9',
    'E2': '#7fb3d3',
    'E3': '#5ba5c9',
}

# Common structural gene patterns (for auto-detection)
STRUCTURAL_PATTERNS = ['C', 'ancC', 'prM', 'pr', 'M', 'E', 'Env', 'E1', 'E2', 'E3']

def parse_tsv(tsv_file):
    """Parse curated TSV and extract gene coordinates"""
    genes = OrderedDict()
    genome_length = 0

    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')

        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            # Pad fields to match header length (some columns may be empty)
            while len(fields) < len(header):
                fields.append('')

            row = dict(zip(header, fields))

            # Only process CDS and UTR features that we're keeping
            if row['action'] not in ['KEEP', 'MODIFY']:
                continue

            feature_type = row['type']

            # Skip gene features (parent polyprotein annotations)
            if feature_type == 'gene':
                continue

            if feature_type not in ['CDS', "5'UTR", "3'UTR"]:
                continue

            gene_name = row.get('gene_name', '').strip()
            start = int(row['start'])
            end = int(row['end'])

            # Track genome length
            if end > genome_length:
                genome_length = end

            # Store coordinates
            if feature_type == 'CDS' and gene_name:
                # Use gene_name as key, store as list [start, end]
                genes[gene_name] = [start, end]
            elif feature_type == "5'UTR":
                genes['5_UTR'] = [start, end]
            elif feature_type == "3'UTR":
                genes['3_UTR'] = [start, end]

    return genes, genome_length

def classify_genes(genes, family):
    """Classify genes into structural and nonstructural"""
    structural = []
    nonstructural = []

    for gene_name in genes.keys():
        # Skip UTRs from classification
        if 'UTR' in gene_name:
            continue

        # Check if it matches structural patterns
        if gene_name in STRUCTURAL_PATTERNS:
            structural.append(gene_name)
        elif gene_name.startswith('ns') or gene_name.startswith('NS'):
            nonstructural.append(gene_name)
        elif family == 'flavivirus':
            # For flaviviruses: ancC, C, prM, pr, M, E are structural
            if gene_name in ['ancC', 'C', 'prM', 'pr', 'M', 'E']:
                structural.append(gene_name)
            else:
                nonstructural.append(gene_name)
        else:
            # Default: if it has "structural" in common name, classify as such
            # Otherwise nonstructural
            nonstructural.append(gene_name)

    return structural, nonstructural

def generate_colors(genes):
    """Generate color palette for genes"""
    colors = {}
    for gene_name in genes.keys():
        if gene_name in DEFAULT_COLORS:
            colors[gene_name] = DEFAULT_COLORS[gene_name]
        else:
            # Assign default color for unknown genes
            colors[gene_name] = '#999999'
    return colors

def add_to_known_viruses(accession, virus_name, family, genes, genome_length,
                          structural_genes, nonstructural_genes, known_viruses_path, force=False):
    """Add virus configuration to known_viruses.json"""

    # Load existing configuration
    if os.path.exists(known_viruses_path):
        with open(known_viruses_path, 'r') as f:
            config = json.load(f, object_pairs_hook=OrderedDict)
    else:
        config = OrderedDict()

    # Check if already exists
    if accession in config and not force:
        print(f"Warning: {accession} already exists in known_viruses.json")
        print(f"Use --force to overwrite without prompting")
        return False

    # Generate color palette
    colors = generate_colors(genes)

    # Create new entry
    new_entry = OrderedDict([
        ("name", virus_name),
        ("family", family),
        ("genome_length", genome_length),
        ("gene_coords", genes),
        ("colors", colors),
        ("structural_genes", structural_genes),
        ("nonstructural_genes", nonstructural_genes)
    ])

    # Add to configuration
    config[accession] = new_entry

    # Write back to file (with pretty formatting)
    with open(known_viruses_path, 'w') as f:
        json.dump(config, f, indent=2)

    print(f"\n✓ Added {accession} to known_viruses.json")
    print(f"  Name: {virus_name}")
    print(f"  Family: {family}")
    print(f"  Genome length: {genome_length}")
    print(f"  Genes: {len(genes)}")
    print(f"  Structural: {', '.join(structural_genes)}")
    print(f"  Nonstructural: {', '.join(nonstructural_genes)}")

    return True

def main():
    parser = argparse.ArgumentParser(
        description='Add curated virus to known_viruses.json for vicast-analyze',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python3 step3_add_to_known_viruses.py NC_001474.2 NC_001474.2_no_polyprotein.tsv \\
        --name "Dengue virus type 2" --family flavivirus
        """
    )

    parser.add_argument('accession', help='Genome accession ID (e.g., NC_001474.2)')
    parser.add_argument('tsv_file', help='Curated TSV file from step1')
    parser.add_argument('--name', required=True, help='Virus name (e.g., "Dengue virus type 2")')
    parser.add_argument('--family', default='unknown',
                        help='Virus family (e.g., flavivirus, alphavirus)')
    parser.add_argument('--known-viruses',
                        help='Path to known_viruses.json (default: ../vicast-analyze/known_viruses.json)')
    parser.add_argument('--force', action='store_true',
                        help='Overwrite existing entry without prompting')

    args = parser.parse_args()

    # Verify TSV file exists
    if not os.path.exists(args.tsv_file):
        print(f"Error: TSV file not found: {args.tsv_file}")
        sys.exit(1)

    # Determine path to known_viruses.json
    if args.known_viruses:
        known_viruses_path = args.known_viruses
    else:
        # Default: assume we're in vicast-annotate/, go to ../vicast-analyze/
        script_dir = os.path.dirname(os.path.abspath(__file__))
        known_viruses_path = os.path.join(script_dir, '..', 'vicast-analyze', 'known_viruses.json')

    if not os.path.exists(known_viruses_path):
        print(f"Error: known_viruses.json not found at: {known_viruses_path}")
        print("Specify path with --known-viruses")
        sys.exit(1)

    print(f"Reading gene coordinates from: {args.tsv_file}")
    genes, genome_length = parse_tsv(args.tsv_file)

    print(f"Found {len(genes)} genes/features")

    if len(genes) == 0:
        print("Error: No genes found in TSV file")
        print("Make sure TSV has KEEP or MODIFY actions for CDS features")
        sys.exit(1)

    # Classify genes
    structural, nonstructural = classify_genes(genes, args.family)

    # Add to known_viruses.json
    success = add_to_known_viruses(
        args.accession,
        args.name,
        args.family,
        genes,
        genome_length,
        structural,
        nonstructural,
        known_viruses_path,
        force=args.force
    )

    if success:
        print(f"\n✓ SUCCESS! {args.accession} is now configured for vicast-analyze")
        print(f"\nYou can now generate individual protein sequences with:")
        print(f"  python3 generate_filtered_consensus.py \\")
        print(f"    --vcf sample.tsv \\")
        print(f"    --reference genome.fasta \\")
        print(f"    --accession {args.accession} \\")
        print(f"    --output-prefix sample_consensus")
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
