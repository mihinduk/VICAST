#!/usr/bin/env python3
"""
Step 3 (Segmented): Add curated segmented virus to known_viruses.json

For segmented viruses (Influenza, etc.) created with vicast_annotate_segmented.py

NOTE: Consensus generation with individual proteins may not work correctly for
segmented genomes as it expects a single linear reference. This script primarily
documents the gene structure for reference.

Usage:
    python3 step3_add_to_known_viruses_segmented.py <genome_id> <tsv_file> \
        --name "Virus Name" --family "virus_family"

Example:
    python3 step3_add_to_known_viruses_segmented.py influenza_pr8 \
        influenza_pr8_annotations.tsv \
        --name "Influenza A virus (A/Puerto Rico/8/1934(H1N1))" \
        --family "orthomyxovirus"
"""

import sys
import os
import json
import argparse
from collections import OrderedDict

def parse_segmented_tsv(tsv_file):
    """Parse segmented virus TSV and extract gene coordinates per segment"""
    segments = OrderedDict()
    segment_accessions = {}

    with open(tsv_file, 'r') as f:
        header = f.readline().strip().split('\t')

        for line in f:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split('\t')
            # Pad fields to match header length
            while len(fields) < len(header):
                fields.append('')

            row = dict(zip(header, fields))

            # Only process CDS features that we're keeping
            if row['action'] not in ['KEEP', 'MODIFY']:
                continue

            feature_type = row['type']

            # Skip gene features (parent annotations)
            if feature_type == 'gene':
                continue

            if feature_type != 'CDS':
                continue

            segment = row.get('segment', '').strip()
            accession = row.get('accession', '').strip()
            gene_name = row.get('gene_name', '').strip()
            start = int(row['start'])
            end = int(row['end'])

            if not segment or not gene_name:
                continue

            # Track segment accessions
            if segment and accession:
                segment_accessions[segment] = accession

            # Store coordinates per segment
            if segment not in segments:
                segments[segment] = OrderedDict()

            # Use segment:gene as key for uniqueness
            key = f"{segment}:{gene_name}"
            segments[segment][gene_name] = [start, end]

    return segments, segment_accessions

def flatten_gene_coords(segments):
    """
    Flatten segmented coordinates into single dict with segment prefixes
    NOTE: This is for documentation only - consensus generation won't work correctly
    """
    flattened = OrderedDict()

    for segment, genes in segments.items():
        for gene_name, coords in genes.items():
            # Use segment:gene format to avoid name collisions
            key = f"{segment}:{gene_name}"
            flattened[key] = coords

    return flattened

def add_to_known_viruses(genome_id, virus_name, family, segments, segment_accessions,
                          known_viruses_path, force=False):
    """Add segmented virus configuration to known_viruses.json"""

    # Load existing configuration
    if os.path.exists(known_viruses_path):
        with open(known_viruses_path, 'r') as f:
            config = json.load(f, object_pairs_hook=OrderedDict)
    else:
        config = OrderedDict()

    # Check if already exists
    if genome_id in config and not force:
        print(f"Warning: {genome_id} already exists in known_viruses.json")
        print(f"Use --force to overwrite without prompting")
        return False

    # Flatten gene coordinates with segment prefixes
    gene_coords = flatten_gene_coords(segments)

    # Calculate total genes
    total_genes = sum(len(genes) for genes in segments.values())

    # Create new entry
    new_entry = OrderedDict([
        ("name", virus_name),
        ("family", family),
        ("genome_type", "segmented"),
        ("segments", list(segments.keys())),
        ("segment_accessions", segment_accessions),
        ("gene_coords", gene_coords),
        ("note", "Segmented genome - individual protein consensus may not work correctly")
    ])

    # Add to configuration
    config[genome_id] = new_entry

    # Write back to file
    with open(known_viruses_path, 'w') as f:
        json.dump(config, f, indent=2)

    print(f"\n✓ Added {genome_id} to known_viruses.json")
    print(f"  Name: {virus_name}")
    print(f"  Family: {family}")
    print(f"  Segments: {len(segments)}")
    print(f"  Total genes: {total_genes}")

    for segment, genes in segments.items():
        print(f"  {segment}: {', '.join(genes.keys())}")

    print(f"\nNOTE: This is a segmented genome. Individual protein consensus generation")
    print(f"      may not work correctly as it expects a single linear reference.")
    print(f"      The configuration is primarily for documentation purposes.")

    return True

def main():
    parser = argparse.ArgumentParser(
        description='Add segmented virus to known_viruses.json',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python3 step3_add_to_known_viruses_segmented.py influenza_pr8 \\
        influenza_pr8_annotations.tsv \\
        --name "Influenza A virus (A/Puerto Rico/8/1934(H1N1))" \\
        --family "orthomyxovirus"
        """
    )

    parser.add_argument('genome_id', help='Genome ID (e.g., influenza_pr8)')
    parser.add_argument('tsv_file', help='Curated TSV file from vicast_annotate_segmented.py')
    parser.add_argument('--name', required=True, help='Virus name')
    parser.add_argument('--family', default='unknown', help='Virus family')
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
        script_dir = os.path.dirname(os.path.abspath(__file__))
        known_viruses_path = os.path.join(script_dir, '..', 'vicast-analyze', 'known_viruses.json')

    if not os.path.exists(known_viruses_path):
        print(f"Error: known_viruses.json not found at: {known_viruses_path}")
        print("Specify path with --known-viruses")
        sys.exit(1)

    print(f"Reading segmented genome from: {args.tsv_file}")
    segments, segment_accessions = parse_segmented_tsv(args.tsv_file)

    total_genes = sum(len(genes) for genes in segments.values())
    print(f"Found {len(segments)} segments with {total_genes} total genes")

    if len(segments) == 0:
        print("Error: No segments found in TSV file")
        sys.exit(1)

    # Add to known_viruses.json
    success = add_to_known_viruses(
        args.genome_id,
        args.name,
        args.family,
        segments,
        segment_accessions,
        known_viruses_path,
        force=args.force
    )

    if success:
        print(f"\n✓ SUCCESS! {args.genome_id} is now documented in known_viruses.json")
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
