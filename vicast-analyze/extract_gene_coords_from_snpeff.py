#!/usr/bin/env python3
"""
Extract gene coordinates from snpEff database and add to known_viruses.json

This script parses the snpEff genes.gbk file for a given accession and extracts:
- Gene coordinates (start, end positions)
- Gene names
- Genome length

Usage:
    python extract_gene_coords_from_snpeff.py <accession> <snpeff_dir> [--virus-name "Name"]

Example:
    python extract_gene_coords_from_snpeff.py AY532665.1 /ref/sahlab/software/snpEff --virus-name "West Nile virus NY99"
"""

import argparse
import json
import sys
from pathlib import Path
from Bio import SeqIO

def extract_gene_coords_from_genbank(gbk_file):
    """Extract gene coordinates from GenBank file"""
    gene_coords = {}
    genome_length = 0

    print(f"Parsing GenBank file: {gbk_file}")

    for record in SeqIO.parse(gbk_file, "genbank"):
        genome_length = len(record.seq)
        print(f"Genome length: {genome_length} bp")

        for feature in record.features:
            if feature.type == "CDS":
                # Get gene name
                gene_name = None
                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    gene_name = feature.qualifiers["product"][0]
                elif "locus_tag" in feature.qualifiers:
                    gene_name = feature.qualifiers["locus_tag"][0]

                if gene_name:
                    # Get coordinates
                    start = int(feature.location.start) + 1  # Convert to 1-based
                    end = int(feature.location.end)

                    gene_coords[gene_name] = [start, end]
                    print(f"  Found gene: {gene_name} [{start}..{end}]")

    return gene_coords, genome_length

def determine_virus_family(gene_coords):
    """Infer virus family from gene names"""
    gene_names = set(gene_coords.keys())

    # Flavivirus genes
    flavivirus_genes = {'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5', 'Env', 'C', 'prM', 'M'}
    if len(gene_names & flavivirus_genes) >= 5:
        return "flavivirus"

    # Alphavirus genes
    alphavirus_genes = {'nsP1', 'nsP2', 'nsP3', 'nsP4', 'E1', 'E2', 'E3', 'C'}
    if len(gene_names & alphavirus_genes) >= 4:
        return "alphavirus"

    return "unknown"

def classify_genes(gene_coords, virus_family):
    """Classify genes as structural or nonstructural"""
    structural_genes = []
    nonstructural_genes = []

    # Flavivirus classification
    if virus_family == "flavivirus":
        structural_keywords = ['C', 'prM', 'pr', 'M', 'Env', 'E', 'ancC']
        nonstructural_keywords = ['NS1', 'NS2', 'NS3', 'NS4', 'NS5', '2K']
    # Alphavirus classification
    elif virus_family == "alphavirus":
        structural_keywords = ['C', 'E1', 'E2', 'E3', '6K']
        nonstructural_keywords = ['nsP1', 'nsP2', 'nsP3', 'nsP4']
    else:
        # Unknown family - make educated guess
        structural_keywords = ['C', 'E', 'M', 'capsid', 'envelope', 'membrane']
        nonstructural_keywords = ['NS', 'nsP', 'P', 'L', 'polymerase', 'protease']

    for gene_name in gene_coords.keys():
        is_structural = any(keyword.lower() in gene_name.lower() for keyword in structural_keywords)
        is_nonstructural = any(keyword.lower() in gene_name.lower() for keyword in nonstructural_keywords)

        if is_structural and not is_nonstructural:
            structural_genes.append(gene_name)
        elif is_nonstructural:
            nonstructural_genes.append(gene_name)
        else:
            # Default to nonstructural if ambiguous
            nonstructural_genes.append(gene_name)

    return structural_genes, nonstructural_genes

def generate_colors(gene_coords, virus_family):
    """Generate color scheme for genes"""
    colors = {}

    # Default color schemes
    if virus_family == "flavivirus":
        default_colors = {
            "C": "#4575b4",
            "prM": "#74add1",
            "pr": "#20b2aa",
            "M": "#74add1",
            "Env": "#abd9e9",
            "NS1": "#fdae61",
            "NS2a": "#f46d43",
            "NS2b": "#d73027",
            "NS3": "#a50026",
            "NS4a": "#762a83",
            "2K": "#e6f5d0",
            "NS4b": "#9970ab",
            "NS5": "#c2a5cf",
            "ancC": "#2166ac"
        }
    elif virus_family == "alphavirus":
        default_colors = {
            "nsP1": "#fdae61",
            "nsP2": "#f46d43",
            "nsP3": "#d73027",
            "nsP4": "#a50026",
            "C": "#4575b4",
            "E3": "#5ba5c9",
            "E2": "#7fb3d3",
            "E1": "#abd9e9"
        }
    else:
        default_colors = {}

    # Assign colors to genes
    for gene_name in gene_coords.keys():
        if gene_name in default_colors:
            colors[gene_name] = default_colors[gene_name]
        else:
            # Generate a default color for unknown genes
            colors[gene_name] = "#999999"

    return colors

def add_to_known_viruses(accession, virus_name, gene_coords, genome_length, known_viruses_file):
    """Add virus entry to known_viruses.json"""
    # Read existing file
    if Path(known_viruses_file).exists():
        with open(known_viruses_file, 'r') as f:
            viruses = json.load(f)
    else:
        viruses = {}

    # Check if already exists
    if accession in viruses:
        print(f"\n⚠️  Warning: {accession} already exists in {known_viruses_file}")
        response = input("Overwrite? (y/n): ")
        if response.lower() != 'y':
            print("Aborted.")
            return False

    # Determine virus family
    virus_family = determine_virus_family(gene_coords)
    print(f"\nDetected virus family: {virus_family}")

    # Classify genes
    structural_genes, nonstructural_genes = classify_genes(gene_coords, virus_family)

    # Generate colors
    colors = generate_colors(gene_coords, virus_family)

    # Create entry
    viruses[accession] = {
        "name": virus_name,
        "family": virus_family,
        "genome_length": genome_length,
        "gene_coords": gene_coords,
        "colors": colors,
        "structural_genes": structural_genes,
        "nonstructural_genes": nonstructural_genes
    }

    # Write back to file
    with open(known_viruses_file, 'w') as f:
        json.dump(viruses, f, indent=2)

    print(f"\n✅ Successfully added {accession} to {known_viruses_file}")
    print(f"   Genome length: {genome_length} bp")
    print(f"   Total genes: {len(gene_coords)}")
    print(f"   Structural genes: {len(structural_genes)}")
    print(f"   Nonstructural genes: {len(nonstructural_genes)}")

    return True

def main():
    parser = argparse.ArgumentParser(
        description='Extract gene coordinates from snpEff database and add to known_viruses.json'
    )
    parser.add_argument('accession', help='Virus accession number (e.g., AY532665.1)')
    parser.add_argument('snpeff_dir', help='Path to snpEff directory (e.g., /ref/sahlab/software/snpEff)')
    parser.add_argument('--virus-name', help='Virus name (e.g., "West Nile virus NY99")', required=False)
    parser.add_argument('--output', help='Output known_viruses.json file', default='known_viruses.json')

    args = parser.parse_args()

    # Construct path to GenBank file in snpEff database
    gbk_file = Path(args.snpeff_dir) / "data" / args.accession / "genes.gbk"

    if not gbk_file.exists():
        print(f"❌ Error: GenBank file not found: {gbk_file}")
        print(f"\nMake sure the snpEff database for {args.accession} has been built.")
        print(f"Expected location: {gbk_file}")
        sys.exit(1)

    # Extract gene coordinates
    print(f"Extracting gene coordinates for {args.accession}...")
    gene_coords, genome_length = extract_gene_coords_from_genbank(gbk_file)

    if not gene_coords:
        print(f"❌ Error: No genes found in {gbk_file}")
        sys.exit(1)

    # Determine virus name
    if not args.virus_name:
        args.virus_name = input(f"\nEnter virus name for {args.accession}: ")

    # Add to known_viruses.json
    add_to_known_viruses(args.accession, args.virus_name, gene_coords, genome_length, args.output)

if __name__ == "__main__":
    main()
