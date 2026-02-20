#!/usr/bin/env python3
"""
Automatically add virus to known_viruses.json if not present

This script checks if a virus accession exists in known_viruses.json.
If not found, it automatically extracts gene coordinates from snpEff database.

Exit codes:
  0 - Virus already exists or successfully added
  1 - Virus not found and snpEff database doesn't exist
  2 - Failed to add virus (but non-fatal)

Usage:
    python auto_add_virus_to_json.py <accession> <snpeff_dir> <known_viruses_file>
"""

import sys
import json
from pathlib import Path
from Bio import SeqIO

def read_known_viruses(json_file):
    """Read existing known_viruses.json"""
    if Path(json_file).exists():
        with open(json_file, 'r') as f:
            return json.load(f)
    return {}

def extract_gene_coords_from_genbank(gbk_file):
    """Extract gene coordinates from GenBank file"""
    gene_coords = {}
    genome_length = 0

    for record in SeqIO.parse(gbk_file, "genbank"):
        genome_length = len(record.seq)

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
                    # Get coordinates (1-based)
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    gene_coords[gene_name] = [start, end]

    return gene_coords, genome_length

def determine_virus_family(gene_coords):
    """Infer virus family from gene names"""
    gene_names = set(gene_coords.keys())

    # Flavivirus genes
    flavivirus_genes = {'NS1', 'NS2a', 'NS2b', 'NS3', 'NS4a', 'NS4b', 'NS5', 'Env', 'E', 'C', 'prM', 'M'}
    if len(gene_names & flavivirus_genes) >= 5:
        return "flavivirus"

    # Alphavirus genes
    alphavirus_genes = {'nsP1', 'nsP2', 'nsP3', 'nsP4', 'E1', 'E2', 'E3', 'C'}
    if len(gene_names & alphavirus_genes) >= 4:
        return "alphavirus"

    return "unknown"

def determine_genome_type(virus_family, gbk_file=None):
    """Infer genome type from virus family or GenBank metadata.

    Returns one of: 'ssRNA', 'dsRNA', 'ssDNA', 'dsDNA', 'unknown'.
    """
    # Family-based inference
    family_map = {
        'flavivirus': 'ssRNA',
        'wnv_lineage2': 'ssRNA',
        'alphavirus': 'ssRNA',
        'coronavirus': 'ssRNA',
        'paramyxovirus': 'ssRNA',
        'orthomyxovirus': 'ssRNA',
        'filovirus': 'ssRNA',
        'rhabdovirus': 'ssRNA',
        'picornavirus': 'ssRNA',
        'calicivirus': 'ssRNA',
        'togavirus': 'ssRNA',
        'reovirus': 'dsRNA',
        'rotavirus': 'dsRNA',
        'orbivirus': 'dsRNA',
        'parvovirus': 'ssDNA',
        'circovirus': 'ssDNA',
        'anellovirus': 'ssDNA',
        'herpesvirus': 'dsDNA',
        'adenovirus': 'dsDNA',
        'poxvirus': 'dsDNA',
    }

    family_lower = virus_family.lower() if virus_family else ''
    for key, gtype in family_map.items():
        if key in family_lower:
            return gtype

    # Try GenBank metadata if available
    if gbk_file and Path(gbk_file).exists():
        try:
            for record in SeqIO.parse(str(gbk_file), "genbank"):
                mol_type = ''
                for feature in record.features:
                    if feature.type == 'source':
                        mol_type = feature.qualifiers.get('mol_type', [''])[0].lower()
                        break
                if 'ssrna' in mol_type or 'ss-rna' in mol_type:
                    return 'ssRNA'
                elif 'dsrna' in mol_type or 'ds-rna' in mol_type:
                    return 'dsRNA'
                elif 'ssdna' in mol_type or 'ss-dna' in mol_type:
                    return 'ssDNA'
                elif 'dsdna' in mol_type or 'ds-dna' in mol_type:
                    return 'dsDNA'
                elif 'rna' in mol_type:
                    return 'ssRNA'
                elif 'dna' in mol_type:
                    return 'dsDNA'
        except Exception:
            pass

    return 'unknown'

def classify_genes(gene_coords, virus_family):
    """Classify genes as structural or nonstructural"""
    structural_genes = []
    nonstructural_genes = []

    if virus_family == "flavivirus":
        structural_keywords = ['C', 'prM', 'pr', 'M', 'Env', 'E', 'ancC']
        nonstructural_keywords = ['NS1', 'NS2', 'NS3', 'NS4', 'NS5', '2K']
    elif virus_family == "alphavirus":
        structural_keywords = ['C', 'E1', 'E2', 'E3', '6K']
        nonstructural_keywords = ['nsP1', 'nsP2', 'nsP3', 'nsP4']
    else:
        structural_keywords = ['C', 'E', 'M', 'capsid', 'envelope', 'membrane']
        nonstructural_keywords = ['NS', 'nsP', 'P', 'L', 'polymerase', 'protease']

    for gene_name in gene_coords.keys():
        is_structural = any(keyword in gene_name for keyword in structural_keywords)
        is_nonstructural = any(keyword in gene_name for keyword in nonstructural_keywords)

        if is_structural and not is_nonstructural:
            structural_genes.append(gene_name)
        elif is_nonstructural:
            nonstructural_genes.append(gene_name)
        else:
            nonstructural_genes.append(gene_name)

    return structural_genes, nonstructural_genes

def generate_colors(gene_coords, virus_family):
    """Generate color scheme for genes"""
    colors = {}

    if virus_family == "flavivirus":
        default_colors = {
            "C": "#4575b4", "prM": "#74add1", "pr": "#20b2aa", "M": "#74add1",
            "Env": "#abd9e9", "E": "#abd9e9", "NS1": "#fdae61", "NS2a": "#f46d43",
            "NS2b": "#d73027", "NS3": "#a50026", "NS4a": "#762a83", "2K": "#e6f5d0",
            "NS4b": "#9970ab", "NS5": "#c2a5cf", "ancC": "#2166ac"
        }
    elif virus_family == "alphavirus":
        default_colors = {
            "nsP1": "#fdae61", "nsP2": "#f46d43", "nsP3": "#d73027", "nsP4": "#a50026",
            "C": "#4575b4", "E3": "#5ba5c9", "E2": "#7fb3d3", "E1": "#abd9e9"
        }
    else:
        default_colors = {}

    for gene_name in gene_coords.keys():
        colors[gene_name] = default_colors.get(gene_name, "#999999")

    return colors

def get_virus_name_from_genbank(gbk_file):
    """Extract virus name from GenBank file"""
    for record in SeqIO.parse(gbk_file, "genbank"):
        if record.description:
            return record.description
        elif record.name:
            return record.name
    return "Unknown virus"

def add_virus_to_json(accession, snpeff_dir, json_file):
    """Add virus to known_viruses.json"""
    # Check if already exists
    viruses = read_known_viruses(json_file)

    if accession in viruses:
        print(f"✓ Virus {accession} already in {json_file}")
        return 0

    # Check if snpEff database exists
    gbk_file = Path(snpeff_dir) / "data" / accession / "genes.gbk"

    if not gbk_file.exists():
        print(f"⚠️  Virus {accession} not in {json_file}")
        print(f"⚠️  snpEff database not found at: {gbk_file}")
        print(f"⚠️  Will generate polyprotein only (no per-gene translation)")
        return 1

    print(f"⚠️  Virus {accession} not in {json_file}")
    print(f"✓ Found snpEff database at: {gbk_file}")
    print(f"⚙️  Auto-extracting gene coordinates...")

    try:
        # Extract gene coordinates
        gene_coords, genome_length = extract_gene_coords_from_genbank(gbk_file)

        if not gene_coords:
            print(f"❌ No genes found in {gbk_file}")
            return 2

        # Get virus name
        virus_name = get_virus_name_from_genbank(gbk_file)

        # Determine family
        virus_family = determine_virus_family(gene_coords)

        # Determine genome type
        genome_type = determine_genome_type(virus_family, gbk_file)

        # Classify genes
        structural_genes, nonstructural_genes = classify_genes(gene_coords, virus_family)

        # Generate colors
        colors = generate_colors(gene_coords, virus_family)

        # Add to dictionary
        viruses[accession] = {
            "name": virus_name,
            "family": virus_family,
            "genome_type": genome_type,
            "genome_length": genome_length,
            "gene_coords": gene_coords,
            "colors": colors,
            "structural_genes": structural_genes,
            "nonstructural_genes": nonstructural_genes
        }

        # Write to file
        with open(json_file, 'w') as f:
            json.dump(viruses, f, indent=2)

        print(f"✅ Successfully added {accession} to {json_file}")
        print(f"   Name: {virus_name}")
        print(f"   Family: {virus_family}")
        print(f"   Genome type: {genome_type}")
        print(f"   Genome length: {genome_length} bp")
        print(f"   Total genes: {len(gene_coords)}")

        return 0

    except Exception as e:
        print(f"❌ Failed to add {accession}: {e}")
        print(f"⚠️  Will generate polyprotein only")
        return 2

def main():
    if len(sys.argv) != 4:
        print("Usage: python auto_add_virus_to_json.py <accession> <snpeff_dir> <known_viruses_file>")
        sys.exit(1)

    accession = sys.argv[1]
    snpeff_dir = sys.argv[2]
    json_file = sys.argv[3]

    exit_code = add_virus_to_json(accession, snpeff_dir, json_file)
    sys.exit(exit_code)

if __name__ == "__main__":
    main()
