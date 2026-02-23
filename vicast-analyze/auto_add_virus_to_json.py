#!/usr/bin/env python3
"""
Automatically add virus to known_viruses.json if not present

This script checks if a virus accession exists in known_viruses.json.
If not found, it automatically extracts gene coordinates from snpEff database.
Supports both single-segment and multi-segment (e.g. influenza) viruses.

Exit codes:
  0 - Virus already exists or successfully added
  1 - Virus not found and snpEff database doesn't exist
  2 - Failed to add virus (but non-fatal)

Usage:
    python auto_add_virus_to_json.py <accession> <snpeff_dir> <known_viruses_file>
"""

import sys
import os
import json
from pathlib import Path
from Bio import SeqIO


def read_known_viruses(json_file):
    """Read existing known_viruses.json"""
    if Path(json_file).exists():
        with open(json_file, 'r') as f:
            return json.load(f)
    return {}


def _find_database_dir(accession, snpeff_dir):
    """Search multiple snpEff data directories for the accession database.

    Returns (gbk_path_or_None, gff_path_or_None) for the first directory found.
    """
    search_dirs = [
        Path(snpeff_dir) / "data" / accession,
    ]

    # Check custom snpEff data directories
    snpeff_data_env = os.environ.get('SNPEFF_DATA', '')
    if snpeff_data_env:
        search_dirs.append(Path(snpeff_data_env) / accession)

    # Docker custom mount point
    custom_mount = Path('/opt/vicast/snpeff_data_custom') / accession
    if custom_mount.exists():
        search_dirs.append(custom_mount)

    for db_dir in search_dirs:
        gbk = db_dir / "genes.gbk"
        gff = db_dir / "genes.gff"
        if gbk.exists() or gff.exists():
            return (gbk if gbk.exists() else None, gff if gff.exists() else None)

    return None, None


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


def extract_gene_coords_from_gff(gff_file):
    """Extract gene coordinates from GFF3 file.

    For segmented viruses (multiple seqids in column 1), returns:
      (gene_coords_by_segment, total_length, True)
    For single-sequence viruses:
      (gene_coords, total_length, False)
    """
    # First pass: collect all CDS features grouped by seqid
    cds_by_seqid = {}
    all_seqids = set()

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.strip() == '':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue

            seqid = parts[0]
            feature_type = parts[2]
            all_seqids.add(seqid)

            if feature_type != 'CDS':
                continue

            start = int(parts[3])  # GFF is 1-based
            end = int(parts[4])

            # Parse attributes for gene name
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attrs[key.strip()] = val.strip()

            gene_name = attrs.get('gene', attrs.get('ID', attrs.get('Name', '')))
            if not gene_name:
                continue

            if seqid not in cds_by_seqid:
                cds_by_seqid[seqid] = {}

            # Use the widest span for this gene (handles spliced genes with multiple CDS lines)
            if gene_name in cds_by_seqid[seqid]:
                existing_start, existing_end = cds_by_seqid[seqid][gene_name]
                cds_by_seqid[seqid][gene_name] = [
                    min(existing_start, start),
                    max(existing_end, end)
                ]
            else:
                cds_by_seqid[seqid][gene_name] = [start, end]

    is_segmented = len(all_seqids) > 1

    if is_segmented:
        # Return per-segment gene coords
        return cds_by_seqid, 0, True
    else:
        # Single sequence — flatten
        all_coords = {}
        for seqid_genes in cds_by_seqid.values():
            all_coords.update(seqid_genes)
        return all_coords, 0, False


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


def determine_virus_family_segmented(gene_coords_by_segment):
    """Infer virus family from segmented gene names."""
    all_genes = set()
    for seg_genes in gene_coords_by_segment.values():
        all_genes.update(seg_genes.keys())

    # Orthomyxovirus (Influenza) genes
    influenza_genes = {'PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'M1', 'M2', 'NS1', 'NEP'}
    if len(all_genes & influenza_genes) >= 5:
        return "orthomyxovirus"

    # Bunyavirus genes
    bunyavirus_genes = {'L', 'M', 'S', 'N', 'Gn', 'Gc', 'RdRp'}
    if len(all_genes & bunyavirus_genes) >= 3:
        return "bunyavirus"

    # Reovirus/Rotavirus
    reovirus_genes = {'VP1', 'VP2', 'VP3', 'VP4', 'VP6', 'VP7', 'NSP1', 'NSP2', 'NSP3', 'NSP4'}
    if len(all_genes & reovirus_genes) >= 4:
        return "reovirus"

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
        'bunyavirus': 'ssRNA',
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
    elif virus_family == "orthomyxovirus":
        structural_keywords = ['HA', 'NA', 'M1', 'M2', 'NP']
        nonstructural_keywords = ['PB2', 'PB1', 'PA', 'NS1', 'NEP', 'PB1-F2', 'PA-X']
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


def classify_genes_segmented(gene_coords_by_segment, virus_family):
    """Classify genes across all segments as structural or nonstructural."""
    all_gene_coords = {}
    for seg_genes in gene_coords_by_segment.values():
        all_gene_coords.update(seg_genes)
    return classify_genes(all_gene_coords, virus_family)


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
    elif virus_family == "orthomyxovirus":
        default_colors = {
            "PB2": "#4575b4", "PB1": "#74add1", "PB1-F2": "#abd9e9",
            "PA": "#fdae61", "PA-X": "#f4a460",
            "HA": "#d73027", "NP": "#a50026", "NA": "#762a83",
            "M1": "#9970ab", "M2": "#c2a5cf",
            "NS1": "#e6f5d0", "NEP": "#1b7837"
        }
    else:
        default_colors = {}

    for gene_name in gene_coords.keys():
        colors[gene_name] = default_colors.get(gene_name, "#999999")

    return colors


def generate_colors_segmented(gene_coords_by_segment, virus_family):
    """Generate colors for all genes across segments."""
    all_gene_coords = {}
    for seg_genes in gene_coords_by_segment.values():
        all_gene_coords.update(seg_genes)
    return generate_colors(all_gene_coords, virus_family)


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
        print(f"Virus {accession} already in {json_file}")
        return 0

    # Search for database files in multiple locations
    gbk_file, gff_file = _find_database_dir(accession, snpeff_dir)

    if not gbk_file and not gff_file:
        print(f"Virus {accession} not in {json_file}")
        print(f"snpEff database not found in any search path")
        print(f"Will generate polyprotein only (no per-gene translation)")
        return 1

    source = "GenBank" if gbk_file else "GFF3"
    source_path = gbk_file or gff_file
    print(f"Virus {accession} not in {json_file}")
    print(f"Found snpEff database ({source}): {source_path}")
    print(f"Auto-extracting gene coordinates...")

    try:
        is_segmented = False

        if gff_file and gff_file.exists():
            # Prefer GFF — it handles segmented viruses natively
            gene_data, genome_length, is_segmented = extract_gene_coords_from_gff(str(gff_file))
        elif gbk_file and gbk_file.exists():
            gene_data, genome_length = extract_gene_coords_from_genbank(str(gbk_file))
        else:
            print(f"No readable database file for {accession}")
            return 2

        if not gene_data:
            print(f"No genes found in {source_path}")
            return 2

        # Get virus name from GenBank if available
        if gbk_file and gbk_file.exists():
            virus_name = get_virus_name_from_genbank(str(gbk_file))
        else:
            virus_name = accession

        if is_segmented:
            # gene_data is gene_coords_by_segment
            virus_family = determine_virus_family_segmented(gene_data)
            genome_type = determine_genome_type(virus_family,
                                                str(gbk_file) if gbk_file else None)
            structural_genes, nonstructural_genes = classify_genes_segmented(
                gene_data, virus_family)
            colors = generate_colors_segmented(gene_data, virus_family)

            n_genes = sum(len(v) for v in gene_data.values())

            viruses[accession] = {
                "name": virus_name,
                "family": virus_family,
                "genome_type": genome_type,
                "segmented": True,
                "gene_coords_by_segment": gene_data,
                "colors": colors,
                "structural_genes": structural_genes,
                "nonstructural_genes": nonstructural_genes
            }

            print(f"Successfully added {accession} to {json_file}")
            print(f"   Name: {virus_name}")
            print(f"   Family: {virus_family}")
            print(f"   Genome type: {genome_type}")
            print(f"   Segmented: {len(gene_data)} segments")
            print(f"   Total genes: {n_genes}")
        else:
            # gene_data is flat gene_coords
            virus_family = determine_virus_family(gene_data)
            genome_type = determine_genome_type(virus_family,
                                                str(gbk_file) if gbk_file else None)
            structural_genes, nonstructural_genes = classify_genes(
                gene_data, virus_family)
            colors = generate_colors(gene_data, virus_family)

            viruses[accession] = {
                "name": virus_name,
                "family": virus_family,
                "genome_type": genome_type,
                "genome_length": genome_length,
                "gene_coords": gene_data,
                "colors": colors,
                "structural_genes": structural_genes,
                "nonstructural_genes": nonstructural_genes
            }

            print(f"Successfully added {accession} to {json_file}")
            print(f"   Name: {virus_name}")
            print(f"   Family: {virus_family}")
            print(f"   Genome type: {genome_type}")
            print(f"   Genome length: {genome_length} bp")
            print(f"   Total genes: {len(gene_data)}")

        # Write to file
        with open(json_file, 'w') as f:
            json.dump(viruses, f, indent=2)

        return 0

    except Exception as e:
        print(f"Failed to add {accession}: {e}")
        print(f"Will generate polyprotein only")
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
