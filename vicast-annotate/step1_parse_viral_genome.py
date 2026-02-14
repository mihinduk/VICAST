#!/usr/bin/env python3
"""
STEP 1: Download and parse viral genome annotations from NCBI.

This is Pathway 2 of VICAST-annotate. Downloads the genome from NCBI,
parses GenBank annotations, automatically skips polyprotein features,
and creates an editable TSV file for manual curation.

Usage:
    # Download and parse a single genome
    step1_parse_viral_genome.py NC_001477

    # Include version number
    step1_parse_viral_genome.py NC_001477.1

    # Keep polyprotein features
    step1_parse_viral_genome.py NC_001477 --keep-polyprotein

    # Specify output directory
    step1_parse_viral_genome.py NC_001477 --output-dir /data/annotations

    # Set custom NCBI email
    step1_parse_viral_genome.py NC_001477 --email your@email.edu
"""

import sys
import os
import argparse
import re
from pathlib import Path
from datetime import datetime

# Add src to path for vicast imports
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

try:
    from Bio import Entrez, SeqIO
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


def download_from_ncbi(accession, email, output_dir, verbose=False):
    """
    Download genome FASTA and GenBank files from NCBI.

    Args:
        accession: NCBI accession number (e.g., NC_001477 or NC_001477.1)
        email: Email address for NCBI Entrez
        output_dir: Directory to save files
        verbose: Print progress

    Returns:
        Tuple of (fasta_path, genbank_path) or (None, None) on failure
    """
    Entrez.email = email

    fasta_path = os.path.join(output_dir, f"{accession}.fasta")
    gb_path = os.path.join(output_dir, f"{accession}.gb")

    # Download FASTA
    if verbose:
        print(f"Downloading FASTA for {accession}...")
    try:
        handle = Entrez.efetch(db='nucleotide', id=accession, rettype='fasta', retmode='text')
        fasta_content = handle.read()
        handle.close()

        if not fasta_content.strip() or 'Error' in fasta_content[:100]:
            print(f"Error: Could not download FASTA for {accession}")
            print(f"Check that the accession exists at: https://www.ncbi.nlm.nih.gov/nuccore/{accession}")
            return None, None

        with open(fasta_path, 'w') as f:
            f.write(fasta_content)
        if verbose:
            print(f"  Saved: {fasta_path}")
    except Exception as e:
        print(f"Error downloading FASTA: {e}")
        return None, None

    # Download GenBank
    if verbose:
        print(f"Downloading GenBank for {accession}...")
    try:
        handle = Entrez.efetch(db='nucleotide', id=accession, rettype='gb', retmode='text')
        gb_content = handle.read()
        handle.close()

        with open(gb_path, 'w') as f:
            f.write(gb_content)
        if verbose:
            print(f"  Saved: {gb_path}")
    except Exception as e:
        print(f"Error downloading GenBank: {e}")
        return None, None

    return fasta_path, gb_path


def is_polyprotein(feature):
    """Check if a CDS feature is a polyprotein."""
    product = feature.qualifiers.get('product', [''])[0].lower()
    gene = feature.qualifiers.get('gene', [''])[0].lower()
    note = ' '.join(feature.qualifiers.get('note', [])).lower()

    polyprotein_terms = ['polyprotein', 'poly protein', 'polypeptide']

    for term in polyprotein_terms:
        if term in product or term in gene or term in note:
            return True

    return False


def extract_feature(feature, seqid, source='NCBI'):
    """
    Extract annotation info from a BioPython feature into a dict.

    Args:
        feature: BioPython SeqFeature
        seqid: Sequence ID
        source: Source annotation (default: NCBI)

    Returns:
        dict with feature information
    """
    gene_name = feature.qualifiers.get('gene', [''])[0]
    protein_id = feature.qualifiers.get('protein_id', [''])[0]
    note = '; '.join(feature.qualifiers.get('note', []))
    db_xref = '; '.join(feature.qualifiers.get('db_xref', []))

    # Set product and type based on feature type
    if feature.type == "5'UTR":
        product = feature.qualifiers.get('product', ["5' untranslated region"])[0]
        feat_type = "5UTR"
        if not gene_name:
            gene_name = "5UTR"
    elif feature.type == "3'UTR":
        product = feature.qualifiers.get('product', ["3' untranslated region"])[0]
        feat_type = "3UTR"
        if not gene_name:
            gene_name = "3UTR"
    elif feature.type == "mat_peptide":
        product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
        feat_type = "CDS"
        if not gene_name:
            gene_name = product.replace(' ', '_')[:20]
    else:
        product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
        feat_type = "CDS"

    # Handle location
    strand = '+' if feature.location.strand == 1 else '-'

    if hasattr(feature.location, 'parts') and len(feature.location.parts) > 1:
        start = min(int(p.start) + 1 for p in feature.location.parts)
        end = max(int(p.end) for p in feature.location.parts)
        location_type = 'join'
        parts_str = ','.join(
            f"{int(p.start)+1}..{int(p.end)}" for p in feature.location.parts
        )
    else:
        start = int(feature.location.start) + 1  # Convert to 1-based
        end = int(feature.location.end)
        location_type = 'simple'
        parts_str = f"{start}..{end}"

    # Generate gene name if missing
    if not gene_name:
        gene_name = product.replace(' ', '_')[:20]

    # Clean gene name for GFF3 compatibility
    gene_name_clean = re.sub(r'[^a-zA-Z0-9_-]', '_', gene_name)

    return {
        'seqid': seqid,
        'source': source,
        'type': feat_type,
        'start': start,
        'end': end,
        'strand': strand,
        'gene_name': gene_name_clean,
        'product': product,
        'protein_id': protein_id,
        'location_type': location_type,
        'location_detail': parts_str,
        'note': note,
        'db_xref': db_xref,
        'action': 'KEEP',
    }


def parse_genbank_to_tsv(gb_path, output_tsv, keep_polyprotein=False, verbose=False):
    """
    Parse GenBank file and create editable TSV with CDS features.

    When polyproteins are skipped (default), automatically extracts mat_peptide
    features as individual protein annotations. This handles viruses like Dengue,
    HCV, and coronaviruses where individual proteins are annotated as mature
    peptides within a polyprotein.

    Args:
        gb_path: Path to GenBank file
        output_tsv: Output TSV file path
        keep_polyprotein: If True, include polyprotein features
        verbose: Print progress

    Returns:
        Tuple of (feature_count, skipped_count, features_list)
    """
    features = []
    skipped = 0
    mat_peptide_count = 0

    # Feature types to extract
    feature_types = {"CDS", "5'UTR", "3'UTR", "mat_peptide"}

    for record in SeqIO.parse(gb_path, "genbank"):
        seqid = record.id
        organism = record.annotations.get('organism', 'unknown')

        if verbose:
            print(f"\nParsing: {seqid} ({organism})")
            print(f"  Sequence length: {len(record.seq)} bp")

        for feature in record.features:
            if feature.type not in feature_types:
                continue

            # Handle mat_peptide features
            if feature.type == "mat_peptide":
                if not keep_polyprotein:
                    # Extract mat_peptides as CDS when polyproteins are skipped
                    feat_dict = extract_feature(feature, seqid, source='NCBI_mat_peptide')
                    features.append(feat_dict)
                    mat_peptide_count += 1
                    if verbose:
                        print(f"  mat_peptide -> CDS: {feat_dict['gene_name']} ({feat_dict['product']})")
                continue

            # Check for polyprotein (CDS only)
            if feature.type == "CDS" and not keep_polyprotein and is_polyprotein(feature):
                product = feature.qualifiers.get('product', ['unknown'])[0]
                if verbose:
                    print(f"  Skipping polyprotein: {product}")
                skipped += 1
                continue

            # Extract standard feature
            feat_dict = extract_feature(feature, seqid)
            features.append(feat_dict)

    if mat_peptide_count > 0:
        print(f"  Extracted {mat_peptide_count} mature peptide(s) from polyprotein(s)")

    # Write TSV
    header = [
        'seqid', 'source', 'type', 'start', 'end',
        'strand', 'gene_name', 'product', 'protein_id',
        'location_type', 'location_detail', 'note', 'db_xref', 'action'
    ]

    with open(output_tsv, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for feat in features:
            row = [str(feat.get(h, '')) for h in header]
            f.write('\t'.join(row) + '\n')

    return len(features), skipped, features


def detect_curation_warnings(features):
    """Detect features that need user attention during curation."""
    warnings = []

    # Check for duplicate gene names
    name_counts = {}
    for f in features:
        name = f['gene_name']
        name_counts[name] = name_counts.get(name, 0) + 1
    duplicates = {name: count for name, count in name_counts.items() if count > 1}
    if duplicates:
        for name, count in duplicates.items():
            products = [f['product'] for f in features if f['gene_name'] == name]
            warnings.append({
                'type': 'DUPLICATE_NAMES',
                'message': f"{count} features share gene_name '{name}' - each must be unique",
                'detail': f"Products: {', '.join(products)}",
                'action': f"Rename to distinct names (e.g., {name}, {name}2, {name}3)"
            })

    # Check for join/spliced features
    joins = [f for f in features if f['location_type'] == 'join']
    if joins:
        for f in joins:
            warnings.append({
                'type': 'JOIN_FEATURE',
                'message': f"'{f['gene_name']}' has a spliced/join location: {f['location_detail']}",
                'detail': f"This may indicate RNA editing, ribosomal frameshift, or splicing.",
                'action': (f"No changes needed. Step2 will automatically generate multiple "
                           f"CDS lines in the GFF3 (one per exon) using the coordinates "
                           f"in location_detail. Review the note column for biological context.")
            })

    # Check for overlapping features on same strand
    for i, f1 in enumerate(features):
        for f2 in features[i+1:]:
            if f1['strand'] == f2['strand'] and f1['seqid'] == f2['seqid']:
                overlap_start = max(f1['start'], f2['start'])
                overlap_end = min(f1['end'], f2['end'])
                if overlap_start <= overlap_end:
                    overlap_bp = overlap_end - overlap_start + 1
                    # Skip if already flagged as duplicate name
                    if f1['gene_name'] != f2['gene_name']:
                        warnings.append({
                            'type': 'OVERLAP',
                            'message': (f"'{f1['gene_name']}' ({f1['start']}..{f1['end']}) "
                                        f"overlaps '{f2['gene_name']}' ({f2['start']}..{f2['end']}) "
                                        f"by {overlap_bp} bp"),
                            'detail': "Overlapping genes can be legitimate in viruses.",
                            'action': "Verify both are real genes. DELETE one if redundant."
                        })

    return warnings


def print_feature_summary(features, skipped, accession):
    """Print a human-readable summary of parsed features."""
    print(f"\n{'='*60}")
    print(f"Genome: {accession}")
    print(f"{'='*60}")
    cds_count = sum(1 for f in features if f['type'] == 'CDS')
    utr_count = sum(1 for f in features if f['type'] in ('5UTR', '3UTR'))
    print(f"  Features found: {len(features)} ({cds_count} CDS, {utr_count} UTR)")
    if skipped > 0:
        print(f"  Polyprotein features skipped: {skipped}")
    print()

    if features:
        print(f"  {'Gene':<20} {'Type':<6} {'Start':>8} {'End':>8} {'Strand':>6} {'Loc':>5} {'Product'}")
        print(f"  {'-'*18}  {'-'*4}  {'-'*8} {'-'*8} {'-'*6} {'-'*5} {'-'*30}")
        for f in features:
            gene = f['gene_name'][:18]
            product = f['product'][:35]
            loc = 'join' if f['location_type'] == 'join' else ''
            print(f"  {gene:<20} {f['type']:<6} {f['start']:>8} {f['end']:>8} {f['strand']:>6} {loc:>5} {product}")

    print()

    # Detect and print warnings
    warnings = detect_curation_warnings(features)
    if warnings:
        print(f"  ATTENTION - {len(warnings)} item(s) need review:")
        print(f"  {'-'*56}")
        for i, w in enumerate(warnings, 1):
            print(f"  {i}. [{w['type']}] {w['message']}")
            print(f"     {w['detail']}")
            print(f"     -> {w['action']}")
            print()


def print_curation_guide(tsv_path):
    """Print instructions for manual curation of the TSV file."""
    print(f"{'='*60}")
    print(f"MANUAL CURATION GUIDE")
    print(f"{'='*60}")
    print()
    print(f"Open the TSV file in a text editor or spreadsheet program:")
    print(f"  nano {tsv_path}")
    print(f"  -or-")
    print(f"  Open in Excel/LibreOffice (tab-delimited)")
    print()
    print(f"Key columns you can edit:")
    print(f"  gene_name  - Short gene identifier (e.g., NP, VP35, L)")
    print(f"  product    - Full protein name (e.g., nucleoprotein)")
    print(f"  start/end  - Coordinates (1-based, inclusive)")
    print(f"  strand     - + or -")
    print(f"  action     - Controls what step2 does with each row")
    print()
    print(f"HOW TO MODIFY, ADD, OR DELETE FEATURES:")
    print(f"{'='*60}")
    print()
    print(f"  KEEP a feature (default):")
    print(f"    All features start with action = KEEP.")
    print(f"    No changes needed if the annotation looks correct.")
    print()
    print(f"  MODIFY a feature:")
    print(f"    Edit gene_name, product, start, end, or strand directly.")
    print(f"    Leave the 'action' column as KEEP.")
    print(f"    Example: Change gene_name from 'hypothetical_protein' to 'NS3'")
    print()
    print(f"  DELETE a feature:")
    print(f"    Change the 'action' column from KEEP to DELETE.")
    print(f"    The row will be skipped when building the SnpEff database.")
    print(f"    Example: Remove a duplicate or spurious annotation.")
    print()
    print(f"  ADD a new feature:")
    print(f"    Copy an existing row and paste it as a new line.")
    print(f"    Edit the seqid, start, end, strand, gene_name, and product.")
    print(f"    Set action to KEEP. Set source to 'manual'.")
    print(f"    Example: Add a missing ORF you found in the literature.")
    print()
    print(f"FEATURES REQUIRING SPECIAL ATTENTION:")
    print(f"{'='*60}")
    print()
    print(f"  DUPLICATE GENE NAMES:")
    print(f"    Each gene_name MUST be unique in the TSV file.")
    print(f"    Step2 uses gene_name as the feature ID in the GFF3.")
    print(f"    If multiple features share a name, rename them.")
    print(f"    Example: Three 'GP' features -> rename to GP, sGP, ssGP")
    print()
    print(f"  JOIN/SPLICED FEATURES (location_type = 'join'):")
    print(f"    These represent genes with non-contiguous coding regions,")
    print(f"    caused by RNA editing, ribosomal frameshifts, or splicing.")
    print(f"    - The 'start' and 'end' columns span the FULL region")
    print(f"    - The 'location_detail' column shows the exact exon coordinates")
    print(f"    - Step2 automatically generates multiple CDS lines in the GFF3")
    print(f"      (one per exon), which SnpEff needs for correct spliced annotation")
    print(f"    - Review the note column for biological context")
    print(f"    Usually no changes needed, but verify the gene is real.")
    print()
    print(f"  OVERLAPPING FEATURES:")
    print(f"    Some viruses have legitimate overlapping genes.")
    print(f"    - Verify both genes are real (check literature)")
    print(f"    - DELETE redundant or spurious annotations")
    print(f"    - Keep both if they encode different proteins")
    print()
    print(f"  Columns you should NOT change:")
    print(f"    seqid  - Must match the FASTA sequence header exactly")
    print(f"    type   - Keep as 'CDS' for protein-coding features")
    print()


def main():
    parser = argparse.ArgumentParser(
        description='Download and parse viral genome annotations from NCBI (Pathway 2)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s NC_001477                          # Dengue-1
  %(prog)s NC_001474.2                        # Dengue-2
  %(prog)s NC_001477 --keep-polyprotein       # Include polyproteins
  %(prog)s NC_001477 --output-dir /data       # Custom output dir
  %(prog)s NC_001477 --email you@uni.edu      # Custom NCBI email

Output files:
  {accession}.fasta                   Genome sequence
  {accession}.gb                      GenBank record
  {accession}_no_polyprotein.tsv      Editable annotation table (default)
  {accession}_all_features.tsv        Full annotation table (--keep-polyprotein)

Next steps:
  1. Edit the TSV file to curate gene names and products
  2. Run: python3 /opt/vicast/vicast-annotate/step2_add_to_snpeff.py {accession} {accession}_no_polyprotein.tsv
        """
    )

    parser.add_argument(
        'accession',
        help='NCBI accession number (e.g., NC_001477 or NC_001477.1)'
    )
    parser.add_argument(
        '--keep-polyprotein', action='store_true',
        help='Include polyprotein features in output'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='.',
        help='Output directory (default: current directory)'
    )
    parser.add_argument(
        '--email',
        default=os.environ.get('NCBI_EMAIL', 'vicast_user@example.com'),
        help='Email for NCBI Entrez (default: $NCBI_EMAIL or vicast_user@example.com)'
    )
    parser.add_argument(
        '--local-gb',
        help='Use a local GenBank file instead of downloading from NCBI'
    )
    parser.add_argument(
        '--local-fasta',
        help='Use a local FASTA file instead of downloading from NCBI'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='Print detailed progress'
    )

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Resolve to absolute path for clear output
    output_dir = os.path.abspath(args.output_dir)

    # HOST_DIR allows displaying host paths when running inside Docker
    # Usage: docker run -e HOST_DIR=/your/host/path ...
    host_dir = os.environ.get('HOST_DIR', '')

    def host_path(container_path):
        """Replace container output_dir with host_dir for display."""
        if host_dir and container_path.startswith(output_dir):
            return container_path.replace(output_dir, host_dir, 1)
        return container_path

    accession = args.accession

    print(f"VICAST Step 1: Parse Viral Genome")
    print(f"{'='*40}")
    print(f"Accession: {accession}")
    print(f"Email: {args.email}")
    print(f"Output dir: {host_path(output_dir)}")
    print(f"Skip polyproteins: {not args.keep_polyprotein}")
    print()

    # Step 1: Get genome files
    if args.local_gb:
        gb_path = os.path.abspath(args.local_gb)
        fasta_path = os.path.abspath(args.local_fasta) if args.local_fasta else None
        if not os.path.isfile(gb_path):
            print(f"Error: GenBank file not found: {host_path(gb_path)}")
            sys.exit(1)
        print(f"Using local GenBank: {host_path(gb_path)}")
        if fasta_path:
            print(f"Using local FASTA: {host_path(fasta_path)}")
    else:
        print("Downloading from NCBI...")
        fasta_path, gb_path = download_from_ncbi(
            accession, args.email, output_dir, verbose=args.verbose
        )
        if gb_path is None:
            print("\nFailed to download genome from NCBI.")
            print(f"Check accession at: https://www.ncbi.nlm.nih.gov/nuccore/{accession}")
            sys.exit(1)
        print(f"  Downloaded: {host_path(fasta_path)}")
        print(f"  Downloaded: {host_path(gb_path)}")

    # Step 2: Parse GenBank annotations
    print("\nParsing annotations...")
    suffix = "_all_features" if args.keep_polyprotein else "_no_polyprotein"
    tsv_path = os.path.join(output_dir, f"{accession}{suffix}.tsv")

    feature_count, skipped_count, features = parse_genbank_to_tsv(
        gb_path, tsv_path,
        keep_polyprotein=args.keep_polyprotein,
        verbose=args.verbose
    )

    # Step 3: Print summary
    print_feature_summary(features, skipped_count, accession)

    print(f"Output files:")
    if fasta_path:
        print(f"  FASTA:   {host_path(fasta_path)}")
    print(f"  GenBank: {host_path(gb_path)}")
    print(f"  TSV:     {host_path(tsv_path)}")

    if feature_count == 0:
        print(f"\n{'='*60}")
        print(f"NO INDIVIDUAL PROTEIN FEATURES FOUND")
        print(f"{'='*60}")
        if skipped_count > 0:
            print(f"\n  {skipped_count} polyprotein(s) were skipped and no mat_peptide")
            print(f"  features exist to define individual protein boundaries.")
            print(f"\n  This is common for alphaviruses, picornaviruses, and other")
            print(f"  viruses where NCBI has not annotated individual mature peptides.")
        else:
            print(f"\n  No CDS or mat_peptide features found in {accession}.")

        print(f"\n  OPTIONS:")
        print(f"  {'='*56}")
        print(f"\n  1. USE PATHWAY 3 (recommended) - BLASTx with a related model:")
        print(f"     Find a related, well-annotated virus to transfer annotations.")
        print(f"     python3 /opt/vicast/vicast-annotate/step1_blastx_annotate.py --subject {accession} --model <MODEL_ACCESSION> --add-model")
        print(f"\n     Pre-built models available: install_prebuilt_database.sh --list")
        print(f"\n  2. KEEP POLYPROTEINS - rerun with --keep-polyprotein:")
        print(f"     python3 /opt/vicast/vicast-annotate/step1_parse_viral_genome.py {accession} --keep-polyprotein")
        print(f"     Then manually split polyproteins in the TSV using literature.")
        print(f"\n{'='*60}")
        sys.exit(1)

    # Step 4: Print curation guide
    print()
    print_curation_guide(host_path(tsv_path))

    print(f"When curation is complete, run step2 inside the container:")
    print(f"  python3 /opt/vicast/vicast-annotate/step2_add_to_snpeff.py {accession} {tsv_path}")


if __name__ == '__main__':
    main()
