#!/usr/bin/env python3
"""
VICAST-annotate for segmented viruses (e.g., Influenza, Rotavirus, Bunyaviruses).
Combines multiple genome segments into a single SnpEff database.
"""

import sys
import os
import argparse
import subprocess
from pathlib import Path
from Bio import SeqIO, Entrez

# Set Entrez email (users should update this)
Entrez.email = "your_email@example.com"

def download_segment(accession, output_dir="."):
    """
    Download GenBank and FASTA for a single segment.
    
    Args:
        accession: NCBI accession (e.g., CY121680)
        output_dir: Output directory
        
    Returns:
        tuple: (gb_file, fasta_file)
    """
    gb_file = os.path.join(output_dir, f"{accession}.gb")
    fasta_file = os.path.join(output_dir, f"{accession}.fasta")
    
    print(f"  Downloading {accession}...")
    
    try:
        # Download GenBank
        if not os.path.exists(gb_file):
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            with open(gb_file, 'w') as f:
                f.write(handle.read())
            handle.close()
        
        # Download FASTA
        if not os.path.exists(fasta_file):
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            with open(fasta_file, 'w') as f:
                f.write(handle.read())
            handle.close()
        
        print(f"    ✓ Downloaded {accession}")
        return gb_file, fasta_file
        
    except Exception as e:
        print(f"    ✗ Error downloading {accession}: {e}")
        return None, None

def combine_fastas(fasta_files, output_fasta, segment_names=None):
    """
    Combine multiple FASTA files into one multi-sequence FASTA.
    
    Args:
        fasta_files: List of FASTA file paths
        output_fasta: Output combined FASTA path
        segment_names: Optional list of segment names (e.g., PB2, PB1, PA...)
        
    Returns:
        int: Number of sequences combined
    """
    print(f"\nCombining {len(fasta_files)} FASTA files...")
    
    sequences = []
    for i, fasta_file in enumerate(fasta_files):
        if not os.path.exists(fasta_file):
            print(f"  ✗ Warning: {fasta_file} not found, skipping")
            continue
            
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Optionally rename sequence IDs to include segment name
            if segment_names and i < len(segment_names):
                original_id = record.id
                record.id = f"{segment_names[i]}"
                record.description = f"{segment_names[i]} {original_id} {record.description}"
            sequences.append(record)
    
    # Write combined FASTA
    SeqIO.write(sequences, output_fasta, "fasta")
    print(f"  ✓ Combined {len(sequences)} sequences → {output_fasta}")
    
    return len(sequences)

def genbank_to_gff(gb_file, output_gff, skip_polyprotein=True, chr_name=None):
    """
    Convert GenBank to GFF3, skipping polyproteins.

    Args:
        gb_file: GenBank file
        output_gff: Output GFF3 file
        skip_polyprotein: Skip polyprotein features
        chr_name: Optional chromosome name override (for segment renaming)

    Returns:
        int: Number of features written
    """
    features_written = 0

    for record in SeqIO.parse(gb_file, "genbank"):
        # Use custom chr_name if provided, otherwise use record.id
        chromosome = chr_name if chr_name else record.id

        for feature in record.features:
            # Skip polyproteins if requested (check both product and note fields)
            if skip_polyprotein and feature.type == "CDS":
                product = feature.qualifiers.get("product", [""])[0].lower()
                notes = " ".join(feature.qualifiers.get("note", [])).lower()
                if "polyprotein" in product or "polyprotein" in notes:
                    continue

            # Skip source features
            if feature.type == "source":
                continue

            # Write relevant features (exclude mRNA to prevent mistranslation)
            if feature.type in ["CDS", "gene", "3'UTR", "5'UTR", "mat_peptide"]:
                start = int(feature.location.start) + 1  # GFF is 1-based
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"

                # Convert mat_peptide to CDS for SnpEff
                output_type = "CDS" if feature.type == "mat_peptide" else feature.type

                # Build attributes - prioritize gene name for ID (more readable, matches FASTA)
                attributes = []
                if "gene" in feature.qualifiers:
                    # Use gene name as ID for consistency with FASTA files
                    gene_id = feature.qualifiers['gene'][0]
                    attributes.append(f"ID={gene_id}")
                    attributes.append(f"gene={gene_id}")
                elif "locus_tag" in feature.qualifiers:
                    attributes.append(f"ID={feature.qualifiers['locus_tag'][0]}")
                elif "protein_id" in feature.qualifiers:
                    attributes.append(f"ID={feature.qualifiers['protein_id'][0]}")
                else:
                    attributes.append(f"ID={output_type}_{features_written}")
                    if "gene" in feature.qualifiers:
                        attributes.append(f"gene={feature.qualifiers['gene'][0]}")

                if "product" in feature.qualifiers:
                    attributes.append(f"product={feature.qualifiers['product'][0]}")
                if "protein_id" in feature.qualifiers:
                    attributes.append(f"protein_id={feature.qualifiers['protein_id'][0]}")

                # Write GFF line with chromosome name
                gff_line = f"{chromosome}\tGenBank\t{output_type}\t{start}\t{end}\t.\t{strand}\t.\t{';'.join(attributes)}\n"

                with open(output_gff, 'a') as f:
                    f.write(gff_line)

                features_written += 1

    return features_written

def combine_gffs(gb_files, output_gff, skip_polyprotein=True, segment_names=None):
    """
    Convert multiple GenBank files to a combined GFF3.

    Args:
        gb_files: List of GenBank file paths
        output_gff: Output combined GFF3 path
        skip_polyprotein: Skip polyprotein features
        segment_names: Optional list of segment names (for chromosome renaming)

    Returns:
        int: Total features written
    """
    print(f"\nCombining {len(gb_files)} GenBank files to GFF3...")

    # Remove existing output file
    if os.path.exists(output_gff):
        os.remove(output_gff)

    # Write GFF3 header
    with open(output_gff, 'w') as f:
        f.write("##gff-version 3\n")

    total_features = 0
    for i, gb_file in enumerate(gb_files):
        if not os.path.exists(gb_file):
            print(f"  ✗ Warning: {gb_file} not found, skipping")
            continue

        # Use segment name if provided
        chr_name = segment_names[i] if segment_names and i < len(segment_names) else None

        features = genbank_to_gff(gb_file, output_gff, skip_polyprotein, chr_name)
        display_name = chr_name if chr_name else Path(gb_file).stem
        print(f"  ✓ {display_name}: {features} features")
        total_features += features

    print(f"  ✓ Total features: {total_features}")
    return total_features

def gff_to_tab_delimited(gff_file, output_tsv, segment_to_accession=None):
    """
    Convert GFF3 to tab-delimited format for manual editing.

    Args:
        gff_file: Path to GFF3 file
        output_tsv: Path to output TSV file
        segment_to_accession: Optional dict mapping segment names to accession numbers

    Returns:
        Path to the TSV file
    """
    import pandas as pd

    data = []

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 9:
                seqid = parts[0]
                source = parts[1]
                feature_type = parts[2]
                start = parts[3]
                end = parts[4]
                score = parts[5]
                strand = parts[6]
                phase = parts[7]
                attributes = parts[8]

                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                # Get accession for this segment if available
                accession = ''
                if segment_to_accession and seqid in segment_to_accession:
                    accession = segment_to_accession[seqid]

                data.append({
                    'action': 'KEEP',  # Default action - user can change to DELETE or MODIFY
                    'segment': seqid,
                    'accession': accession,
                    'source': source,
                    'type': feature_type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_name': attr_dict.get('gene', ''),
                    'product': attr_dict.get('product', ''),
                    'ID': attr_dict.get('ID', ''),
                    'protein_id': attr_dict.get('protein_id', ''),
                    'notes': ''  # For user to add custom notes
                })

    # Create DataFrame and save to TSV
    df = pd.DataFrame(data)
    df.to_csv(output_tsv, sep='\t', index=False)

    print(f"\n✓ Created annotation TSV: {output_tsv}")
    print(f"  Total features: {len(df)}")

    return output_tsv

def generate_cds_protein_fasta(genome_fasta, gff_file, genome_id, snpeff_data_dir):
    """
    Generate CDS and protein FASTA files for SnpEff validation.
    Uses BioPython to ensure full-length sequences with proper line wrapping.

    Args:
        genome_fasta: Path to combined genome FASTA file
        gff_file: Path to combined GFF3 file
        genome_id: Genome identifier
        snpeff_data_dir: SnpEff data directory

    Returns:
        tuple: (cds_file, protein_file) paths
    """
    from Bio.Seq import Seq

    print(f"\nGenerating CDS and protein FASTA files...")

    # Load all genome sequences into a dictionary using BioPython
    genome_seqs = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        genome_seqs[record.id] = record.seq

    print(f"  Loaded {len(genome_seqs)} genome segments")

    # Parse GFF3 to extract CDS features
    cds_features = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                seqid = parts[0]
                start = int(parts[3]) - 1  # Convert to 0-based
                end = int(parts[4])
                strand = parts[6]
                attributes = parts[8]

                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value

                cds_features.append({
                    'seqid': seqid,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'id': attr_dict.get('ID', ''),
                    'gene': attr_dict.get('gene', ''),
                    'product': attr_dict.get('product', 'hypothetical protein')
                })

    print(f"  Found {len(cds_features)} CDS features")

    # Generate CDS and protein sequences
    cds_records = []
    protein_records = []

    for feature in cds_features:
        seqid = feature['seqid']

        if seqid not in genome_seqs:
            print(f"    Warning: Sequence {seqid} not found in FASTA")
            continue

        # Extract FULL CDS sequence using BioPython slicing
        cds_seq = genome_seqs[seqid][feature['start']:feature['end']]

        # Reverse complement if on minus strand
        if feature['strand'] == '-':
            cds_seq = cds_seq.reverse_complement()

        # Generate sequence ID
        seq_id = feature['gene'] if feature['gene'] else feature['id']

        # Create CDS record
        cds_record = SeqIO.SeqRecord(
            cds_seq,
            id=seq_id,
            description=feature['product']
        )
        cds_records.append(cds_record)

        # Translate to protein (FULL sequence)
        try:
            protein_seq = cds_seq.translate(to_stop=True)

            # If doesn't start with M, try alternate start
            if not str(protein_seq).startswith('M') and len(cds_seq) >= 3:
                protein_seq = Seq('M') + cds_seq[3:].translate(to_stop=True)

            protein_record = SeqIO.SeqRecord(
                protein_seq,
                id=seq_id,
                description=feature['product']
            )
            protein_records.append(protein_record)

        except Exception as e:
            print(f"    Warning: Could not translate {seq_id}: {e}")
            protein_record = SeqIO.SeqRecord(
                Seq(''),
                id=seq_id,
                description=f"{feature['product']} [translation failed]"
            )
            protein_records.append(protein_record)

    # Write files using BioPython (handles full-length sequences with 60 char/line wrapping)
    genome_dir = os.path.join(snpeff_data_dir, genome_id)
    os.makedirs(genome_dir, exist_ok=True)

    cds_file = os.path.join(genome_dir, 'cds.fa')
    protein_file = os.path.join(genome_dir, 'protein.fa')

    SeqIO.write(cds_records, cds_file, "fasta")
    SeqIO.write(protein_records, protein_file, "fasta")

    print(f"  Generated CDS FASTA: {cds_file} ({len(cds_records)} sequences)")
    print(f"  Generated protein FASTA: {protein_file} ({len(protein_records)} sequences)")

    return cds_file, protein_file

def add_to_snpeff(genome_id, fasta_file, gff_file, snpeff_data_dir=None):
    """
    Add combined segmented genome to SnpEff.
    
    Args:
        genome_id: Genome identifier
        fasta_file: Combined FASTA file
        gff_file: Combined GFF3 file
        snpeff_data_dir: SnpEff data directory
        
    Returns:
        bool: Success status
    """
    # Determine SnpEff data directory
    if not snpeff_data_dir:
        snpeff_home = os.environ.get('SNPEFF_HOME')
        if not snpeff_home:
            print("ERROR: SNPEFF_HOME not set. Please source snpEff environment.")
            return False
        snpeff_data_dir = os.path.join(snpeff_home, 'data')
    
    # Create genome directory
    genome_dir = os.path.join(snpeff_data_dir, genome_id)
    os.makedirs(genome_dir, exist_ok=True)
    
    # Copy files to SnpEff structure
    import shutil
    sequences_file = os.path.join(genome_dir, 'sequences.fa')
    genes_file = os.path.join(genome_dir, 'genes.gff')
    
    print(f"\nAdding to SnpEff...")
    shutil.copy(fasta_file, sequences_file)
    shutil.copy(gff_file, genes_file)
    
    print(f"  Genome directory: {genome_dir}")
    print(f"  Sequences: {sequences_file}")
    print(f"  Annotations: {genes_file}")

    # Generate CDS and protein FASTA files for validation
    try:
        generate_cds_protein_fasta(fasta_file, gff_file, genome_id, snpeff_data_dir)
    except Exception as e:
        print(f"  Warning: Could not generate CDS/protein files: {e}")
        print("  Continuing without validation files...")

    # Update snpEff.config
    config_file = os.path.join(os.path.dirname(snpeff_data_dir), 'snpEff.config')
    
    genome_in_config = False
    if os.path.exists(config_file):
        with open(config_file, 'r') as f:
            if f"{genome_id}.genome" in f.read():
                genome_in_config = True
    
    if not genome_in_config:
        print(f"\nAdding {genome_id} to snpEff.config...")
        with open(config_file, 'a') as f:
            f.write(f"\n# {genome_id} - Segmented virus genome\n")
            f.write(f"{genome_id}.genome : {genome_id}\n")
        print(f"  ✓ Added to config")
    
    # Build database
    print(f"\nBuilding SnpEff database...")
    try:
        snpeff_jar = os.path.join(snpeff_home, 'snpEff.jar')
        java_home = os.environ.get('JAVA_HOME')
        
        if java_home and os.path.exists(snpeff_jar):
            java_cmd = os.path.join(java_home, 'bin', 'java')
            cmd = [java_cmd, '-jar', snpeff_jar, 'build', '-gff3', '-v',
                   '-noCheckCds', '-noCheckProtein', genome_id]
        else:
            cmd = ['snpeff', 'build', '-gff3', '-v',
                   '-noCheckCds', '-noCheckProtein', genome_id]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print(f"\n✓ Successfully built SnpEff database for {genome_id}")
            return True
        else:
            print(f"\n✗ Error building database:")
            print(result.stderr)
            return False
            
    except Exception as e:
        print(f"\n✗ Error running SnpEff: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(
        description='VICAST-annotate for segmented viruses',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script handles segmented viruses (e.g., Influenza, Rotavirus) by:
1. Downloading all segments from NCBI
2. Combining sequences into single multi-sequence FASTA
3. Combining annotations into single multi-feature GFF3
4. Building a unified SnpEff database

Examples:

  # Influenza A (H1N1) with 8 segments
  python3 vicast_annotate_segmented.py influenza_h1n1_2009 \\
    --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \\
    --names PB2,PB1,PA,HA,NP,NA,M,NS

  # Rotavirus with 11 segments
  python3 vicast_annotate_segmented.py rotavirus_strain_X \\
    --segments NC_XXX1,NC_XXX2,NC_XXX3,NC_XXX4,NC_XXX5,NC_XXX6,NC_XXX7,NC_XXX8,NC_XXX9,NC_XXX10,NC_XXX11 \\
    --names VP1,VP2,VP3,VP4,NSP1,VP6,NSP3,NSP2,VP7,NSP4,NSP5

  # Use local files instead of downloading
  python3 vicast_annotate_segmented.py my_virus \\
    --fasta-files seg1.fasta,seg2.fasta,seg3.fasta \\
    --gb-files seg1.gb,seg2.gb,seg3.gb \\
    --names Seg1,Seg2,Seg3

After successful completion:
  snpeff influenza_h1n1_2009 variants.vcf > annotated.vcf
        """
    )
    
    parser.add_argument('genome_id',
                       help='Genome identifier for SnpEff database (e.g., influenza_h1n1_2009)')
    parser.add_argument('--segments', 
                       help='Comma-separated list of NCBI accessions (e.g., CY121680,CY121681,...)')
    parser.add_argument('--names',
                       help='Comma-separated list of segment names (e.g., PB2,PB1,PA,...)')
    parser.add_argument('--fasta-files',
                       help='Use local FASTA files instead of downloading (comma-separated)')
    parser.add_argument('--gb-files',
                       help='Use local GenBank files instead of downloading (comma-separated)')
    parser.add_argument('--output-dir', default='.',
                       help='Output directory for intermediate files (default: current)')
    parser.add_argument('--skip-polyprotein', action='store_true', default=True,
                       help='Skip polyprotein features (default: True)')
    parser.add_argument('--keep-polyprotein', action='store_false', dest='skip_polyprotein',
                       help='Keep polyprotein features')
    parser.add_argument('--no-review', action='store_true',
                       help='Skip TSV review step and proceed directly to snpEff build')
    
    args = parser.parse_args()
    
    print("="*60)
    print("VICAST-annotate: Segmented Virus Pipeline")
    print("="*60)
    print(f"\nGenome ID: {args.genome_id}")
    
    # Parse segment names
    segment_names = None
    if args.names:
        segment_names = [name.strip() for name in args.names.split(',')]
        print(f"Segment names: {', '.join(segment_names)}")
    
    # Determine input method
    if args.fasta_files and args.gb_files:
        # Use local files
        print("\nUsing local files...")
        fasta_files = [f.strip() for f in args.fasta_files.split(',')]
        gb_files = [f.strip() for f in args.gb_files.split(',')]
        
        # Verify files exist
        for f in fasta_files + gb_files:
            if not os.path.exists(f):
                print(f"Error: File not found: {f}")
                sys.exit(1)
        
    elif args.segments:
        # Download from NCBI
        print("\nDownloading segments from NCBI...")
        accessions = [acc.strip() for acc in args.segments.split(',')]
        print(f"Segments: {len(accessions)}")
        
        os.makedirs(args.output_dir, exist_ok=True)
        
        fasta_files = []
        gb_files = []
        
        for accession in accessions:
            gb_file, fasta_file = download_segment(accession, args.output_dir)
            if gb_file and fasta_file:
                gb_files.append(gb_file)
                fasta_files.append(fasta_file)
        
        if not fasta_files:
            print("\nError: No segments downloaded successfully")
            sys.exit(1)
            
    else:
        print("\nError: Must provide either --segments OR (--fasta-files AND --gb-files)")
        sys.exit(1)
    
    # Combine FASTA files
    combined_fasta = os.path.join(args.output_dir, f"{args.genome_id}_combined.fasta")
    num_seqs = combine_fastas(fasta_files, combined_fasta, segment_names)
    
    if num_seqs == 0:
        print("\nError: No sequences to combine")
        sys.exit(1)
    
    # Combine GFF files
    combined_gff = os.path.join(args.output_dir, f"{args.genome_id}_combined.gff3")
    num_features = combine_gffs(gb_files, combined_gff, args.skip_polyprotein, segment_names)

    if num_features == 0:
        print("\nError: No features to combine")
        sys.exit(1)

    # Create segment-to-accession mapping for TSV
    segment_to_accession = {}
    if segment_names and 'accessions' in locals():
        # Map segment names to accessions if both are available
        if len(segment_names) == len(accessions):
            segment_to_accession = dict(zip(segment_names, accessions))

    # Generate TSV for manual review
    combined_tsv = os.path.join(args.output_dir, f"{args.genome_id}_annotations.tsv")
    gff_to_tab_delimited(combined_gff, combined_tsv, segment_to_accession)

    # Check if review is needed
    if not args.no_review:
        print("\n" + "="*60)
        print("REVIEW ANNOTATIONS")
        print("="*60)
        print(f"\n1. REVIEW the annotation file:")
        print(f"   {combined_tsv}")
        print("\n   You can open this in Excel or a text editor.")
        print("   - Check gene names and products")
        print("   - Mark unwanted features with action='DELETE'")
        print("   - Add missing annotations")
        print("   - Correct any errors")
        print(f"\n2. SAVE your edits (keep as TSV format)")
        print(f"\n3. PRESS ENTER to continue building snpEff database...")
        print("   (or Ctrl+C to abort)")

        try:
            input()
        except KeyboardInterrupt:
            print("\n\nAborted by user.")
            sys.exit(0)

        print("\nContinuing with snpEff build...")

    # Add to SnpEff
    success = add_to_snpeff(args.genome_id, combined_fasta, combined_gff)
    
    if success:
        print("\n" + "="*60)
        print("SUCCESS!")
        print("="*60)
        print(f"\nSegmented genome {args.genome_id} added to SnpEff")
        print(f"\nUsage:")
        print(f"  snpeff {args.genome_id} variants.vcf > annotated.vcf")
        print("\nNote: VCF CHROM field must match segment names in FASTA")
        if segment_names:
            print(f"  Expected CHROM values: {', '.join(segment_names)}")
        print("\n" + "="*60)
    else:
        print("\nFailed to add genome to SnpEff")
        sys.exit(1)

if __name__ == '__main__':
    main()
