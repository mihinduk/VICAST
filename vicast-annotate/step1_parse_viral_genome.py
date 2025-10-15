#!/usr/bin/env python3
"""
STEP 1: Parse viral genome GenBank file and create editable TSV.
This script parses GenBank files, skips polyproteins, and creates a TSV for manual review.
"""

import sys
import os
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import re

def parse_gene_name_from_product(product):
    """
    Parse gene name from product description.
    Extracts the most relevant part for gene naming.
    
    Args:
        product: Product description from GenBank
        
    Returns:
        Parsed gene name (user can edit as needed)
    """
    if not product:
        return ''
    
    # Clean up the product string
    product_clean = product.strip()
    
    # Look for common patterns and extract gene-like names
    # Pattern 1: "protein XYZ" -> XYZ
    match = re.search(r'protein\s+([A-Za-z0-9_-]+)', product_clean, re.IGNORECASE)
    if match:
        return match.group(1)
    
    # Pattern 2: "XYZ protein" -> XYZ  
    match = re.search(r'([A-Za-z0-9_-]+)\s+protein', product_clean, re.IGNORECASE)
    if match:
        return match.group(1)
    
    # Pattern 3: Look for NS patterns (NS1, NS2A, etc.)
    match = re.search(r'\b(NS\d+[A-Za-z]?)\b', product_clean, re.IGNORECASE)
    if match:
        return match.group(1).upper()
    
    # Pattern 4: Single letter/short names at end
    match = re.search(r'\b([A-Za-z]\w?)\s*$', product_clean)
    if match:
        return match.group(1)
    
    # If no clear pattern, return first word or abbreviation
    words = product_clean.split()
    if words:
        first_word = words[0]
        # If first word looks like an abbreviation, use it
        if len(first_word) <= 4 or first_word.isupper():
            return first_word
    
    # Default: return empty for user to fill
    return ''

def extract_frameshift_position(note, feature_start, feature_end):
    """
    Extract frameshift position from note text.
    
    Args:
        note: Note text that mentions frameshift
        feature_start: Start position of the feature
        feature_end: End position of the feature
        
    Returns:
        Position of frameshift or None if not found
    """
    # Look for position patterns in the note
    import re
    
    # Pattern: "position 1234" or "at 1234" or "nucleotide 1234"
    pos_match = re.search(r'(?:position|at|nucleotide)\s+(\d+)', note, re.IGNORECASE)
    if pos_match:
        return int(pos_match.group(1))
    
    # Pattern: just a number in the note
    num_match = re.search(r'\b(\d+)\b', note)
    if num_match:
        pos = int(num_match.group(1))
        # Check if it's within the feature range
        if feature_start <= pos <= feature_end:
            return pos
    
    # Default: assume middle of feature
    return (feature_start + feature_end) // 2

def parse_genbank_skip_polyprotein(gb_file, output_gff):
    """
    Parse GenBank file and create GFF3, skipping polyprotein features.
    
    Args:
        gb_file: Path to GenBank file
        output_gff: Path to output GFF3 file
    
    Returns:
        Number of features written (excluding polyproteins)
    """
    features_written = 0
    polyproteins_skipped = 0
    
    with open(output_gff, 'w') as gff_out:
        # Write GFF3 header
        gff_out.write("##gff-version 3\n")
        
        for record in SeqIO.parse(gb_file, "genbank"):
            # Write sequence region
            gff_out.write(f"##sequence-region {record.id} 1 {len(record.seq)}\n")
            
            for feature in record.features:
                # Skip polyprotein features
                if feature.type == "CDS":
                    product = feature.qualifiers.get("product", [""])[0].lower()
                    if "polyprotein" in product:
                        polyproteins_skipped += 1
                        print(f"  Skipping polyprotein: {feature.qualifiers.get('product', [''])[0]}")
                        continue
                
                # Skip source features (they're not needed for snpEff)
                if feature.type == "source":
                    continue
                
                # Write feature to GFF3
                if feature.type in ["CDS", "gene", "mRNA", "3'UTR", "5'UTR", "mat_peptide"]:
                    start = int(feature.location.start) + 1  # GFF is 1-based
                    end = int(feature.location.end)
                    strand = "+" if feature.location.strand == 1 else "-"
                    
                    # Convert mat_peptide to CDS for better snpEff annotation
                    output_type = "CDS" if feature.type == "mat_peptide" else feature.type
                    
                    # Build attributes
                    attributes = []
                    if "locus_tag" in feature.qualifiers:
                        attributes.append(f"ID={feature.qualifiers['locus_tag'][0]}")
                    if "gene" in feature.qualifiers:
                        attributes.append(f"gene={feature.qualifiers['gene'][0]}")
                    if "product" in feature.qualifiers:
                        attributes.append(f"product={feature.qualifiers['product'][0]}")
                    if "protein_id" in feature.qualifiers:
                        attributes.append(f"protein_id={feature.qualifiers['protein_id'][0]}")
                    
                    if not attributes:
                        attributes.append(f"ID={output_type}_{features_written}")
                    
                    # Write GFF line with converted type
                    gff_line = f"{record.id}\tGenBank\t{output_type}\t{start}\t{end}\t.\t{strand}\t.\t{';'.join(attributes)}\n"
                    gff_out.write(gff_line)
                    features_written += 1
                    
                    # Check for programmed frameshifts in this feature
                    if "note" in feature.qualifiers:
                        for note in feature.qualifiers["note"]:
                            if "frameshift" in note.lower():
                                # Extract frameshift position if mentioned
                                frameshift_pos = extract_frameshift_position(note, start, end)
                                if frameshift_pos:
                                    # Add frameshift as separate feature
                                    fs_attributes = [f"ID=frameshift_{features_written}", 
                                                   f"product=programmed frameshift", 
                                                   f"note={note}"]
                                    fs_line = f"{record.id}\tGenBank\tframeshift\t{frameshift_pos}\t{frameshift_pos}\t.\t{strand}\t.\t{';'.join(fs_attributes)}\n"
                                    gff_out.write(fs_line)
                                    features_written += 1
    
    print(f"\nProcessed GenBank file: {gb_file}")
    print(f"  Features written: {features_written}")
    print(f"  Polyproteins skipped: {polyproteins_skipped}")
    
    return features_written

def gff_to_tab_delimited(gff_file, output_tsv):
    """
    Convert GFF3 to tab-delimited format for manual editing.
    
    Args:
        gff_file: Path to GFF3 file
        output_tsv: Path to output TSV file
    
    Returns:
        Path to the TSV file
    """
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
                
                # Parse gene name from product description
                product_name = attr_dict.get('product', '')
                gene_name = parse_gene_name_from_product(product_name)
                
                # Convert mat_peptide to CDS for consistency
                display_type = "CDS" if feature_type == "mat_peptide" else feature_type
                
                data.append({
                    'action': 'KEEP',  # Default action - user can change to DELETE or MODIFY
                    'seqid': seqid,
                    'source': source,
                    'type': display_type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_name': gene_name,
                    'product': product_name,
                    'ID': attr_dict.get('ID', ''),
                    'gene': attr_dict.get('gene', ''),
                    'protein_id': attr_dict.get('protein_id', ''),
                    'notes': ''  # For user to add custom notes
                })
    
    # Create DataFrame and save to TSV
    df = pd.DataFrame(data)
    df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"\nCreated tab-delimited file: {output_tsv}")
    print(f"  Total features: {len(df)}")
    
    return output_tsv

def main():
    parser = argparse.ArgumentParser(
        description='STEP 1: Parse viral genome GenBank file and create editable TSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process downloaded genome
  python3 step1_parse_viral_genome.py NC_009942.1
  
  # Process with custom output names
  python3 step1_parse_viral_genome.py NC_009942.1 --output my_genome
  
After running this script:
  1. Review the generated TSV file
  2. Make any necessary edits (add/remove/modify features)
  3. Run step2_add_to_snpeff.py with the edited TSV
        """
    )
    
    parser.add_argument('genome_id', 
                       help='Genome ID (e.g., NC_009942.1) or path to GenBank file')
    parser.add_argument('--output', '-o',
                       help='Output base name (default: genome_id)')
    
    args = parser.parse_args()
    
    # Determine input files
    if os.path.exists(args.genome_id) and args.genome_id.endswith('.gb'):
        # Direct path to GB file
        gb_file = args.genome_id
        genome_id = Path(gb_file).stem
    else:
        # Genome ID provided, look for downloaded files
        genome_id = args.genome_id
        gb_file = f"{genome_id}.gb"
        
        if not os.path.exists(gb_file):
            print(f"Error: GenBank file not found: {gb_file}")
            print(f"\nPlease run: download_ncbi_genome {genome_id}")
            sys.exit(1)
    
    # Determine output names
    output_base = args.output if args.output else genome_id
    gff_file = f"{output_base}_no_polyprotein.gff3"
    tsv_file = f"{output_base}_no_polyprotein.tsv"
    
    print("="*60)
    print("STEP 1: Parse Viral Genome and Create Editable TSV")
    print("="*60)
    
    # Parse GenBank and skip polyproteins
    print(f"\nParsing GenBank file and skipping polyproteins...")
    parse_genbank_skip_polyprotein(gb_file, gff_file)
    
    # Convert to TSV for editing
    print(f"\nConverting to tab-delimited format for editing...")
    gff_to_tab_delimited(gff_file, tsv_file)
    
    print("\n" + "="*60)
    print("NEXT STEPS:")
    print("="*60)
    print(f"\n1. REVIEW AND EDIT the TSV file:")
    print(f"   {tsv_file}")
    print("\n   You can open this in Excel, a text editor, or any spreadsheet program.")
    print("   - Check gene names and products")
    print("   - Add missing annotations")
    print("   - Remove unwanted features")
    print("   - Correct any errors")
    
    print(f"\n2. SAVE your edited file (keep it as TSV format)")
    
    print(f"\n3. RUN STEP 2 to add to snpEff:")
    print(f"   python3 step2_add_to_snpeff.py {genome_id} {tsv_file}")
    print("\n" + "="*60)

if __name__ == '__main__':
    main()