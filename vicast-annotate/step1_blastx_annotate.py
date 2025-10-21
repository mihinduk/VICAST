#!/usr/bin/env python3
"""
STEP 1 (Pathway 3): BLASTx-based annotation for poorly annotated genomes.
This script uses homology-based annotation when NCBI annotations are poor/missing.
"""

import sys
import os
import argparse
import pandas as pd
import subprocess
from pathlib import Path
from Bio import SeqIO
try:
    from Bio.Blast import NCBIXML
except ImportError:
    print("Warning: BioPython BLAST XML parser not available")
import re

def run_blastx(fasta_file, blast_db='nr', evalue=1e-5, max_target_seqs=5):
    """
    Run BLASTx on genome against protein database.
    
    Args:
        fasta_file: Path to genome FASTA
        blast_db: BLAST database to use
        evalue: E-value threshold
        max_target_seqs: Maximum number of target sequences
        
    Returns:
        Path to BLAST tabular output file
    """
    output_file = f"{Path(fasta_file).stem}_blastx.txt"
    
    print(f"\nRunning BLASTx annotation...")
    print(f"  Database: {blast_db}")
    print(f"  E-value: {evalue}")
    print("  This may take several minutes...")
    
    try:
        cmd = [
            'blastx',
            '-query', fasta_file,
            '-db', blast_db,
            '-evalue', str(evalue),
            '-max_target_seqs', str(max_target_seqs),
            '-outfmt', '6 qseqid qstart qend sseqid pident evalue stitle',
            '-out', output_file
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        
        if result.returncode == 0:
            print(f"  ✓ BLASTx completed: {output_file}")
            return output_file
        else:
            print(f"  ✗ BLASTx failed:")
            print(result.stderr)
            return None
            
    except subprocess.TimeoutExpired:
        print("  ✗ BLASTx timed out (>60 minutes)")
        return None
    except FileNotFoundError:
        print("  ✗ BLASTx not found. Please install BLAST+:")
        print("     conda install -c bioconda blast")
        return None
    except Exception as e:
        print(f"  ✗ Error running BLASTx: {e}")
        return None

def parse_blastx_tabular(blast_file):
    """
    Parse BLASTx tabular output and extract annotations.
    
    Args:
        blast_file: Path to BLAST tabular file
        
    Returns:
        List of annotation dictionaries
    """
    annotations = []
    
    print(f"\nParsing BLASTx results...")
    
    try:
        with open(blast_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue
                
                qseqid, qstart, qend, sseqid, pident, evalue, stitle = parts[:7]
                
                # Extract gene/product info from hit title
                product = extract_product_from_title(stitle)
                gene = extract_gene_from_title(stitle)
                
                # Determine strand based on coordinates
                qstart, qend = int(qstart), int(qend)
                strand = '+' if qstart < qend else '-'
                if strand == '-':
                    qstart, qend = qend, qstart
                
                annotations.append({
                    'query': qseqid,
                    'start': qstart,
                    'end': qend,
                    'strand': strand,
                    'product': product,
                    'gene': gene,
                    'evalue': float(evalue),
                    'identity': float(pident),
                    'hit_accession': sseqid,
                    'hit_title': stitle[:100]
                })
        
        print(f"  Found {len(annotations)} BLAST hits")
        return annotations
        
    except Exception as e:
        print(f"  ✗ Error parsing BLAST results: {e}")
        return []

def extract_product_from_title(title):
    """Extract product name from BLAST hit title."""
    # Remove accession at start
    match = re.search(r'\s+([A-Za-z].*?)(?:\[|$)', title)
    if match:
        return match.group(1).strip()
    return title[:50]

def extract_gene_from_title(title):
    """Extract gene name from BLAST hit title."""
    # Look for common gene name patterns
    match = re.search(r'\b([A-Z]\d+[A-Z]?|NS\d+[A-Z]?)\b', title)
    if match:
        return match.group(1)
    return ''

def merge_overlapping_hits(annotations, overlap_threshold=0.8):
    """
    Merge overlapping BLAST hits, keeping the best one.
    
    Args:
        annotations: List of annotation dicts
        overlap_threshold: Minimum overlap fraction to merge
        
    Returns:
        Filtered list of annotations
    """
    if not annotations:
        return []
    
    # Sort by start position and identity
    sorted_annots = sorted(annotations, key=lambda x: (x['start'], -x['identity']))
    
    merged = []
    current = sorted_annots[0].copy()
    
    for annot in sorted_annots[1:]:
        # Check for overlap
        overlap_start = max(current['start'], annot['start'])
        overlap_end = min(current['end'], annot['end'])
        overlap_len = max(0, overlap_end - overlap_start)
        
        current_len = current['end'] - current['start']
        annot_len = annot['end'] - annot['start']
        
        overlap_frac = overlap_len / min(current_len, annot_len) if min(current_len, annot_len) > 0 else 0
        
        if overlap_frac >= overlap_threshold:
            # Keep the better hit (higher identity)
            if annot['identity'] > current['identity']:
                current = annot.copy()
        else:
            # No significant overlap, keep current and move to next
            merged.append(current)
            current = annot.copy()
    
    # Add the last one
    merged.append(current)
    
    print(f"  Merged to {len(merged)} non-overlapping features")
    return merged

def create_annotation_tsv(annotations, output_tsv, seqid):
    """
    Create TSV file from BLASTx annotations.
    
    Args:
        annotations: List of annotation dicts
        output_tsv: Path to output TSV
        seqid: Sequence ID
    """
    data = []
    
    for i, annot in enumerate(annotations, 1):
        data.append({
            'action': 'KEEP',
            'seqid': seqid,
            'source': 'BLASTx',
            'type': 'CDS',
            'start': annot['start'],
            'end': annot['end'],
            'strand': annot['strand'],
            'gene_name': annot.get('gene', f'ORF{i}'),
            'product': annot['product'],
            'ID': f"cds_{i}",
            'gene': annot.get('gene', ''),
            'protein_id': '',
            'notes': f"BLASTx hit: {annot['hit_accession']} (E={annot['evalue']:.2e}, ID={annot['identity']:.1f}%)"
        })
    
    df = pd.DataFrame(data)
    df.to_csv(output_tsv, sep='\t', index=False)
    
    print(f"\nCreated annotation TSV: {output_tsv}")
    print(f"  Total features: {len(df)}")
    
    return output_tsv

def main():
    parser = argparse.ArgumentParser(
        description='STEP 1 (Pathway 3): BLASTx-based annotation for poorly annotated genomes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Use this pathway when:
  - NCBI has no or very poor annotations
  - You need homology-based annotation

Workflow:
  1. Run BLASTx against protein database
  2. Parse and extract best hits
  3. Merge overlapping hits
  4. Create editable TSV for manual curation

Requirements:
  - BLAST+ installed (conda install -c bioconda blast)
  - Access to BLAST database (nr or viral protein database)

Examples:
  # Use local BLAST database
  python3 step1_blastx_annotate.py NC_XXXXXX.fasta --blast-db /path/to/nr
  
  # Quick annotation with viral database
  python3 step1_blastx_annotate.py NC_XXXXXX.fasta --blast-db refseq_protein
  
  # Adjust sensitivity
  python3 step1_blastx_annotate.py NC_XXXXXX.fasta --evalue 1e-10

After running:
  Review and edit the TSV file, then:
    python3 step2_add_to_snpeff.py genome_id output.tsv
        """
    )
    
    parser.add_argument('fasta_file',
                       help='Genome FASTA file')
    parser.add_argument('--blast-db', default='nr',
                       help='BLAST database to use (default: nr)')
    parser.add_argument('--evalue', type=float, default=1e-5,
                       help='E-value threshold (default: 1e-5)')
    parser.add_argument('--output', '-o',
                       help='Output base name (default: from fasta filename)')
    
    args = parser.parse_args()
    
    # Check input file
    if not os.path.exists(args.fasta_file):
        print(f"Error: FASTA file not found: {args.fasta_file}")
        sys.exit(1)
    
    # Determine output names
    genome_id = Path(args.fasta_file).stem
    output_base = args.output if args.output else genome_id
    output_tsv = f"{output_base}_blastx.tsv"
    
    print("="*60)
    print("STEP 1 (Pathway 3): BLASTx-based Annotation")
    print("="*60)
    print(f"\nGenome: {args.fasta_file}")
    print(f"BLAST DB: {args.blast_db}")
    
    # Run BLASTx
    blast_file = run_blastx(args.fasta_file, args.blast_db, args.evalue)
    
    if not blast_file or not os.path.exists(blast_file):
        print("\nError: BLASTx failed. Cannot proceed.")
        sys.exit(1)
    
    # Parse results
    annotations = parse_blastx_tabular(blast_file)
    
    if not annotations:
        print("\nWarning: No BLAST hits found!")
        print("Consider:")
        print("  - Using a different BLAST database")
        print("  - Relaxing the E-value threshold (--evalue 1e-3)")
        print("  - Checking if genome is truly novel")
        sys.exit(1)
    
    # Merge overlapping hits
    print("\nMerging overlapping annotations...")
    merged_annotations = merge_overlapping_hits(annotations)
    
    # Get sequence ID from FASTA
    seqid = genome_id
    for record in SeqIO.parse(args.fasta_file, "fasta"):
        seqid = record.id
        break
    
    # Create TSV
    create_annotation_tsv(merged_annotations, output_tsv, seqid)
    
    print("\n" + "="*60)
    print("NEXT STEPS:")
    print("="*60)
    print(f"\n1. REVIEW AND EDIT the TSV file:")
    print(f"   {output_tsv}")
    print("\n   BLASTx annotations may need refinement:")
    print("   - Verify gene boundaries")
    print("   - Check for missed ORFs")
    print("   - Refine gene names and products")
    print("   - Remove false positives")
    
    print(f"\n2. SAVE your edited file (keep as TSV)")
    
    print(f"\n3. RUN STEP 2 to add to snpEff:")
    print(f"   python3 step2_add_to_snpeff.py {genome_id} {output_tsv}")
    print("\n" + "="*60)

if __name__ == '__main__':
    main()
