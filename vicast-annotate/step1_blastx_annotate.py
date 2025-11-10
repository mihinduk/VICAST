#!/usr/bin/env python3
"""
STEP 1 (Pathway 3): Model-based homology annotation for poorly annotated genomes.

This script annotates poorly-annotated viral genomes using well-annotated model viruses.
Three modes of operation:
  1. Model already in SnpEff database
  2. Model needs to be added to SnpEff first (auto-run Pathway 2)
  3. Custom model proteins provided by user

Systematic workflow ensures high-quality annotation transfer from related viruses.
"""

import sys
import os
import argparse
import pandas as pd
import subprocess
from pathlib import Path
from Bio import SeqIO, Entrez
import re
import tempfile
import shutil

# Set email for NCBI Entrez
Entrez.email = "vicast@example.com"

#=============================================================================
# SNPEFF DATABASE FUNCTIONS
#=============================================================================

def get_snpeff_paths():
    """Get SnpEff paths from environment variables."""
    snpeff_jar = os.environ.get('SNPEFF_JAR')
    snpeff_data = os.environ.get('SNPEFF_DATA')

    if not snpeff_jar or not snpeff_data:
        print("\n" + "="*60)
        print("ERROR: SnpEff environment variables not set!")
        print("="*60)
        print("\nPlease set these variables:")
        print("  export SNPEFF_JAR=/path/to/snpEff.jar")
        print("  export SNPEFF_DATA=/path/to/snpEff/data")
        print("\nOn HTCF:")
        print("  export SNPEFF_JAR=/ref/sahlab/software/snpEff/snpEff.jar")
        print("  export SNPEFF_DATA=/ref/sahlab/software/snpEff/data")
        print("\nAdd to ~/.bashrc for permanent setup.")
        print("="*60)
        sys.exit(1)

    return snpeff_jar, snpeff_data

def check_model_in_snpeff(model_id):
    """
    Check if model virus is already in SnpEff database.

    Args:
        model_id: Model genome ID (e.g., NC_001477)

    Returns:
        tuple: (is_available, database_path)
    """
    snpeff_jar, snpeff_data = get_snpeff_paths()

    # Check if model directory exists with required files
    model_dir = os.path.join(snpeff_data, model_id)

    if os.path.exists(model_dir):
        protein_file = os.path.join(model_dir, 'protein.fa')
        cds_file = os.path.join(model_dir, 'cds.fa')
        genes_file = os.path.join(model_dir, 'genes.gff')

        # Must have at least protein or CDS file
        if os.path.exists(protein_file) or os.path.exists(cds_file):
            return True, model_dir

    return False, None

def extract_proteins_from_snpeff(model_id):
    """
    Extract protein sequences from SnpEff database.

    Args:
        model_id: Model genome ID

    Returns:
        Path to protein FASTA file, or None if not found
    """
    is_available, model_dir = check_model_in_snpeff(model_id)

    if not is_available:
        print(f"  ✗ Model {model_id} not found in SnpEff database")
        return None

    protein_file = os.path.join(model_dir, 'protein.fa')

    if os.path.exists(protein_file):
        print(f"  ✓ Found protein sequences: {protein_file}")
        return protein_file

    # Try CDS file as fallback
    cds_file = os.path.join(model_dir, 'cds.fa')
    if os.path.exists(cds_file):
        print(f"  ✓ Found CDS sequences (will use for BLASTx): {cds_file}")
        return cds_file

    print(f"  ✗ No protein/CDS sequences found in {model_dir}")
    return None

#=============================================================================
# GENOME DOWNLOAD FUNCTIONS
#=============================================================================

def download_genome(accession):
    """
    Download genome from NCBI.

    Args:
        accession: NCBI accession (e.g., NC_XXXXXX)

    Returns:
        tuple: (genbank_file, fasta_file) or (None, None) if failed
    """
    gb_file = f"{accession}.gb"
    fasta_file = f"{accession}.fasta"

    # Check if already downloaded
    if os.path.exists(gb_file) and os.path.exists(fasta_file):
        print(f"  ✓ Using existing files: {gb_file}, {fasta_file}")
        return gb_file, fasta_file

    print(f"  Downloading {accession} from NCBI...")

    try:
        # Download GenBank
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="gbwithparts", retmode="text")
        with open(gb_file, 'w') as f:
            f.write(handle.read())
        handle.close()

        # Download FASTA
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="fasta", retmode="text")
        with open(fasta_file, 'w') as f:
            f.write(handle.read())
        handle.close()

        print(f"  ✓ Downloaded: {gb_file}, {fasta_file}")
        return gb_file, fasta_file

    except Exception as e:
        print(f"  ✗ Download failed: {e}")
        return None, None

#=============================================================================
# PATHWAY 2 INTEGRATION
#=============================================================================

def add_model_to_snpeff(model_accession, script_dir):
    """
    Add model virus to SnpEff using Pathway 2.

    Args:
        model_accession: NCBI accession for model
        script_dir: Directory containing VICAST scripts

    Returns:
        True if successful, False otherwise
    """
    print("\n" + "="*60)
    print("Adding model virus to SnpEff (Pathway 2)")
    print("="*60)

    # Find step1 and step2 scripts
    step1_script = os.path.join(script_dir, 'step1_parse_viral_genome.py')
    step2_script = os.path.join(script_dir, 'step2_add_to_snpeff.py')

    if not os.path.exists(step1_script) or not os.path.exists(step2_script):
        print(f"  ✗ Cannot find Pathway 2 scripts in {script_dir}")
        return False

    try:
        # Run step1
        print(f"\nRunning step1 on {model_accession}...")
        result = subprocess.run(
            ['python3', step1_script, model_accession],
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            print(f"  ✗ Step1 failed:")
            print(result.stderr)
            return False

        print("  ✓ Step1 completed")

        # Find the TSV file created by step1
        tsv_file = f"{model_accession}_no_polyprotein.tsv"
        if not os.path.exists(tsv_file):
            print(f"  ✗ Expected TSV not found: {tsv_file}")
            return False

        # Run step2
        print(f"\nRunning step2 on {model_accession}...")
        result = subprocess.run(
            ['python3', step2_script, model_accession, tsv_file],
            capture_output=True,
            text=True,
            timeout=300
        )

        if result.returncode != 0:
            print(f"  ✗ Step2 failed:")
            print(result.stderr)
            return False

        print("  ✓ Step2 completed")
        print(f"  ✓ Model {model_accession} added to SnpEff!")

        return True

    except subprocess.TimeoutExpired:
        print("  ✗ Pathway 2 timed out")
        return False
    except Exception as e:
        print(f"  ✗ Error running Pathway 2: {e}")
        return False

#=============================================================================
# BLAST DATABASE FUNCTIONS
#=============================================================================

def create_blast_database(protein_fasta, output_name='model_proteins'):
    """
    Create BLAST database from protein sequences.

    Args:
        protein_fasta: Path to protein FASTA file
        output_name: Name for BLAST database

    Returns:
        Path to BLAST database (basename), or None if failed
    """
    print(f"\nCreating BLAST database from {protein_fasta}...")

    try:
        cmd = [
            'makeblastdb',
            '-in', protein_fasta,
            '-dbtype', 'prot',
            '-out', output_name,
            '-parse_seqids'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        if result.returncode == 0:
            print(f"  ✓ BLAST database created: {output_name}")
            return output_name
        else:
            print(f"  ✗ makeblastdb failed:")
            print(result.stderr)
            return None

    except FileNotFoundError:
        print("  ✗ makeblastdb not found. Please install BLAST+")
        return None
    except Exception as e:
        print(f"  ✗ Error creating BLAST database: {e}")
        return None

#=============================================================================
# BLASTX FUNCTIONS
#=============================================================================

def run_blastx(query_fasta, blast_db, evalue=0.05, max_target_seqs=100):
    """
    Run BLASTx against protein database.

    Args:
        query_fasta: Subject genome FASTA
        blast_db: BLAST database path
        evalue: E-value threshold (default 0.05, optimized for single-genome annotation transfer)
        max_target_seqs: Maximum hits per query (default 100, ensures all model proteins are considered)

    Returns:
        Path to BLAST output file, or None if failed

    Note:
        Parameters optimized for viral genome annotation transfer (not database searching):
        - E-value 0.05: More permissive to find divergent homologs in related strains
        - Word size 5: Balanced for normal-sized viral proteins (50-800 aa)
        - Gap costs 11,1: Standard for BLOSUM62
        - Window size 40: Composition-based statistics

        Since we're searching ONE viral genome (not millions), we prioritize:
        - Hit length over stringency
        - Finding full-length proteins over perfect short matches
        - Sensitivity over specificity
    """
    output_file = f"{Path(query_fasta).stem}_blastx.txt"

    print(f"\nRunning BLASTx annotation...")
    print(f"  Query: {query_fasta}")
    print(f"  Database: {blast_db}")
    print(f"  E-value: {evalue} (optimized for single-genome annotation transfer)")
    print(f"  Word size: 5, Gap costs: 11,1, Window: 40")
    print("  This may take several minutes...")

    try:
        cmd = [
            'blastx',
            '-query', query_fasta,
            '-db', blast_db,
            '-evalue', str(evalue),
            '-word_size', '5',              # Balanced for normal viral proteins (50-800 aa)
            '-gapopen', '11',               # Gap open penalty (BLOSUM62 standard)
            '-gapextend', '1',              # Gap extension penalty
            '-window_size', '40',           # Window size for composition-based stats
            '-seg', 'no',                   # Turn off low-complexity filtering
            '-max_target_seqs', str(max_target_seqs),
            '-max_hsps', '20',              # Allow more HSPs per hit
            '-outfmt', '6 qseqid qstart qend sseqid pident evalue stitle qlen slen',
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
        print("  ✗ blastx not found. Please install BLAST+ or add to PATH:")
        print("     export PATH=\"/ref/sahlab/software/ncbi-blast-2.10.0+/bin:$PATH\"")
        return None
    except Exception as e:
        print(f"  ✗ Error running BLASTx: {e}")
        return None

def parse_blastx_results(blast_file):
    """
    Parse BLASTx tabular output.

    Args:
        blast_file: Path to BLAST output

    Returns:
        List of hit dictionaries
    """
    hits = []

    print(f"\nParsing BLASTx results...")

    try:
        with open(blast_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue

                # Parse fields - now includes qlen and slen at end
                qseqid, qstart, qend, sseqid, pident, evalue, stitle = parts[:7]
                qlen = parts[7] if len(parts) > 7 else None
                slen = parts[8] if len(parts) > 8 else None

                # Determine strand
                qstart, qend = int(qstart), int(qend)
                strand = '+' if qstart < qend else '-'
                if strand == '-':
                    qstart, qend = qend, qstart

                # Extract gene/product from title
                gene = extract_gene_from_title(stitle)
                product = extract_product_from_title(stitle)

                hit_dict = {
                    'query': qseqid,
                    'start': qstart,
                    'end': qend,
                    'strand': strand,
                    'hit_id': sseqid,
                    'identity': float(pident),
                    'evalue': float(evalue),
                    'gene': gene,
                    'product': product,
                    'hit_title': stitle[:100]
                }

                # Add optional fields if available
                if qlen:
                    hit_dict['qlen'] = int(qlen)
                if slen:
                    hit_dict['slen'] = int(slen)

                hits.append(hit_dict)

        print(f"  Found {len(hits)} BLAST hits")
        return hits

    except Exception as e:
        print(f"  ✗ Error parsing BLAST results: {e}")
        return []

def extract_gene_from_title(title):
    """Extract gene name from BLAST hit title."""
    # Pattern 1: gene=XXXX in title
    match = re.search(r'gene[=:]([A-Za-z0-9_-]+)', title, re.IGNORECASE)
    if match:
        return match.group(1)

    # Pattern 2: Common viral gene patterns (NS1, VP1, etc.)
    match = re.search(r'\b(NS\d+[A-Z]?|VP\d+|[EML]\d+)\b', title)
    if match:
        return match.group(1)

    return ''

def extract_product_from_title(title):
    """Extract product description from BLAST hit title."""
    # Remove accession at start
    match = re.search(r'\w+\.\d+\s+(.+)', title)
    if match:
        product = match.group(1).strip()
    else:
        product = title.strip()

    # Remove organism info in brackets at end
    product = re.sub(r'\s*\[.*?\]\s*$', '', product)

    return product[:100] if product else 'hypothetical protein'

def merge_overlapping_hits(hits, overlap_threshold=0.5):
    """
    Merge overlapping BLAST hits, keeping the best.

    Args:
        hits: List of hit dictionaries
        overlap_threshold: Minimum overlap fraction to merge

    Returns:
        Filtered list of non-overlapping hits
    """
    if not hits:
        return []

    print("\nMerging overlapping hits...")

    # Sort by start position, then by identity (descending)
    sorted_hits = sorted(hits, key=lambda x: (x['start'], -x['identity']))

    merged = []
    current = sorted_hits[0].copy()

    for hit in sorted_hits[1:]:
        # Calculate overlap
        overlap_start = max(current['start'], hit['start'])
        overlap_end = min(current['end'], hit['end'])
        overlap_len = max(0, overlap_end - overlap_start)

        current_len = current['end'] - current['start']
        hit_len = hit['end'] - hit['start']

        if current_len > 0 and hit_len > 0:
            overlap_frac = overlap_len / min(current_len, hit_len)
        else:
            overlap_frac = 0

        if overlap_frac >= overlap_threshold:
            # Overlapping - keep the better hit
            if hit['identity'] > current['identity']:
                current = hit.copy()
        else:
            # No significant overlap - save current and move to next
            merged.append(current)
            current = hit.copy()

    # Add the last one
    merged.append(current)

    print(f"  Merged to {len(merged)} non-overlapping features")
    return merged

#=============================================================================
# TSV CREATION
#=============================================================================

def create_annotation_tsv(hits, output_tsv, seqid):
    """
    Create TSV file from BLAST hits for manual curation.

    Args:
        hits: List of hit dictionaries
        output_tsv: Output TSV path
        seqid: Sequence ID

    Returns:
        Output TSV path
    """
    print(f"\nCreating annotation TSV...")

    data = []

    for i, hit in enumerate(hits, 1):
        data.append({
            'action': 'KEEP',
            'seqid': seqid,
            'source': 'BLASTx',
            'type': 'CDS',
            'start': hit['start'],
            'end': hit['end'],
            'strand': hit['strand'],
            'gene_name': hit.get('gene', f'ORF{i}'),
            'product': hit['product'],
            'ID': f"cds_{i}",
            'gene': hit.get('gene', ''),
            'protein_id': '',
            'notes': f"Model:{hit['hit_id']} E={hit['evalue']:.2e} ID={hit['identity']:.1f}%"
        })

    df = pd.DataFrame(data)
    df.to_csv(output_tsv, sep='\t', index=False)

    print(f"  ✓ Created: {output_tsv}")
    print(f"  Total features: {len(df)}")

    return output_tsv


def generate_protein_fasta(hits, subject_fasta, output_fasta, seqid):
    """
    Generate protein FASTA file from BLASTx hits.

    Args:
        hits: List of hit dictionaries with genomic coordinates
        subject_fasta: Path to subject genome FASTA
        output_fasta: Output protein FASTA path
        seqid: Sequence ID

    Returns:
        Output FASTA path
    """
    from Bio import SeqIO
    from Bio.Seq import Seq

    print(f"\nGenerating protein FASTA...")

    # Read subject genome
    genome_seq = None
    for record in SeqIO.parse(subject_fasta, "fasta"):
        if record.id == seqid or record.id.split()[0] == seqid:
            genome_seq = str(record.seq)
            break

    if not genome_seq:
        print(f"  ⚠ Warning: Could not find sequence {seqid} in {subject_fasta}")
        return None

    # Generate protein sequences
    proteins = []
    for i, hit in enumerate(hits, 1):
        start = hit['start'] - 1  # Convert to 0-based
        end = hit['end']
        strand = hit['strand']

        # Extract nucleotide sequence
        nt_seq = genome_seq[start:end]

        # Reverse complement if minus strand
        if strand == '-':
            nt_seq = str(Seq(nt_seq).reverse_complement())

        # Translate
        try:
            protein_seq = str(Seq(nt_seq).translate(to_stop=False))

            # Create FASTA header
            gene_name = hit.get('gene', f'ORF{i}')
            product = hit['product']
            header = f">{seqid}_{gene_name} {product} [BLASTx:{hit['hit_id']} E={hit['evalue']:.2e} ID={hit['identity']:.1f}%] [{start+1}:{end}({strand})]"

            proteins.append(f"{header}\n{protein_seq}")

        except Exception as e:
            print(f"  ⚠ Warning: Could not translate {gene_name}: {e}")
            continue

    # Write to file
    with open(output_fasta, 'w') as f:
        f.write('\n'.join(proteins) + '\n')

    print(f"  ✓ Created: {output_fasta}")
    print(f"  Total proteins: {len(proteins)}")

    return output_fasta


#=============================================================================
# MAIN WORKFLOW
#=============================================================================

def main():
    parser = argparse.ArgumentParser(
        description='STEP 1 (Pathway 3): Model-based homology annotation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Use this pathway when NCBI has poor/no annotations for your virus.

THREE MODES OF OPERATION:

1. Model already in SnpEff:
   python3 step1_blastx_annotate.py --subject NC_XXXXXX --model NC_001477

2. Add model to SnpEff first (auto-runs Pathway 2):
   python3 step1_blastx_annotate.py --subject NC_XXXXXX --model NC_001477 --add-model

3. Use custom model proteins:
   python3 step1_blastx_annotate.py --subject NC_XXXXXX.fasta --model-proteins my_proteins.fa

WORKFLOW:
  1. Check/prepare model virus
  2. Extract protein sequences from model
  3. Create BLAST database
  4. Download subject genome (if accession provided)
  5. Run BLASTx: subject → model proteins
  6. Transfer annotations from model to subject
  7. Create TSV for manual curation

NEXT STEPS:
  After reviewing TSV, run step2 to add to SnpEff:
    python3 step2_add_to_snpeff.py subject_id output.tsv
        """
    )

    # Subject genome (to be annotated)
    parser.add_argument('--subject', required=True,
                       help='Subject genome to annotate (NCBI accession or FASTA file)')

    # Model virus (three modes)
    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument('--model',
                            help='Model virus ID (NCBI accession, e.g., NC_001477)')
    model_group.add_argument('--model-proteins',
                            help='Custom model protein FASTA file')

    # Model handling
    parser.add_argument('--add-model', action='store_true',
                       help='Add model to SnpEff first (auto-run Pathway 2)')

    # BLAST parameters
    parser.add_argument('--evalue', type=float, default=0.05,
                       help='E-value threshold (default: 0.05, optimized for single-genome annotation transfer)')
    parser.add_argument('--max-target-seqs', type=int, default=100,
                       help='Maximum BLAST hits to report (default: 100, ensures all model proteins found)')
    parser.add_argument('--overlap-threshold', type=float, default=0.5,
                       help='Overlap threshold for merging hits (default: 0.5)')
    parser.add_argument('--no-merge', action='store_true',
                       help='Skip merging overlapping hits (keep all hits)')

    # Output
    parser.add_argument('--output', '-o',
                       help='Output base name (default: from subject filename)')

    args = parser.parse_args()

    print("="*60)
    print("STEP 1 (Pathway 3): Model-based Homology Annotation")
    print("="*60)

    # Get script directory for Pathway 2 integration
    script_dir = os.path.dirname(os.path.abspath(__file__))

    #=========================================================================
    # STEP 1: Handle subject genome
    #=========================================================================

    print("\n" + "="*60)
    print("STEP 1: Preparing subject genome")
    print("="*60)

    # Check if subject is file or accession
    if os.path.exists(args.subject):
        subject_fasta = args.subject
        subject_id = Path(subject_fasta).stem
        print(f"  ✓ Using subject FASTA: {subject_fasta}")
    else:
        # Assume it's an accession - download it
        print(f"  Subject: {args.subject} (NCBI accession)")
        gb_file, subject_fasta = download_genome(args.subject)
        if not subject_fasta:
            print("\n✗ Failed to download subject genome")
            sys.exit(1)
        subject_id = args.subject

    # Get sequence ID from FASTA
    try:
        for record in SeqIO.parse(subject_fasta, "fasta"):
            seqid = record.id
            genome_length = len(record.seq)
            print(f"  Sequence ID: {seqid}")
            print(f"  Length: {genome_length:,} bp")
            break
    except Exception as e:
        print(f"  ✗ Error reading FASTA: {e}")
        sys.exit(1)

    #=========================================================================
    # STEP 2: Handle model virus
    #=========================================================================

    print("\n" + "="*60)
    print("STEP 2: Preparing model virus")
    print("="*60)

    if args.model_proteins:
        # Mode 3: Custom proteins
        print(f"  Mode: Custom model proteins")
        if not os.path.exists(args.model_proteins):
            print(f"  ✗ Custom protein file not found: {args.model_proteins}")
            sys.exit(1)

        model_proteins = args.model_proteins
        print(f"  ✓ Using custom proteins: {model_proteins}")

    else:
        # Mode 1 or 2: Model from SnpEff
        model_id = args.model
        print(f"  Mode: Model virus from SnpEff")
        print(f"  Model: {model_id}")

        # Check if model is in SnpEff
        is_in_snpeff, model_dir = check_model_in_snpeff(model_id)

        if is_in_snpeff:
            print(f"  ✓ Model found in SnpEff: {model_dir}")
        elif args.add_model:
            # Mode 2: Add model first
            print(f"  Model not in SnpEff - will add using Pathway 2")
            success = add_model_to_snpeff(model_id, script_dir)
            if not success:
                print("\n✗ Failed to add model to SnpEff")
                sys.exit(1)
        else:
            print(f"  ✗ Model {model_id} not in SnpEff")
            print("\n  Use --add-model to automatically add it:")
            print(f"    python3 {sys.argv[0]} --subject {args.subject} --model {model_id} --add-model")
            sys.exit(1)

        # Extract proteins from SnpEff
        model_proteins = extract_proteins_from_snpeff(model_id)
        if not model_proteins:
            print(f"\n✗ Could not extract proteins from model {model_id}")
            sys.exit(1)

    #=========================================================================
    # STEP 3: Create BLAST database
    #=========================================================================

    print("\n" + "="*60)
    print("STEP 3: Creating BLAST database")
    print("="*60)

    blast_db = create_blast_database(model_proteins, 'model_proteins')
    if not blast_db:
        print("\n✗ Failed to create BLAST database")
        sys.exit(1)

    #=========================================================================
    # STEP 4: Run BLASTx
    #=========================================================================

    print("\n" + "="*60)
    print("STEP 4: Running BLASTx annotation")
    print("="*60)

    blast_output = run_blastx(subject_fasta, blast_db, args.evalue, args.max_target_seqs)
    if not blast_output:
        print("\n✗ BLASTx failed")
        sys.exit(1)

    #=========================================================================
    # STEP 5: Parse and merge hits
    #=========================================================================

    print("\n" + "="*60)
    print("STEP 5: Processing BLAST hits")
    print("="*60)

    hits = parse_blastx_results(blast_output)

    if not hits:
        print("\n⚠ WARNING: No BLAST hits found!")
        print("\nPossible reasons:")
        print("  - Subject and model are not related")
        print("  - E-value threshold too stringent (try --evalue 1e-2 or higher)")
        print("  - Subject genome is truly novel")
        sys.exit(1)

    # Merge overlapping hits (unless --no-merge specified)
    if args.no_merge:
        print("\nSkipping merge (--no-merge specified)")
        print(f"  Keeping all {len(hits)} hits")
        merged_hits = hits
    else:
        merged_hits = merge_overlapping_hits(hits, args.overlap_threshold)

    #=========================================================================
    # STEP 6: Create annotation TSV
    #=========================================================================

    print("\n" + "="*60)
    print("STEP 6: Creating annotation TSV")
    print("="*60)

    output_base = args.output if args.output else subject_id
    output_tsv = f"{output_base}_blastx.tsv"
    output_fasta = f"{output_base}_blastx_proteins.faa"

    create_annotation_tsv(merged_hits, output_tsv, seqid)

    # Generate protein FASTA from hits
    generate_protein_fasta(merged_hits, subject_fasta, output_fasta, seqid)

    #=========================================================================
    # DONE
    #=========================================================================

    print("\n" + "="*60)
    print("✓ PATHWAY 3 COMPLETE")
    print("="*60)

    print(f"\nAnnotated {len(merged_hits)} features using model-based homology")

    print("\n" + "="*60)
    print("NEXT STEPS:")
    print("="*60)
    print(f"\n1. REVIEW the annotation files:")
    print(f"   TSV:     {output_tsv}")
    print(f"   Proteins: {output_fasta}")
    print("\n   Use the protein FASTA to:")
    print("   - Verify protein sequences look correct")
    print("   - Check for premature stops (indicated by * in sequence)")
    print("   - BLAST against NCBI to confirm identities")
    print("   - Assess annotation quality")

    print(f"\n2. EDIT the TSV file if needed:")
    print("   - Verify gene boundaries are correct")
    print("   - Check for missed features")
    print("   - Refine gene names and products")
    print("   - Remove any false positives")
    print("   - Mark frameshifts if present")

    print(f"\n3. SAVE your edited file (keep as TSV)")

    print(f"\n4. RUN STEP 2 to add to SnpEff:")
    print(f"   python3 step2_add_to_snpeff.py {subject_id} {output_tsv}")

    print("\n" + "="*60)

    # Cleanup BLAST database files
    for ext in ['.phr', '.pin', '.psq', '.pdb', '.pot', '.ptf', '.pto']:
        db_file = f"model_proteins{ext}"
        if os.path.exists(db_file):
            os.remove(db_file)

if __name__ == '__main__':
    main()
