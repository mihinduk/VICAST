#!/usr/bin/env python3
"""
Check Variant Co-occurrence in BAM Files
==========================================

This module provides read-level phasing analysis for viral variants, enabling
TRUE haplotype reconstruction by checking if variants appear on the same reads.

Scientific Context
------------------
Current limitation of frequency-based haplotype inference:
- Frequency-based methods (like generate_realistic_haplotype_consensus.py) can
  infer that high-frequency variants likely co-occur based on statistics
- However, they CANNOT prove that mutations actually appear together on the
  same viral genomes without examining sequencing reads

This tool provides:
- Direct evidence of variant co-occurrence on the same reads/read pairs
- Read-level phasing for variants within sequencing distance (~150bp read length
  or ~500bp insert size for paired-end)
- Quantitative co-occurrence statistics (not just inference)

Use Cases
---------
1. Validate frequency-based haplotype predictions
2. Prove linkage for nearby variants (< insert size)
3. Identify recombination events (variants that don't co-occur as expected)
4. Distinguish true haplotypes from technical artifacts

Limitations
-----------
- Only works for variants within read length or insert size
- Cannot phase distant variants (>500bp for typical paired-end sequencing)
- Requires sufficient coverage at both variant positions
- Low-frequency variants may have insufficient reads for statistical confidence

Usage Examples
--------------
Basic usage (check all variant pairs in VCF):
    python check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf

With custom maximum distance:
    python check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf --max-distance 300

Filter by minimum coverage:
    python check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf --min-coverage 100

Check only high-confidence variants:
    python check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf --min-qual 1000 --min-depth 200

Output to specific file:
    python check_read_cooccurrence.py --bam aligned.bam --vcf variants.vcf --output cooccurrence.tsv

Requirements
------------
- pysam: BAM file parsing
- pandas: VCF parsing and data manipulation
- numpy: Statistical calculations (optional, for future enhancements)

Author: Claude/VICAST Team
Date: 2026-02-05
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Optional, Set
import logging

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Check for required dependencies
try:
    import pysam
except ImportError:
    logger.error("pysam is not installed. Install with: conda install -c bioconda pysam")
    logger.error("Or: pip install pysam")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    logger.error("pandas is not installed. Install with: conda install pandas")
    sys.exit(1)


class VariantCooccurrenceAnalyzer:
    """
    Analyzer for checking variant co-occurrence in BAM files.

    This class provides methods to:
    1. Parse BAM files to extract reads spanning multiple variant positions
    2. Determine which allele (reference or alternate) is present at each position
    3. Calculate co-occurrence statistics for variant pairs
    4. Report proven linkage vs. inferred linkage
    """

    def __init__(self, bam_file: str, min_base_quality: int = 20, min_mapping_quality: int = 20):
        """
        Initialize the analyzer.

        Parameters
        ----------
        bam_file : str
            Path to sorted, indexed BAM file
        min_base_quality : int
            Minimum base quality to consider (default: 20)
        min_mapping_quality : int
            Minimum mapping quality to consider (default: 20)
        """
        self.bam_file = Path(bam_file)
        self.min_base_quality = min_base_quality
        self.min_mapping_quality = min_mapping_quality

        # Validate BAM file
        if not self.bam_file.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_file}")

        # Check for index
        index_file = Path(str(self.bam_file) + ".bai")
        if not index_file.exists():
            logger.warning(f"BAM index not found at {index_file}")
            logger.info("Attempting to create index...")
            try:
                pysam.index(str(self.bam_file))
                logger.info("Successfully created BAM index")
            except Exception as e:
                raise RuntimeError(f"Failed to index BAM file: {e}")

        # Open BAM file
        try:
            self.bamfile = pysam.AlignmentFile(str(self.bam_file), "rb")
        except Exception as e:
            raise RuntimeError(f"Failed to open BAM file: {e}")

        # Get reference names
        self.references = self.bamfile.references
        if len(self.references) == 0:
            raise ValueError("BAM file contains no reference sequences")

        logger.info(f"Opened BAM file: {self.bam_file}")
        logger.info(f"Reference sequences: {', '.join(self.references)}")

    def __del__(self):
        """Close BAM file when object is destroyed."""
        if hasattr(self, 'bamfile'):
            self.bamfile.close()

    def get_allele_at_position(self, read, position: int, ref_allele: str, alt_allele: str) -> Optional[str]:
        """
        Determine which allele (ref or alt) is present at a position in a read.

        Parameters
        ----------
        read : pysam.AlignedSegment
            The read to analyze
        position : int
            Genomic position (0-based)
        ref_allele : str
            Reference allele
        alt_allele : str
            Alternate allele

        Returns
        -------
        str or None
            'ref' if reference allele is present
            'alt' if alternate allele is present
            None if position is not covered or quality is too low

        Notes
        -----
        This function handles:
        - Matches and mismatches
        - Insertions and deletions
        - Soft clips (ignored)
        - Base quality filtering
        """
        # Check if read is mapped
        if read.is_unmapped:
            return None

        # Check mapping quality
        if read.mapping_quality < self.min_mapping_quality:
            return None

        # Check if position is covered by this read
        if position < read.reference_start or position >= read.reference_end:
            return None

        # Get aligned pairs (read position, reference position)
        # This handles insertions, deletions, and soft clips
        try:
            aligned_pairs = read.get_aligned_pairs(matches_only=False, with_seq=False)
        except Exception as e:
            logger.debug(f"Failed to get aligned pairs for read {read.query_name}: {e}")
            return None

        # Find the query position corresponding to this reference position
        query_pos = None
        for qpos, rpos in aligned_pairs:
            if rpos == position:
                query_pos = qpos
                break

        if query_pos is None:
            # Position is in a deletion in this read
            if len(ref_allele) > len(alt_allele):
                # This is a deletion variant - check if this read has the deletion
                # For simplicity, we'll consider this as missing data for now
                return None
            return None

        # Check base quality
        if read.query_qualities[query_pos] < self.min_base_quality:
            return None

        # Get the base at this position
        base = read.query_sequence[query_pos]

        # Compare to ref and alt alleles
        # Handle SNPs (most common case)
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            if base == ref_allele:
                return 'ref'
            elif base == alt_allele:
                return 'alt'
            else:
                # Base doesn't match either - could be another variant or sequencing error
                return None

        # Handle insertions and deletions (more complex)
        # For now, we'll use a simplified approach
        # TODO: Implement proper indel handling
        if len(ref_allele) != len(alt_allele):
            logger.debug(f"Indel at position {position}: {ref_allele}>{alt_allele} - simplified handling")
            # For deletions (ref longer than alt), check if base matches ref
            # For insertions (alt longer than ref), check if base matches ref at this position
            if base == ref_allele[0]:
                return 'ref'
            elif base == alt_allele[0]:
                return 'alt'
            else:
                return None

        return None

    def get_spanning_reads(self, chrom: str, pos1: int, pos2: int,
                          use_pairs: bool = True) -> List:
        """
        Get reads that span both variant positions.

        Parameters
        ----------
        chrom : str
            Chromosome/reference name
        pos1 : int
            First variant position (0-based)
        pos2 : int
            Second variant position (0-based)
        use_pairs : bool
            If True, consider read pairs (use insert size for distance)
            If False, only use individual reads (use read length for distance)

        Returns
        -------
        list of pysam.AlignedSegment
            Reads that span both positions
        """
        spanning_reads = []

        # Fetch reads overlapping the region
        min_pos = min(pos1, pos2)
        max_pos = max(pos1, pos2)

        try:
            for read in self.bamfile.fetch(chrom, min_pos, max_pos + 1):
                # Skip unmapped reads
                if read.is_unmapped:
                    continue

                # Skip reads with low mapping quality
                if read.mapping_quality < self.min_mapping_quality:
                    continue

                # Skip secondary and supplementary alignments
                if read.is_secondary or read.is_supplementary:
                    continue

                # Check if this read spans both positions
                if use_pairs and read.is_paired and read.is_proper_pair:
                    # For paired reads, check if the insert spans both positions
                    # Insert is from leftmost position to rightmost position of the pair
                    leftmost = min(read.reference_start, read.next_reference_start)
                    rightmost = max(read.reference_end,
                                   read.next_reference_start + abs(read.template_length))

                    if leftmost <= min_pos and rightmost >= max_pos:
                        spanning_reads.append(read)
                else:
                    # For single reads, check if the read itself spans both positions
                    if read.reference_start <= min_pos and read.reference_end >= max_pos:
                        spanning_reads.append(read)

        except Exception as e:
            logger.error(f"Error fetching reads for {chrom}:{min_pos}-{max_pos}: {e}")
            return []

        return spanning_reads

    def analyze_variant_pair(self, chrom: str, variant1: Dict, variant2: Dict,
                            use_pairs: bool = True) -> Dict:
        """
        Analyze co-occurrence of two variants.

        Parameters
        ----------
        chrom : str
            Chromosome/reference name
        variant1 : dict
            First variant: {'pos': int (0-based), 'ref': str, 'alt': str, 'id': str}
        variant2 : dict
            Second variant: {'pos': int (0-based), 'ref': str, 'alt': str, 'id': str}
        use_pairs : bool
            Whether to use read pairs for distant variants

        Returns
        -------
        dict
            Co-occurrence statistics including:
            - variant_pair: tuple of variant IDs
            - distance: distance between variants (bp)
            - total_spanning_reads: number of reads spanning both positions
            - both_alt: reads with both alternate alleles
            - variant1_only: reads with only variant1 alternate
            - variant2_only: reads with only variant2 alternate
            - both_ref: reads with both reference alleles
            - undetermined: reads where alleles couldn't be determined
            - cooccurrence_rate: proportion of spanning reads with both variants
            - linkage_proven: whether we have evidence for linkage
            - evidence_level: 'read-level' or 'insert-level'
        """
        pos1 = variant1['pos']
        pos2 = variant2['pos']
        distance = abs(pos2 - pos1)

        # Get reads spanning both positions
        spanning_reads = self.get_spanning_reads(chrom, pos1, pos2, use_pairs)

        if len(spanning_reads) == 0:
            return {
                'variant_pair': (variant1['id'], variant2['id']),
                'distance': distance,
                'total_spanning_reads': 0,
                'both_alt': 0,
                'variant1_only': 0,
                'variant2_only': 0,
                'both_ref': 0,
                'undetermined': 0,
                'cooccurrence_rate': None,
                'linkage_proven': False,
                'evidence_level': 'none',
                'message': 'No spanning reads found'
            }

        # Classify each read
        both_alt = 0
        variant1_only = 0
        variant2_only = 0
        both_ref = 0
        undetermined = 0

        for read in spanning_reads:
            allele1 = self.get_allele_at_position(read, pos1, variant1['ref'], variant1['alt'])
            allele2 = self.get_allele_at_position(read, pos2, variant2['ref'], variant2['alt'])

            if allele1 is None or allele2 is None:
                undetermined += 1
            elif allele1 == 'alt' and allele2 == 'alt':
                both_alt += 1
            elif allele1 == 'alt' and allele2 == 'ref':
                variant1_only += 1
            elif allele1 == 'ref' and allele2 == 'alt':
                variant2_only += 1
            elif allele1 == 'ref' and allele2 == 'ref':
                both_ref += 1

        # Calculate co-occurrence rate
        # Only count reads where we could determine alleles at both positions
        informative_reads = both_alt + variant1_only + variant2_only + both_ref

        if informative_reads == 0:
            cooccurrence_rate = None
            linkage_proven = False
        else:
            cooccurrence_rate = both_alt / informative_reads
            linkage_proven = both_alt > 0  # We have direct evidence of co-occurrence

        # Determine evidence level
        if distance <= 150:  # Typical read length
            evidence_level = 'read-level'
        elif use_pairs and distance <= 500:  # Typical insert size
            evidence_level = 'insert-level'
        else:
            evidence_level = 'fragment-level'

        return {
            'variant_pair': (variant1['id'], variant2['id']),
            'variant1_pos': pos1 + 1,  # Convert to 1-based for output
            'variant2_pos': pos2 + 1,
            'distance': distance,
            'total_spanning_reads': len(spanning_reads),
            'informative_reads': informative_reads,
            'both_alt': both_alt,
            'variant1_only': variant1_only,
            'variant2_only': variant2_only,
            'both_ref': both_ref,
            'undetermined': undetermined,
            'cooccurrence_rate': cooccurrence_rate,
            'linkage_proven': linkage_proven,
            'evidence_level': evidence_level
        }

    def analyze_all_variant_pairs(self, variants: List[Dict], max_distance: int = 500,
                                  use_pairs: bool = True) -> List[Dict]:
        """
        Analyze co-occurrence for all variant pairs within maximum distance.

        Parameters
        ----------
        variants : list of dict
            List of variants, each with 'chrom', 'pos' (0-based), 'ref', 'alt', 'id'
        max_distance : int
            Maximum distance between variants to check (default: 500)
        use_pairs : bool
            Whether to use read pairs for distant variants

        Returns
        -------
        list of dict
            Co-occurrence statistics for all variant pairs
        """
        results = []

        # Group variants by chromosome
        variants_by_chrom = defaultdict(list)
        for variant in variants:
            variants_by_chrom[variant['chrom']].append(variant)

        # Sort variants within each chromosome
        for chrom in variants_by_chrom:
            variants_by_chrom[chrom].sort(key=lambda v: v['pos'])

        # Analyze pairs within each chromosome
        for chrom, chrom_variants in variants_by_chrom.items():
            logger.info(f"Analyzing {len(chrom_variants)} variants on {chrom}")

            for i in range(len(chrom_variants)):
                for j in range(i + 1, len(chrom_variants)):
                    variant1 = chrom_variants[i]
                    variant2 = chrom_variants[j]

                    distance = variant2['pos'] - variant1['pos']

                    # Skip if too far apart
                    if distance > max_distance:
                        break  # No need to check further variants

                    # Analyze this pair
                    result = self.analyze_variant_pair(chrom, variant1, variant2, use_pairs)
                    results.append(result)

                    # Log progress
                    if result['linkage_proven']:
                        logger.info(f"  {variant1['id']} + {variant2['id']}: "
                                  f"{result['both_alt']}/{result['informative_reads']} reads "
                                  f"({result['cooccurrence_rate']:.2%}) - LINKED")

        return results


def parse_vcf(vcf_file: str, min_qual: float = 0, min_depth: int = 0,
              min_freq: float = 0) -> List[Dict]:
    """
    Parse VCF file and extract variant information.

    Parameters
    ----------
    vcf_file : str
        Path to VCF or TSV file
    min_qual : float
        Minimum quality score to include
    min_depth : int
        Minimum depth to include
    min_freq : float
        Minimum allele frequency to include

    Returns
    -------
    list of dict
        List of variants with 'chrom', 'pos' (0-based), 'ref', 'alt', 'id'
    """
    logger.info(f"Parsing variants from {vcf_file}")

    # Try to read as TSV (VICAST format) or VCF
    try:
        with open(vcf_file, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('#CHROM') or 'POS' in first_line:
                # TSV format (VICAST filtered output)
                df = pd.read_csv(vcf_file, sep='\t', comment='#')

                # Apply filters if columns exist
                if 'QUAL' in df.columns:
                    df = df[df['QUAL'] >= min_qual]
                if 'Total_Depth' in df.columns:
                    df = df[df['Total_Depth'] >= min_depth]
                if 'Allele_Frequency' in df.columns:
                    df = df[df['Allele_Frequency'] >= min_freq]

                # Convert to variant list
                variants = []
                for _, row in df.iterrows():
                    chrom = str(row.get('CHROM', row.get('#CHROM', 'unknown')))
                    pos = int(row['POS']) - 1  # Convert to 0-based
                    ref = str(row['REF'])
                    alt = str(row['ALT'])
                    var_id = f"{row['POS']}{ref}>{alt}"

                    variants.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'id': var_id
                    })

                logger.info(f"Loaded {len(variants)} variants from TSV")
                return variants

            else:
                # Try VCF format
                logger.info("Attempting to parse as VCF format")
                # Re-open and skip header lines
                with open(vcf_file, 'r') as f:
                    lines = [line for line in f if not line.startswith('##')]

                df = pd.read_csv(pd.io.common.StringIO(''.join(lines)), sep='\t')

                # VCF format
                variants = []
                for _, row in df.iterrows():
                    chrom = str(row['#CHROM'])
                    pos = int(row['POS']) - 1  # Convert to 0-based
                    ref = str(row['REF'])
                    alt = str(row['ALT'])
                    qual = float(row['QUAL']) if row['QUAL'] != '.' else 0

                    # Apply quality filter
                    if qual < min_qual:
                        continue

                    var_id = f"{row['POS']}{ref}>{alt}"

                    variants.append({
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'id': var_id
                    })

                logger.info(f"Loaded {len(variants)} variants from VCF")
                return variants

    except Exception as e:
        logger.error(f"Failed to parse VCF file: {e}")
        raise


def write_results(results: List[Dict], output_file: str):
    """
    Write co-occurrence results to TSV file.

    Parameters
    ----------
    results : list of dict
        Co-occurrence statistics
    output_file : str
        Output file path
    """
    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Reorder columns for better readability
    column_order = [
        'variant1_pos', 'variant2_pos', 'distance', 'evidence_level',
        'total_spanning_reads', 'informative_reads',
        'both_alt', 'variant1_only', 'variant2_only', 'both_ref', 'undetermined',
        'cooccurrence_rate', 'linkage_proven'
    ]

    # Only include columns that exist
    column_order = [col for col in column_order if col in df.columns]
    df = df[column_order]

    # Format cooccurrence_rate as percentage
    if 'cooccurrence_rate' in df.columns:
        df['cooccurrence_rate'] = df['cooccurrence_rate'].apply(
            lambda x: f"{x:.4f}" if x is not None else "NA"
        )

    # Write to file
    df.to_csv(output_file, sep='\t', index=False)
    logger.info(f"Results written to {output_file}")


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Check variant co-occurrence in BAM files for read-level phasing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  %(prog)s --bam aligned.bam --vcf variants.vcf

  # With custom distance and quality filters
  %(prog)s --bam aligned.bam --vcf variants.vcf --max-distance 300 --min-qual 1000

  # High-quality variants only
  %(prog)s --bam aligned.bam --vcf variants.vcf --min-qual 1000 --min-depth 200 --min-freq 0.05

  # Use only individual reads (not pairs)
  %(prog)s --bam aligned.bam --vcf variants.vcf --no-pairs

For detailed information, see module docstring or VICAST documentation.
        """
    )

    # Required arguments
    parser.add_argument('--bam', required=True, help='Input BAM file (sorted and indexed)')
    parser.add_argument('--vcf', required=True, help='Input VCF or filtered TSV file')

    # Optional arguments
    parser.add_argument('--output', '-o', help='Output TSV file (default: cooccurrence.tsv)',
                       default='cooccurrence.tsv')
    parser.add_argument('--max-distance', type=int, default=500,
                       help='Maximum distance between variants to check (default: 500)')
    parser.add_argument('--min-qual', type=float, default=0,
                       help='Minimum variant quality score (default: 0)')
    parser.add_argument('--min-depth', type=int, default=0,
                       help='Minimum variant depth (default: 0)')
    parser.add_argument('--min-freq', type=float, default=0,
                       help='Minimum variant allele frequency (default: 0)')
    parser.add_argument('--min-base-quality', type=int, default=20,
                       help='Minimum base quality for allele calling (default: 20)')
    parser.add_argument('--min-mapping-quality', type=int, default=20,
                       help='Minimum mapping quality for reads (default: 20)')
    parser.add_argument('--no-pairs', action='store_true',
                       help='Use only individual reads, not read pairs')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose output')

    args = parser.parse_args()

    # Set logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)

    # Print banner
    logger.info("=" * 80)
    logger.info("VARIANT CO-OCCURRENCE ANALYSIS - Read-Level Phasing")
    logger.info("=" * 80)
    logger.info(f"BAM file: {args.bam}")
    logger.info(f"VCF file: {args.vcf}")
    logger.info(f"Maximum distance: {args.max_distance} bp")
    logger.info(f"Use read pairs: {not args.no_pairs}")
    logger.info(f"Filters: QUAL>={args.min_qual}, Depth>={args.min_depth}, Freq>={args.min_freq}")
    logger.info("=" * 80)

    try:
        # Parse variants
        variants = parse_vcf(args.vcf, args.min_qual, args.min_depth, args.min_freq)

        if len(variants) == 0:
            logger.warning("No variants found in VCF file after filtering")
            sys.exit(0)

        if len(variants) == 1:
            logger.warning("Only one variant found - need at least 2 for co-occurrence analysis")
            sys.exit(0)

        logger.info(f"Analyzing {len(variants)} variants")

        # Initialize analyzer
        analyzer = VariantCooccurrenceAnalyzer(
            args.bam,
            min_base_quality=args.min_base_quality,
            min_mapping_quality=args.min_mapping_quality
        )

        # Analyze all pairs
        results = analyzer.analyze_all_variant_pairs(
            variants,
            max_distance=args.max_distance,
            use_pairs=not args.no_pairs
        )

        if len(results) == 0:
            logger.warning("No variant pairs within maximum distance")
            sys.exit(0)

        # Write results
        write_results(results, args.output)

        # Summary statistics
        logger.info("")
        logger.info("=" * 80)
        logger.info("SUMMARY")
        logger.info("=" * 80)
        logger.info(f"Total variant pairs analyzed: {len(results)}")

        linked_pairs = sum(1 for r in results if r['linkage_proven'])
        logger.info(f"Pairs with proven linkage: {linked_pairs}")

        if linked_pairs > 0:
            avg_cooccurrence = sum(r['cooccurrence_rate'] for r in results
                                  if r['cooccurrence_rate'] is not None) / linked_pairs
            logger.info(f"Average co-occurrence rate: {avg_cooccurrence:.2%}")

        logger.info(f"\nResults saved to: {args.output}")
        logger.info("\nNext steps:")
        logger.info("  1. Review cooccurrence.tsv for proven variant linkages")
        logger.info("  2. Compare with frequency-based haplotype predictions")
        logger.info("  3. Update haplotype consensus with validated linkages")

    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
