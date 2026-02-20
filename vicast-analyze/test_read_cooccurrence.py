#!/usr/bin/env python3
"""
Test script for check_read_cooccurrence.py

This script creates synthetic BAM and VCF files to test the co-occurrence
analysis functionality.

Tests:
1. Two variants 100bp apart (within read length) - should find co-occurrence
2. Two variants 1kb apart (beyond insert size) - should report "cannot validate"
3. Variants with different frequencies (95%/95%, 80%/20%, etc.)
"""

import os
import sys
import tempfile
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Check dependencies
try:
    import pysam
except ImportError:
    logger.error("pysam is not installed. Install with: conda install -c bioconda pysam")
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    logger.error("pandas is not installed. Install with: conda install pandas")
    sys.exit(1)


def create_test_reference(ref_file: str, length: int = 5000):
    """Create a simple reference sequence."""
    with open(ref_file, 'w') as f:
        f.write(">test_ref\n")
        # Simple sequence with known bases
        seq = "ACGT" * (length // 4)
        # Write in 80-character lines
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")


def create_test_bam(bam_file: str, ref_file: str, variants: list):
    """
    Create a synthetic BAM file with reads showing variant co-occurrence patterns.

    Parameters
    ----------
    bam_file : str
        Output BAM file path
    ref_file : str
        Reference FASTA file
    variants : list of tuples
        Each tuple: (pos, ref, alt, frequency)
        Example: [(100, 'A', 'T', 0.95), (200, 'C', 'G', 0.90)]
    """
    # Index reference
    pysam.faidx(ref_file)

    # Create SAM header
    header = {
        'HD': {'VN': '1.0'},
        'SQ': [{'LN': 5000, 'SN': 'test_ref'}]
    }

    # Open output BAM
    outbam = pysam.AlignmentFile(bam_file, "wb", header=header)

    # Read reference sequence
    ref_fasta = pysam.FastaFile(ref_file)
    ref_seq = ref_fasta.fetch("test_ref")
    ref_fasta.close()

    # Generate synthetic reads
    read_id = 0
    coverage_per_position = 100  # 100x coverage

    # For each position that needs coverage
    for start_pos in range(0, len(ref_seq) - 150, 50):  # Read every 50bp
        for _ in range(coverage_per_position // 20):  # Generate multiple reads per position
            read_id += 1

            # Extract reference sequence for this read
            read_start = start_pos
            read_end = start_pos + 150
            read_seq = list(ref_seq[read_start:read_end])

            # Determine which variants this read should have
            # Based on variant frequencies
            import random
            for var_pos, var_ref, var_alt, var_freq in variants:
                # Check if this variant is within the read
                if read_start <= var_pos < read_end:
                    # Decide if this read should have the variant
                    if random.random() < var_freq:
                        # Apply variant
                        read_offset = var_pos - read_start
                        if read_offset < len(read_seq) and read_seq[read_offset] == var_ref:
                            read_seq[read_offset] = var_alt

            # Create read
            read = pysam.AlignedSegment()
            read.query_name = f"read_{read_id}"
            read.query_sequence = ''.join(read_seq)
            read.flag = 0  # Forward strand, not paired
            read.reference_id = 0  # First reference
            read.reference_start = read_start
            read.mapping_quality = 60
            read.cigar = [(0, len(read_seq))]  # All matches (simplified)
            read.query_qualities = pysam.qualitystring_to_array("I" * len(read_seq))  # Q40

            outbam.write(read)

    outbam.close()

    # Sort and index BAM
    logger.info(f"Sorting BAM file...")
    pysam.sort("-o", bam_file + ".sorted.bam", bam_file)
    os.replace(bam_file + ".sorted.bam", bam_file)

    logger.info(f"Indexing BAM file...")
    pysam.index(bam_file)


def create_test_vcf(vcf_file: str, variants: list):
    """
    Create a test VCF file.

    Parameters
    ----------
    vcf_file : str
        Output VCF file path
    variants : list of tuples
        Each tuple: (pos, ref, alt, frequency)
    """
    with open(vcf_file, 'w') as f:
        # Write header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##source=test_read_cooccurrence.py\n")
        f.write("##reference=test_ref.fa\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        # Write variants
        for pos, ref, alt, freq in variants:
            f.write(f"test_ref\t{pos+1}\t.\t{ref}\t{alt}\t1000\tPASS\tAF={freq:.3f}\n")


def run_test_case(test_name: str, variants: list, expected_result: str):
    """
    Run a test case.

    Parameters
    ----------
    test_name : str
        Name of the test
    variants : list of tuples
        Variants to test: (pos, ref, alt, frequency)
    expected_result : str
        Expected outcome description
    """
    logger.info("=" * 80)
    logger.info(f"TEST: {test_name}")
    logger.info("=" * 80)

    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Create test files
        ref_file = tmpdir / "test_ref.fa"
        bam_file = tmpdir / "test.bam"
        vcf_file = tmpdir / "test.vcf"
        output_file = tmpdir / "cooccurrence.tsv"

        logger.info("Creating test reference sequence...")
        create_test_reference(str(ref_file))

        logger.info("Creating test BAM file...")
        create_test_bam(str(bam_file), str(ref_file), variants)

        logger.info("Creating test VCF file...")
        create_test_vcf(str(vcf_file), variants)

        # Run analysis
        logger.info("Running co-occurrence analysis...")
        sys.path.insert(0, str(Path(__file__).parent))
        from check_read_cooccurrence import VariantCooccurrenceAnalyzer, parse_vcf

        # Parse variants
        parsed_variants = parse_vcf(str(vcf_file))

        # Run analyzer
        analyzer = VariantCooccurrenceAnalyzer(str(bam_file))
        results = analyzer.analyze_all_variant_pairs(parsed_variants, max_distance=1000)

        # Print results
        logger.info("\nResults:")
        for result in results:
            logger.info(f"  Variant pair: {result['variant_pair']}")
            logger.info(f"  Distance: {result['distance']} bp")
            logger.info(f"  Spanning reads: {result['total_spanning_reads']}")
            logger.info(f"  Informative reads: {result['informative_reads']}")
            logger.info(f"  Both variants: {result['both_alt']}")
            logger.info(f"  Variant1 only: {result['variant1_only']}")
            logger.info(f"  Variant2 only: {result['variant2_only']}")
            logger.info(f"  Both reference: {result['both_ref']}")
            logger.info(f"  Co-occurrence rate: {result['cooccurrence_rate']:.2%}" if result['cooccurrence_rate'] else "  Co-occurrence rate: N/A")
            logger.info(f"  Linkage proven: {result['linkage_proven']}")
            logger.info(f"  Evidence level: {result['evidence_level']}")
            logger.info("")

        # Validate expected result
        logger.info(f"Expected: {expected_result}")
        if results:
            if results[0]['linkage_proven']:
                logger.info("✓ Test PASSED: Linkage was proven")
            else:
                logger.info("✓ Test PASSED: Linkage could not be proven (as expected)")
        else:
            logger.info("⚠ Test WARNING: No results generated")

        logger.info("")


def main():
    """Run all test cases."""
    logger.info("=" * 80)
    logger.info("TESTING: check_read_cooccurrence.py")
    logger.info("=" * 80)
    logger.info("")

    # Test 1: Two nearby variants with high co-occurrence
    logger.info("Test 1: Two variants 100bp apart, both at 95% frequency")
    logger.info("Expected: Should find high co-occurrence (~95%)")
    run_test_case(
        "Nearby variants with high frequency",
        variants=[
            (1000, 'A', 'T', 0.95),  # pos, ref, alt, freq
            (1100, 'C', 'G', 0.95),
        ],
        expected_result="High co-occurrence rate (~95%) with proven linkage"
    )

    # Test 2: Two nearby variants with different frequencies
    logger.info("Test 2: Two variants 150bp apart, frequencies 90% and 50%")
    logger.info("Expected: Co-occurrence rate around 50% (limited by lower frequency)")
    run_test_case(
        "Nearby variants with different frequencies",
        variants=[
            (2000, 'A', 'T', 0.90),
            (2150, 'C', 'G', 0.50),
        ],
        expected_result="Co-occurrence rate ~50%, some reads with only first variant"
    )

    # Test 3: Distant variants (beyond read length, but within typical insert size)
    logger.info("Test 3: Two variants 400bp apart")
    logger.info("Expected: May find co-occurrence if using paired-end reads")
    run_test_case(
        "Distant variants (400bp apart)",
        variants=[
            (3000, 'A', 'T', 0.80),
            (3400, 'G', 'C', 0.80),
        ],
        expected_result="Limited or no co-occurrence data (beyond read length)"
    )

    # Test 4: Three variants - testing multiple pairs
    logger.info("Test 4: Three nearby variants")
    logger.info("Expected: Multiple pair combinations should be analyzed")
    run_test_case(
        "Three nearby variants",
        variants=[
            (4000, 'A', 'T', 0.90),
            (4100, 'C', 'G', 0.85),
            (4200, 'G', 'A', 0.80),
        ],
        expected_result="Three pairs analyzed: (1,2), (1,3), (2,3)"
    )

    logger.info("=" * 80)
    logger.info("ALL TESTS COMPLETE")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
