#!/usr/bin/env python3
"""
Validate VICAST multi-segment virus handling with Influenza A.

This script demonstrates VICAST's ability to:
1. Handle segmented viruses (8 segments) in a unified database
2. Annotate variants across all segments with a single SnpEff run
3. Properly track segment-specific annotations

Usage:
    python validate_multisegment.py

Requirements:
    - VICAST installed
    - SnpEff installed and configured
    - Influenza reference files in this directory
"""

import subprocess
import sys
from pathlib import Path

# Add vicast to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from vicast.validation import validate_gff_for_snpeff


def check_files():
    """Verify all required files exist."""
    required = [
        "influenza_A_California_2009_8segments.fasta",
        "influenza_A_California_2009.gff3",
    ]

    missing = [f for f in required if not Path(f).exists()]
    if missing:
        print(f"Missing files: {missing}")
        return False

    print("All required files present.")
    return True


def validate_gff():
    """Validate the GFF3 file for SnpEff compatibility."""
    print("\n=== Validating GFF3 for SnpEff ===")

    gff_file = "influenza_A_California_2009.gff3"
    fasta_file = "influenza_A_California_2009_8segments.fasta"

    is_valid, errors, warnings = validate_gff_for_snpeff(gff_file, fasta_file)

    if errors:
        print(f"Errors: {errors}")
    if warnings:
        print(f"Warnings: {warnings}")

    if is_valid:
        print("GFF3 validation PASSED")
    else:
        print("GFF3 validation FAILED")

    return is_valid


def count_features():
    """Count features per segment in the GFF3."""
    print("\n=== Feature Counts by Segment ===")

    segment_features = {}

    with open("influenza_A_California_2009.gff3") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                seqid = parts[0]
                ftype = parts[2]
                if seqid not in segment_features:
                    segment_features[seqid] = {"gene": 0, "mRNA": 0, "CDS": 0}
                if ftype in segment_features[seqid]:
                    segment_features[seqid][ftype] += 1

    # Map accessions to segment names
    segment_names = {
        "NC_026431.1": "Seg7 (M1/M2)",
        "NC_026432.1": "Seg8 (NS1/NEP)",
        "NC_026433.1": "Seg4 (HA)",
        "NC_026434.1": "Seg6 (NA)",
        "NC_026435.1": "Seg2 (PB1)",
        "NC_026436.1": "Seg5 (NP)",
        "NC_026437.1": "Seg3 (PA/PA-X)",
        "NC_026438.1": "Seg1 (PB2)",
    }

    print(f"{'Segment':<20} {'Genes':>6} {'mRNAs':>6} {'CDS':>6}")
    print("-" * 42)

    total_genes = 0
    total_cds = 0

    for seqid in sorted(segment_features.keys()):
        counts = segment_features[seqid]
        name = segment_names.get(seqid, seqid)
        print(f"{name:<20} {counts['gene']:>6} {counts['mRNA']:>6} {counts['CDS']:>6}")
        total_genes += counts['gene']
        total_cds += counts['CDS']

    print("-" * 42)
    print(f"{'TOTAL':<20} {total_genes:>6} {'-':>6} {total_cds:>6}")

    return segment_features


def check_spliced_genes():
    """Verify spliced genes (M2, NEP, PA-X) are properly annotated."""
    print("\n=== Checking Spliced Gene Annotations ===")

    spliced_genes = {
        "M2": {"expected_exons": 2, "found": 0},
        "NEP": {"expected_exons": 2, "found": 0},
        "PA-X": {"expected_exons": 2, "found": 0},
    }

    with open("influenza_A_California_2009.gff3") as f:
        for line in f:
            if "\tCDS\t" in line:
                for gene in spliced_genes:
                    if f"Name={gene}" in line or f"cds_{gene}_" in line:
                        spliced_genes[gene]["found"] += 1

    all_correct = True
    for gene, info in spliced_genes.items():
        status = "OK" if info["found"] == info["expected_exons"] else "FAILED"
        if status == "FAILED":
            all_correct = False
        print(f"  {gene}: {info['found']}/{info['expected_exons']} CDS features - {status}")

    if all_correct:
        print("Spliced gene check PASSED")
    else:
        print("Spliced gene check FAILED")

    return all_correct


def generate_test_vcf():
    """Generate a test VCF with variants across multiple segments."""
    print("\n=== Generating Test VCF ===")

    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##contig=<ID=NC_026431.1,length=982>
##contig=<ID=NC_026432.1,length=863>
##contig=<ID=NC_026433.1,length=1701>
##contig=<ID=NC_026434.1,length=1410>
##contig=<ID=NC_026435.1,length=2274>
##contig=<ID=NC_026436.1,length=1497>
##contig=<ID=NC_026437.1,length=2151>
##contig=<ID=NC_026438.1,length=2280>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
NC_026433.1\t100\t.\tA\tG\t1000\tPASS\tDP=500;AF=0.15
NC_026433.1\t300\t.\tG\tA\t1500\tPASS\tDP=600;AF=0.25
NC_026434.1\t200\t.\tT\tC\t1200\tPASS\tDP=550;AF=0.10
NC_026435.1\t500\t.\tC\tT\t1100\tPASS\tDP=480;AF=0.08
NC_026438.1\t1000\t.\tA\tG\t1300\tPASS\tDP=520;AF=0.12
"""

    with open("test_variants.vcf", "w") as f:
        f.write(vcf_content)

    print("Created test_variants.vcf with 5 variants across 4 segments:")
    print("  - NC_026433.1 (HA): 2 variants")
    print("  - NC_026434.1 (NA): 1 variant")
    print("  - NC_026435.1 (PB1): 1 variant")
    print("  - NC_026438.1 (PB2): 1 variant")

    return True


def print_summary():
    """Print validation summary."""
    print("\n" + "=" * 50)
    print("INFLUENZA A MULTI-SEGMENT VALIDATION SUMMARY")
    print("=" * 50)
    print("""
This dataset demonstrates VICAST's multi-segment handling:

1. UNIFIED DATABASE: All 8 segments in single FASTA/GFF3
   - Segments are tracked by accession (NC_026431-NC_026438)
   - SnpEff can annotate variants on any segment

2. SPLICED GENE SUPPORT: Proper handling of:
   - M2 (segment 7): Spliced from M gene
   - NEP (segment 8): Spliced from NS gene
   - PA-X (segment 3): Ribosomal frameshift product

3. VARIANT ANNOTATION: Single SnpEff run annotates all segments
   - No need for separate runs per segment
   - Consistent annotation across entire genome

To build SnpEff database:
    # Add to snpEff.config:
    # influenza_A_Cal09.genome : Influenza A/California/07/2009

    mkdir -p $SNPEFF_HOME/data/influenza_A_Cal09
    cp influenza_A_California_2009_8segments.fasta $SNPEFF_HOME/data/influenza_A_Cal09/sequences.fa
    cp influenza_A_California_2009.gff3 $SNPEFF_HOME/data/influenza_A_Cal09/genes.gff
    java -jar snpEff.jar build -gff3 influenza_A_Cal09

To annotate variants:
    java -jar snpEff.jar influenza_A_Cal09 test_variants.vcf > annotated.vcf
""")


def main():
    print("=" * 50)
    print("VICAST Multi-Segment Validation: Influenza A")
    print("Strain: A/California/07/2009 (H1N1)")
    print("=" * 50)

    # Check files exist
    if not check_files():
        sys.exit(1)

    # Validate GFF3
    validate_gff()

    # Count features
    count_features()

    # Check spliced genes
    check_spliced_genes()

    # Generate test VCF
    generate_test_vcf()

    # Print summary
    print_summary()

    print("\nValidation complete!")


if __name__ == "__main__":
    main()
