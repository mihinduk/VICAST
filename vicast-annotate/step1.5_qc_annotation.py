#!/usr/bin/env python3
"""
STEP 1.5 (Pathway 2 QC): Quality Control for Model Virus Annotations

This script validates that a virus annotation is suitable for use as a model
in Pathway 3 annotation transfer. It detects gaps in polyprotein annotations
that would result in missing proteins during transfer.

MANUAL REVIEW WORKFLOW:
This is a manual review step. The script generates reports and provides
recommendations, but YOU decide whether to proceed or fix issues first.

Usage:
    python step1.5_qc_annotation.py <genome_id> [options]

Example:
    python step1.5_qc_annotation.py NC_038433.1

Input Files (auto-detected):
    - <genome_id>.fasta (nucleotide sequence)
    - <genome_id>_no_polyprotein.tsv (annotation from step1)

Output:
    - QC report printed to stdout
    - <genome_id>_pathway2_qc_report.md (Markdown report)
    - <genome_id>_pathway2_qc_report.pdf (PDF report, if dependencies available)
    - Exit code 0 if approved or --force used, 1 if gaps found

Integration:
    Run this between step1 and step2 to validate model quality:
    1. python step1_parse_viral_genome.py NC_038433.1
    2. python step1.5_qc_annotation.py NC_038433.1  ← NEW (manual review)
    3. Review generated reports, fix issues if needed
    4. python step2_add_to_snpeff.py NC_038433.1 NC_038433.1_no_polyprotein.tsv
"""

import sys
import os
import argparse
from pathlib import Path

# Import the QC module
from viral_gap_qc import Pathway2QC, GapSeverity


def find_input_files(genome_id):
    """
    Find the FASTA and TSV files for a genome ID.

    Args:
        genome_id: Genome accession (e.g., NC_038433.1)

    Returns:
        tuple: (fasta_file, tsv_file) or (None, None) if not found
    """
    # Look for files in current directory
    fasta_file = f"{genome_id}.fasta"
    tsv_file = f"{genome_id}_no_polyprotein.tsv"

    if not os.path.exists(fasta_file):
        print(f"ERROR: FASTA file not found: {fasta_file}")
        print(f"Expected output from step1_parse_viral_genome.py")
        return None, None

    if not os.path.exists(tsv_file):
        print(f"ERROR: TSV file not found: {tsv_file}")
        print(f"Expected output from step1_parse_viral_genome.py")
        return None, None

    return fasta_file, tsv_file


def main():
    parser = argparse.ArgumentParser(
        description='STEP 1.5: QC model virus annotation before use in Pathway 3',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This quality control step validates that a virus annotation is complete enough
to serve as a model for annotation transfer in Pathway 3.

WHY THIS MATTERS:
Pathway 3 transfers annotations from a model virus to a poorly-annotated virus
using BLASTx homology. If the model has gaps (unannotated regions containing
proteins), those proteins CANNOT be transferred, resulting in incomplete
annotations for the target virus.

WHAT IT CHECKS:
- Detects gaps between consecutive CDS features
- Classifies severity (MINOR → CRITICAL)
- Attempts to identify what proteins are missing
- Generates Markdown and PDF reports for documentation

MANUAL REVIEW WORKFLOW:
This is NOT an automatic gate. You review the generated reports and decide
whether to fix issues or proceed. Use --force if you want to proceed despite gaps.

EXAMPLE WORKFLOW:

Step 1: Parse genome
  python step1_parse_viral_genome.py NC_038433.1
  → Creates NC_038433.1.fasta, NC_038433.1_no_polyprotein.tsv

Step 1.5: QC annotation (MANUAL REVIEW)
  python step1.5_qc_annotation.py NC_038433.1
  → Generates reports: NC_038433.1_pathway2_qc_report.{md,pdf}
  → You review the reports
  → If gaps found: you decide whether to fix or proceed

  If fixing:
    - Use getorf, BLAST, or synteny to identify missing proteins
    - Add to TSV manually
    - Re-run step1.5 to validate

  If proceeding anyway:
    - Re-run with --force flag

Step 2: Add to SnpEff
  python step2_add_to_snpeff.py NC_038433.1 NC_038433.1_no_polyprotein.tsv

INTEGRATION OPTIONS:

A. Manual review (recommended):
  python step1.5_qc_annotation.py $GENOME_ID
  # Review reports, make decision

B. Automated gate (strict):
  if ! python step1.5_qc_annotation.py $GENOME_ID; then
      echo "QC failed - review reports and fix"
      exit 1
  fi

C. Warning only (soft fail):
  python step1.5_qc_annotation.py $GENOME_ID || echo "WARNING: Review QC reports"
        """
    )

    parser.add_argument('genome_id',
                       help='Genome ID (e.g., NC_038433.1)')

    parser.add_argument('--auto-repair', action='store_true',
                       help='Attempt automatic repair using ORF finding (experimental)')

    parser.add_argument('--force', action='store_true',
                       help='Proceed even if gaps found (for manual review workflow)')

    parser.add_argument('--no-reports', action='store_true',
                       help='Skip saving Markdown/PDF reports (only print to stdout)')

    args = parser.parse_args()

    print("="*70)
    print("STEP 1.5: Pathway 2 QC - Model Annotation Validation")
    print("="*70)
    print(f"\nGenome ID: {args.genome_id}")
    print(f"Workflow: Manual Review (generate reports, you decide)")

    # Find input files
    fasta_file, tsv_file = find_input_files(args.genome_id)

    if not fasta_file or not tsv_file:
        print("\nERROR: Input files not found")
        print("\nExpected workflow:")
        print(f"  1. python step1_parse_viral_genome.py {args.genome_id}")
        print(f"  2. python step1.5_qc_annotation.py {args.genome_id}  ← you are here")
        sys.exit(1)

    print(f"  FASTA: {fasta_file}")
    print(f"  TSV:   {tsv_file}")

    # Run QC
    print("\nRunning quality control checks...")
    print("="*70)

    qc = Pathway2QC()

    try:
        is_approved, report = qc.validate_model_for_transfer(
            fasta_file,
            tsv_file,
            auto_repair=args.auto_repair,
            save_reports=(not args.no_reports),
            genome_id=args.genome_id
        )
    except Exception as e:
        print(f"\nERROR: QC failed with exception:")
        print(f"  {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    # Print report
    print("\n" + report)
    print("\n" + "="*70)

    # Determine outcome
    if is_approved:
        print("✓ QC PASSED: Annotation approved for use as Pathway 3 model")
        print("="*70)
        print("\nNext step:")
        print(f"  python step2_add_to_snpeff.py {args.genome_id} {tsv_file}")
        sys.exit(0)

    elif args.force:
        print("⚠ QC IDENTIFIED GAPS but --force specified: Proceeding")
        print("="*70)
        print("\nMANUAL REVIEW DECISION: Proceed despite gaps")
        print("  - You've reviewed the reports and decided to proceed")
        print("  - Document these gaps in your methods section")
        print("  - Be aware: Missing proteins cannot be transferred in Pathway 3")
        print("\nNext step:")
        print(f"  python step2_add_to_snpeff.py {args.genome_id} {tsv_file}")
        print("="*70)
        sys.exit(0)

    else:
        print("✗ QC IDENTIFIED GAPS: Review required")
        print("="*70)
        print("\nMANUAL REVIEW OPTIONS:")
        print("\n1. FIX THE GAPS (recommended for Pathway 3 models):")
        print("   a. Review the generated reports")
        print("   b. Identify missing proteins using suggested methods:")
        print("      - ORF finding (getorf)")
        print("      - BLAST gap region to NCBI")
        print("      - Compare to related well-annotated viruses (synteny)")
        print("   c. Add missing CDS features to TSV file")
        print("   d. Re-run this QC step to validate")
        print("\n2. PROCEED ANYWAY (use with caution):")
        print(f"   python step1.5_qc_annotation.py {args.genome_id} --force")
        print("   → Document gaps in methods section")
        print("   → Missing proteins cannot be transferred to target viruses")
        print("\n3. CHOOSE DIFFERENT MODEL:")
        print("   → Find a more completely annotated virus as model")
        print("   → Or use Pathway 4 (de novo annotation) instead")
        print("="*70)
        sys.exit(1)


if __name__ == '__main__':
    main()
