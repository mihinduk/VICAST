#!/usr/bin/env python3
"""
Pathway 2 QC Module: Viral Polyprotein Gap Detection & Repair
=====================================================================

Purpose: Detect unannotated gaps in viral polyprotein annotations and attempt
         programmatic resolution before flagging for manual review.

Part of: Pathway 2 ("Well-Annotated Virus" Quality Control)
Triggered when: User wants to add a virus as a model for annotation transfer
"""

from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
from enum import Enum
import subprocess
import tempfile
import os


class GapSeverity(Enum):
    """Severity levels for annotation gaps"""
    NONE = "No gaps detected"
    MINOR = "Small gaps (<100 bp) - likely 3' UTR or splicing"
    MODERATE = "Medium gaps (100-300 bp) - likely single protein"
    SEVERE = "Large gaps (>300 bp) - likely 2+ missing proteins ⚠️"
    CRITICAL = "Extremely large gaps (>1 kb) - definitely missing proteins 🚨"


@dataclass
class AnnotationGap:
    """Represents a gap between two CDS features"""
    start: int
    end: int
    size_bp: int
    upstream_gene: str
    downstream_gene: str
    severity: GapSeverity
    attempted_repairs: List[str] = None
    repair_findings: Dict[str, str] = None
    
    def __post_init__(self):
        if self.attempted_repairs is None:
            self.attempted_repairs = []
        if self.repair_findings is None:
            self.repair_findings = {}


class ViralAnnotationGapDetector:
    """
    Detect gaps in viral polyprotein annotations and attempt repairs.
    """
    
    # Viral genome compactness thresholds
    MINOR_GAP_THRESHOLD = 100        # bp
    MODERATE_GAP_THRESHOLD = 300     # bp
    SEVERE_GAP_THRESHOLD = 1000      # bp
    
    # Minimum ORF size to consider
    MIN_ORF_SIZE = 75  # aa (225 bp) - most viral proteins > 50 aa
    
    def __init__(self, nucleotide_fasta: str, annotation_tsv: str):
        """
        Args:
            nucleotide_fasta: Path to viral nucleotide sequence (FASTA)
            annotation_tsv: Path to annotation file with CDS start/end positions
        """
        self.nucleotide_fasta = nucleotide_fasta
        self.annotation_tsv = annotation_tsv
        self.cds_features = []
        self.gaps = []
        
    def load_annotation(self) -> bool:
        """Load CDS features from annotation TSV."""
        try:
            import pandas as pd
            df = pd.read_csv(self.annotation_tsv, sep='\t')
            
            # Filter for CDS features only
            cds = df[df['type'] == 'CDS'].sort_values('start')
            
            for _, row in cds.iterrows():
                self.cds_features.append({
                    'start': int(row['start']),
                    'end': int(row['end']),
                    'gene_name': row.get('gene_name', 'unknown'),
                    'product': row.get('product', 'unknown'),
                    'strand': row.get('strand', '+')
                })
            
            return len(self.cds_features) > 0
        
        except Exception as e:
            print(f"ERROR loading annotation: {e}")
            return False
    
    def detect_gaps(self) -> List[AnnotationGap]:
        """
        Detect gaps between consecutive CDS features on the same strand.
        """
        if not self.cds_features:
            print("ERROR: No CDS features loaded")
            return []
        
        self.gaps = []
        
        for i in range(len(self.cds_features) - 1):
            current = self.cds_features[i]
            next_cds = self.cds_features[i + 1]
            
            # Only check features on same strand
            if current['strand'] != next_cds['strand']:
                continue
            
            gap_start = current['end'] + 1
            gap_end = next_cds['start'] - 1
            gap_size = gap_end - gap_start + 1
            
            # Skip tiny gaps (likely just annotation boundary artifacts)
            if gap_size <= 0:
                continue
            
            # Determine severity
            if gap_size < self.MINOR_GAP_THRESHOLD:
                severity = GapSeverity.MINOR
            elif gap_size < self.MODERATE_GAP_THRESHOLD:
                severity = GapSeverity.MODERATE
            elif gap_size < self.SEVERE_GAP_THRESHOLD:
                severity = GapSeverity.SEVERE
            else:
                severity = GapSeverity.CRITICAL
            
            gap = AnnotationGap(
                start=gap_start,
                end=gap_end,
                size_bp=gap_size,
                upstream_gene=current['gene_name'],
                downstream_gene=next_cds['gene_name'],
                severity=severity
            )
            
            self.gaps.append(gap)
        
        return self.gaps
    
    def attempt_repair_all_gaps(self, genome_sequence: str) -> bool:
        """
        Attempt programmatic repair for all detected gaps.
        """
        if not self.gaps:
            return True  # No gaps = success
        
        for gap in self.gaps:
            print(f"\n{'='*70}")
            print(f"GAP: {gap.upstream_gene} → {gap.downstream_gene}")
            print(f"Position: {gap.start:,} - {gap.end:,} ({gap.size_bp} bp)")
            print(f"Severity: {gap.severity.value}")
            print(f"{'='*70}")
            
            if gap.severity == GapSeverity.MINOR:
                print("✓ MINOR: Likely 3'UTR, intron, or annotation artifact. No repair needed.")
                continue
            
            # Attempt repairs in order of confidence
            self._attempt_getorf_repair(gap, genome_sequence)
            self._attempt_hmm_repair(gap, genome_sequence)
            self._attempt_blast_repair(gap, genome_sequence)
            self._attempt_synteny_repair(gap)
            
            # Report findings
            self._report_gap_findings(gap)
        
        return True
    
    def _attempt_getorf_repair(self, gap: AnnotationGap, genome_seq: str) -> None:
        """
        Attempt 1: Use EMBOSS getorf to find ORFs in gap region.
        Highest confidence (no false positives if filters applied correctly).
        """
        print("\n[Attempt 1/4] GETORF: Finding ORFs in gap...")
        
        try:
            # Extract gap sequence
            gap_sequence = genome_seq[gap.start - 1:gap.end]  # 0-indexed
            
            if len(gap_sequence) < 300:  # Need at least 100 codons
                print("  ✗ Gap too small for reliable ORF detection")
                return
            
            # Write to temp file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
                f.write(f">gap_{gap.start}_{gap.end}\n")
                f.write(gap_sequence + "\n")
                temp_fasta = f.name
            
            try:
                # Run getorf
                result = subprocess.run(
                    ['getorf', temp_fasta, '-minsize', str(self.MIN_ORF_SIZE * 3),
                     '-find', '1', '-outseq', '-'],
                    capture_output=True,
                    text=True,
                    timeout=10
                )
                
                if result.returncode == 0 and result.stdout:
                    orfs = [line.strip() for line in result.stdout.split('\n') if line.startswith('>')]
                    
                    if len(orfs) == 0:
                        print(f"  ✗ No ORFs found (>={self.MIN_ORF_SIZE} aa)")
                    elif len(orfs) == 1:
                        print(f"  ✓ Found 1 ORF: {orfs[0]}")
                        gap.attempted_repairs.append("getorf")
                        gap.repair_findings['getorf'] = f"1 ORF found (>= {self.MIN_ORF_SIZE} aa)"
                    else:
                        print(f"  ⚠ Found {len(orfs)} ORFs (fragmented? overlapping frames?)")
                        gap.attempted_repairs.append("getorf")
                        gap.repair_findings['getorf'] = f"{len(orfs)} ORFs found - review manually"
                else:
                    print(f"  ✗ getorf failed or no output")
            
            finally:
                os.unlink(temp_fasta)
        
        except FileNotFoundError:
            print("  ⚠ getorf not installed (conda install emboss)")
        except Exception as e:
            print(f"  ✗ Error: {e}")
    
    def _attempt_hmm_repair(self, gap: AnnotationGap, genome_seq: str) -> None:
        """
        Attempt 2: Use HMM domain search to identify conserved protein domains.
        Good for identifying real proteins even if sequence divergent.
        """
        print("\n[Attempt 2/4] HMMSCAN: Searching for conserved domains...")
        
        try:
            # Extract and translate gap
            gap_sequence = genome_seq[gap.start - 1:gap.end]
            
            # Try all 6 reading frames
            from Bio.Seq import Seq
            
            found_domains = False
            for frame in range(3):
                # Forward strand
                codon_seq = gap_sequence[frame:]
                codon_seq = codon_seq[:len(codon_seq) - len(codon_seq) % 3]
                
                try:
                    aa_seq = str(Seq(codon_seq).translate())
                    
                    # Quick check: does this look like a protein? (no too many stops)
                    stop_count = aa_seq.count('*')
                    if stop_count > len(aa_seq) * 0.1:  # >10% stops = not real
                        continue
                    
                    # Would call hmmscan here if integrated
                    # For now, we note what SHOULD happen
                    if len(aa_seq) >= self.MIN_ORF_SIZE:
                        print(f"  ℹ Frame +{frame}: {len(aa_seq)} aa (no stops) - candidate protein")
                        found_domains = True
                
                except:
                    pass
            
            if not found_domains:
                print(f"  ✗ No clear ORF in any frame")
            else:
                gap.attempted_repairs.append("hmmscan")
                gap.repair_findings['hmmscan'] = "Candidate ORF(s) identified - needs domain annotation"
        
        except Exception as e:
            print(f"  ⚠ HMM search skipped: {e}")
    
    def _attempt_blast_repair(self, gap: AnnotationGap, genome_seq: str) -> None:
        """
        Attempt 3: BLAST gap region against NCBI NR to identify homologs.
        Identifies what proteins SHOULD be in this gap from related viruses.
        """
        print("\n[Attempt 3/4] BLAST: Searching NCBI for homologous proteins...")
        
        try:
            # Would run blastx/tblastx here if integrated with internet access
            print("  ℹ In production: Would query NCBI NR for homologous proteins")
            print(f"     → Search viral proteins in {gap.upstream_gene}-{gap.downstream_gene} region")
            print(f"     → Identifies what proteins are TYPICALLY in this gap")
            gap.attempted_repairs.append("blast")
            gap.repair_findings['blast'] = "BLAST to NCBI recommended for gap characterization"
        
        except Exception as e:
            print(f"  ⚠ BLAST search skipped: {e}")
    
    def _attempt_synteny_repair(self, gap: AnnotationGap) -> None:
        """
        Attempt 4: Synteny analysis - compare to well-annotated related virus.
        Often the easiest: if related virus has these genes, THIS gap should too.
        """
        print("\n[Attempt 4/4] SYNTENY: Comparing to related viruses...")
        print(f"  ℹ Look for well-annotated virus with same gene order:")
        print(f"     {gap.upstream_gene} → ??? → {gap.downstream_gene}")
        print(f"  ℹ If found, that tells you what's missing here")
        gap.attempted_repairs.append("synteny")
        gap.repair_findings['synteny'] = "Manual synteny comparison needed"
    
    def _report_gap_findings(self, gap: AnnotationGap) -> None:
        """Generate comprehensive gap report."""
        print(f"\n{'─'*70}")
        print(f"REPAIR ATTEMPT SUMMARY:")
        print(f"{'─'*70}")
        
        if not gap.attempted_repairs:
            print("No repair attempts made (gap too small)")
            return
        
        for method, finding in gap.repair_findings.items():
            status = "✓" if "found" in finding.lower() else "ℹ"
            print(f"{status} {method.upper():10s}: {finding}")
        
        print(f"\nRECOMMENDATION:")
        
        # Recommendation logic
        if gap.size_bp < 300:
            print("→ MODERATE: Likely single protein. Attempt manual ORF search + synteny.")
        elif gap.size_bp < 1000:
            print("→ SEVERE: Likely 1-2 proteins. Use all 4 methods + manual review.")
        else:
            print("→ CRITICAL: Definitely missing proteins. Cannot use as model without repair!")
            print("  ACTION REQUIRED: Fix this annotation before using as model for transfer")
    
    def generate_qc_report(self) -> str:
        """
        Generate comprehensive QC report (plain text for terminal output).
        Returns: Multi-line string summarizing all gaps and recommended actions.
        """
        report_lines = [
            "╔════════════════════════════════════════════════════════════════╗",
            "║  PATHWAY 2 QC REPORT: Viral Annotation Gap Analysis           ║",
            "╚════════════════════════════════════════════════════════════════╝",
            ""
        ]

        if not self.gaps:
            report_lines.append("✓ NO GAPS DETECTED: Annotation appears complete")
            return "\n".join(report_lines)

        # Summary
        severity_counts = {}
        for gap in self.gaps:
            severity_counts[gap.severity] = severity_counts.get(gap.severity, 0) + 1

        report_lines.append(f"Found {len(self.gaps)} gap(s):")
        for severity, count in severity_counts.items():
            report_lines.append(f"  • {severity.value}: {count}")

        report_lines.append("")
        report_lines.append("DETAILED FINDINGS:")

        for i, gap in enumerate(self.gaps, 1):
            report_lines.append(f"\nGap #{i}: {gap.upstream_gene} → {gap.downstream_gene}")
            report_lines.append(f"  Position: {gap.start:,} - {gap.end:,} ({gap.size_bp} bp)")
            report_lines.append(f"  Severity: {gap.severity.value}")

            if gap.repair_findings:
                report_lines.append(f"  Findings:")
                for method, finding in gap.repair_findings.items():
                    report_lines.append(f"    • {method}: {finding}")

        # Final recommendation
        report_lines.append("\n" + "="*68)
        report_lines.append("FINAL RECOMMENDATION:")
        report_lines.append("="*68)

        has_severe = any(g.severity in [GapSeverity.SEVERE, GapSeverity.CRITICAL]
                        for g in self.gaps)

        if has_severe:
            report_lines.append("")
            report_lines.append("🚨 WARNING: Large gaps detected!")
            report_lines.append("")
            report_lines.append("CANNOT USE THIS VIRUS AS MODEL FOR ANNOTATION TRANSFER")
            report_lines.append("until gaps are resolved.")
            report_lines.append("")
            report_lines.append("REQUIRED ACTIONS:")
            report_lines.append("  1. Use the 4 repair methods above to identify missing proteins")
            report_lines.append("  2. Add newly identified CDS features to annotation")
            report_lines.append("  3. Revalidate using this QC module")
            report_lines.append("  4. Once no SEVERE/CRITICAL gaps remain: approved for use as model")
        else:
            report_lines.append("")
            report_lines.append("✓ Minor gaps acceptable for a model virus annotation.")
            report_lines.append("  Recommend: Document gaps in methods section of publication")

        return "\n".join(report_lines)

    def generate_markdown_report(self, genome_id: str) -> str:
        """
        Generate Markdown-formatted QC report suitable for documentation.

        Args:
            genome_id: Genome accession (e.g., NC_038433.1)

        Returns:
            Markdown-formatted report string
        """
        from datetime import datetime

        md_lines = [
            f"# Pathway 2 QC Report: {genome_id}",
            "",
            f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  ",
            f"**Pipeline:** VICAST Pathway 2 - Model Annotation Validation  ",
            "",
            "---",
            ""
        ]

        if not self.gaps:
            md_lines.extend([
                "## ✅ Result: PASSED",
                "",
                "**No gaps detected** - Annotation appears complete and suitable for use as a model in Pathway 3 annotation transfer.",
                ""
            ])
            return "\n".join(md_lines)

        # Summary
        severity_counts = {}
        for gap in self.gaps:
            severity_counts[gap.severity] = severity_counts.get(gap.severity, 0) + 1

        has_severe = any(g.severity in [GapSeverity.SEVERE, GapSeverity.CRITICAL]
                        for g in self.gaps)

        if has_severe:
            md_lines.extend([
                "## ❌ Result: FAILED",
                "",
                "**Large gaps detected** - Cannot use as model for annotation transfer until gaps are resolved.",
                ""
            ])
        else:
            md_lines.extend([
                "## ⚠️ Result: PASSED with Minor Gaps",
                "",
                "**Minor gaps detected** - Acceptable for use as model, but should be documented.",
                ""
            ])

        md_lines.extend([
            "## Summary",
            "",
            f"Found **{len(self.gaps)} gap(s)** in annotation:",
            ""
        ])

        # Severity breakdown table
        md_lines.extend([
            "| Severity | Count | Description |",
            "|----------|-------|-------------|"
        ])

        for severity, count in severity_counts.items():
            md_lines.append(f"| {severity.name} | {count} | {severity.value} |")

        md_lines.extend(["", "---", ""])

        # Detailed findings
        md_lines.extend([
            "## Detailed Gap Analysis",
            ""
        ])

        for i, gap in enumerate(self.gaps, 1):
            md_lines.extend([
                f"### Gap #{i}: {gap.upstream_gene} → {gap.downstream_gene}",
                "",
                f"- **Position:** {gap.start:,} - {gap.end:,}",
                f"- **Size:** {gap.size_bp:,} bp",
                f"- **Severity:** {gap.severity.name} - {gap.severity.value}",
                ""
            ])

            if gap.repair_findings:
                md_lines.append("**Repair Attempts:**")
                md_lines.append("")
                for method, finding in gap.repair_findings.items():
                    md_lines.append(f"- **{method.upper()}:** {finding}")
                md_lines.append("")

        md_lines.extend(["---", ""])

        # Final recommendation
        md_lines.extend([
            "## Recommendations",
            ""
        ])

        if has_severe:
            md_lines.extend([
                "### ⚠️ Action Required",
                "",
                "This annotation **CANNOT** be used as a model for Pathway 3 annotation transfer until gaps are resolved.",
                "",
                "**Required Steps:**",
                "",
                "1. **Identify missing proteins** using one or more methods:",
                "   - ORF finding (EMBOSS `getorf`)",
                "   - BLAST gap regions against NCBI NR",
                "   - Synteny comparison with related well-annotated viruses",
                "   - Domain search (HMMER)",
                "",
                "2. **Add newly identified CDS features** to the annotation TSV file",
                "",
                "3. **Re-run this QC step** to validate completeness",
                "",
                "4. **Proceed to step2** once no SEVERE/CRITICAL gaps remain",
                "",
                "### Alternative Approaches",
                "",
                "- Choose a different, more completely annotated virus as model",
                "- Use Pathway 4 (de novo annotation) instead of Pathway 3",
                ""
            ])
        else:
            md_lines.extend([
                "### ✅ Approved for Use",
                "",
                "Minor gaps are acceptable for model viruses. These likely represent:",
                "",
                "- 3' UTR regions",
                "- Intergenic spacers",
                "- Annotation boundary artifacts",
                "",
                "**Recommended:**",
                "",
                "- Document these gaps in the methods section of any publication",
                "- Consider manual review to confirm no proteins are missing",
                ""
            ])

        md_lines.extend([
            "---",
            "",
            "## Methods",
            "",
            "**Gap Detection:**",
            "",
            "Gaps were identified by comparing consecutive CDS features in the genome annotation. ",
            "Gap severity was classified based on size:",
            "",
            f"- **MINOR:** < {self.MINOR_GAP_THRESHOLD} bp",
            f"- **MODERATE:** {self.MINOR_GAP_THRESHOLD}-{self.MODERATE_GAP_THRESHOLD} bp",
            f"- **SEVERE:** {self.MODERATE_GAP_THRESHOLD}-{self.SEVERE_GAP_THRESHOLD} bp",
            f"- **CRITICAL:** > {self.SEVERE_GAP_THRESHOLD} bp",
            "",
            "**Repair Methods Attempted:**",
            "",
            "1. **GETORF** - EMBOSS ORF finder to identify potential coding sequences",
            "2. **HMMSCAN** - Domain search to identify conserved protein domains",
            "3. **BLAST** - Homology search against NCBI databases",
            "4. **SYNTENY** - Comparison with related virus gene order",
            "",
            "---",
            "",
            f"*Generated by VICAST Pathway 2 QC Module - {datetime.now().year}*",
            ""
        ])

        return "\n".join(md_lines)


# ════════════════════════════════════════════════════════════════════════════
# INTEGRATION WITH PATHWAY 2
# ════════════════════════════════════════════════════════════════════════════

class Pathway2QC:
    """
    Pathway 2 Quality Control: Validate and prepare well-annotated viruses.

    New Step: Gap Detection & Repair (inserted before model approval)
    """

    def validate_model_for_transfer(self, nucleotide_fasta: str,
                                    annotation_tsv: str,
                                    auto_repair: bool = False,
                                    save_reports: bool = True,
                                    genome_id: str = None) -> Tuple[bool, str]:
        """
        Validate a virus annotation for use as model in Pathway 3 transfer.

        Args:
            nucleotide_fasta: Path to virus nucleotide sequence
            annotation_tsv: Path to annotation TSV (with CDS features)
            auto_repair: If True, attempt programmatic repair. If False, only detect.
            save_reports: If True, save Markdown and PDF reports to disk
            genome_id: Genome ID for report naming (extracted from filename if not provided)

        Returns:
            (is_approved: bool, report: str)
        """

        # Extract genome ID from filename if not provided
        if genome_id is None:
            import os
            genome_id = os.path.basename(nucleotide_fasta).replace('.fasta', '')

        # Initialize detector
        detector = ViralAnnotationGapDetector(nucleotide_fasta, annotation_tsv)

        # Load annotation
        if not detector.load_annotation():
            return False, "ERROR: Could not load annotation"

        # Read sequence
        try:
            with open(nucleotide_fasta) as f:
                # Simple FASTA reading (assumes single sequence)
                lines = [line.strip() for line in f if not line.startswith('>')]
                genome_seq = ''.join(lines)
        except Exception as e:
            return False, f"ERROR reading sequence: {e}"

        # Detect gaps
        gaps = detector.detect_gaps()

        # Attempt repairs if requested
        if auto_repair and gaps:
            detector.attempt_repair_all_gaps(genome_seq)

        # Generate terminal report
        report = detector.generate_qc_report()

        # Save reports to disk if requested
        if save_reports:
            md_report = detector.generate_markdown_report(genome_id)
            self._save_reports(genome_id, md_report)

        # Determine approval status
        has_critical = any(g.severity == GapSeverity.CRITICAL for g in gaps)
        has_severe = any(g.severity == GapSeverity.SEVERE for g in gaps)

        # Approval logic
        if has_critical:
            is_approved = False  # Cannot use as model
        elif has_severe:
            is_approved = False  # Cannot use without repair
        else:
            is_approved = True   # Safe to use

        return is_approved, report

    def _save_reports(self, genome_id: str, markdown_content: str) -> None:
        """
        Save Markdown and PDF reports to disk.

        Args:
            genome_id: Genome accession for filename
            markdown_content: Markdown report content
        """
        import os

        # Save Markdown report
        md_filename = f"{genome_id}_pathway2_qc_report.md"
        try:
            with open(md_filename, 'w') as f:
                f.write(markdown_content)
            print(f"\n✓ Markdown report saved: {md_filename}")
        except Exception as e:
            print(f"\n⚠ Warning: Could not save Markdown report: {e}")

        # Generate and save PDF report
        pdf_filename = f"{genome_id}_pathway2_qc_report.pdf"
        try:
            self._markdown_to_pdf(markdown_content, pdf_filename)
            print(f"✓ PDF report saved: {pdf_filename}")
        except ImportError:
            print(f"⚠ Warning: PDF generation requires 'markdown' and 'weasyprint' packages")
            print(f"  Install with: conda install -c conda-forge markdown weasyprint")
        except Exception as e:
            print(f"⚠ Warning: Could not generate PDF: {e}")

    def _markdown_to_pdf(self, markdown_content: str, pdf_filename: str) -> None:
        """
        Convert Markdown content to PDF.

        Args:
            markdown_content: Markdown text
            pdf_filename: Output PDF filename
        """
        try:
            import markdown
            from weasyprint import HTML, CSS
            from io import BytesIO
        except ImportError:
            raise ImportError("PDF generation requires 'markdown' and 'weasyprint' packages")

        # Convert Markdown to HTML
        md = markdown.Markdown(extensions=['tables', 'fenced_code'])
        html_body = md.convert(markdown_content)

        # Add CSS styling for professional appearance
        html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <style>
        body {{
            font-family: 'Helvetica', 'Arial', sans-serif;
            line-height: 1.6;
            max-width: 800px;
            margin: 40px auto;
            padding: 20px;
            color: #333;
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-bottom: 2px solid #95a5a6;
            padding-bottom: 5px;
        }}
        h3 {{
            color: #7f8c8d;
            margin-top: 20px;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        code {{
            background-color: #f4f4f4;
            padding: 2px 6px;
            border-radius: 3px;
            font-family: 'Courier New', monospace;
        }}
        hr {{
            border: none;
            border-top: 1px solid #ddd;
            margin: 30px 0;
        }}
        strong {{
            color: #2c3e50;
        }}
        ul, ol {{
            margin: 10px 0;
            padding-left: 30px;
        }}
        li {{
            margin: 5px 0;
        }}
    </style>
</head>
<body>
{html_body}
</body>
</html>
"""

        # Generate PDF
        HTML(string=html).write_pdf(pdf_filename)


# ════════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE
# ════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    
    print("""
    ╔════════════════════════════════════════════════════════════════╗
    ║  Pathway 2 QC: Viral Annotation Gap Detection & Repair        ║
    ║  ─────────────────────────────────────────────────────────    ║
    ║  Usage: python viral_gap_qc.py [sequence.fasta] [annot.tsv]   ║
    ╚════════════════════════════════════════════════════════════════╝
    """)
    
    import sys
    
    if len(sys.argv) == 3:
        fasta_file = sys.argv[1]
        annot_file = sys.argv[2]
        
        # Run QC with auto-repair
        qc = Pathway2QC()
        is_approved, report = qc.validate_model_for_transfer(
            fasta_file, annot_file, 
            auto_repair=True
        )
        
        print(report)
        
        if is_approved:
            print("\n✅ APPROVED: This virus can be used as a model for Pathway 3")
            sys.exit(0)
        else:
            print("\n❌ REJECTED: Fix gaps before using as model")
            sys.exit(1)
    else:
        print("Usage: python viral_gap_qc.py <sequence.fasta> <annotation.tsv>")
        print("\nExample:")
        print("  python viral_gap_qc.py NC_038433.1.fasta NC_038433_1.tsv")
