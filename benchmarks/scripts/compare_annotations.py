#!/usr/bin/env python3
"""
Compare annotation outputs from different tools.

Calculates metrics for comparing GFF3 annotation files:
- Gene overlap (Jaccard similarity)
- Boundary accuracy
- Feature counts
- Product name matching

Usage:
    python compare_annotations.py --vicast vicast.gff3 --reference ncbi.gff3
    python compare_annotations.py --vicast vicast.gff3 --vadr vadr.gff3 --output comparison.json
"""

import argparse
import json
import re
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


@dataclass
class Feature:
    """Represents a GFF3 feature."""
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    strand: str
    attributes: Dict[str, str]

    @property
    def gene_name(self) -> Optional[str]:
        return self.attributes.get("gene") or self.attributes.get("Name")

    @property
    def product(self) -> Optional[str]:
        return self.attributes.get("product")

    @property
    def feature_id(self) -> Optional[str]:
        return self.attributes.get("ID")

    def overlaps(self, other: "Feature", min_overlap: float = 0.5) -> bool:
        """Check if this feature overlaps with another."""
        if self.seqid != other.seqid or self.strand != other.strand:
            return False

        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)

        if overlap_start >= overlap_end:
            return False

        overlap_len = overlap_end - overlap_start
        self_len = self.end - self.start
        other_len = other.end - other.start

        # Check if overlap is significant for both features
        return (overlap_len / self_len >= min_overlap or
                overlap_len / other_len >= min_overlap)

    def boundary_match(self, other: "Feature", tolerance: int = 3) -> Tuple[bool, bool]:
        """Check if boundaries match within tolerance."""
        start_match = abs(self.start - other.start) <= tolerance
        end_match = abs(self.end - other.end) <= tolerance
        return start_match, end_match


def parse_gff3(filepath: Path) -> List[Feature]:
    """Parse a GFF3 file and return list of features."""
    features = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            # Parse attributes
            attributes = {}
            for attr in parts[8].split(";"):
                if "=" in attr:
                    key, value = attr.split("=", 1)
                    attributes[key] = value

            feature = Feature(
                seqid=parts[0],
                source=parts[1],
                feature_type=parts[2],
                start=int(parts[3]),
                end=int(parts[4]),
                strand=parts[6],
                attributes=attributes
            )
            features.append(feature)

    return features


def parse_vadr_tbl(filepath: Path) -> List[Feature]:
    """Parse VADR feature table format."""
    features = []
    current_feature = None
    current_coords = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip()

            # New feature definition
            if line.startswith(">Feature"):
                seqid = line.split()[1] if len(line.split()) > 1 else "unknown"
                continue

            # Coordinate line (starts with numbers)
            coord_match = re.match(r"^(\d+)\s+(\d+)\s+(\w+)", line)
            if coord_match:
                if current_feature:
                    features.append(current_feature)

                start, end, ftype = coord_match.groups()
                strand = "+" if int(start) < int(end) else "-"
                if strand == "-":
                    start, end = end, start

                current_feature = Feature(
                    seqid=seqid,
                    source="vadr",
                    feature_type=ftype,
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    attributes={}
                )
                continue

            # Attribute line (starts with tabs)
            attr_match = re.match(r"^\t+(\w+)\s+(.+)", line)
            if attr_match and current_feature:
                key, value = attr_match.groups()
                current_feature.attributes[key] = value

    if current_feature:
        features.append(current_feature)

    return features


@dataclass
class ComparisonMetrics:
    """Metrics from comparing two annotation sets."""
    file1: str
    file2: str
    file1_label: str
    file2_label: str

    # Feature counts
    file1_total_features: int
    file2_total_features: int
    file1_cds_count: int
    file2_cds_count: int
    file1_gene_count: int
    file2_gene_count: int

    # Overlap metrics
    shared_cds: int
    file1_only_cds: int
    file2_only_cds: int
    cds_jaccard_similarity: float

    # Boundary accuracy
    exact_boundary_matches: int
    start_only_matches: int
    end_only_matches: int
    no_boundary_matches: int

    # Name matching
    gene_name_matches: int
    product_matches: int


def compare_annotations(
    file1_path: Path,
    file2_path: Path,
    file1_label: str = "file1",
    file2_label: str = "file2"
) -> ComparisonMetrics:
    """Compare two annotation files."""

    # Parse files
    if file1_path.suffix == ".tbl":
        features1 = parse_vadr_tbl(file1_path)
    else:
        features1 = parse_gff3(file1_path)

    if file2_path.suffix == ".tbl":
        features2 = parse_vadr_tbl(file2_path)
    else:
        features2 = parse_gff3(file2_path)

    # Filter to CDS features
    cds1 = [f for f in features1 if f.feature_type == "CDS"]
    cds2 = [f for f in features2 if f.feature_type == "CDS"]

    genes1 = [f for f in features1 if f.feature_type == "gene"]
    genes2 = [f for f in features2 if f.feature_type == "gene"]

    # Find overlapping CDS
    matched1 = set()
    matched2 = set()
    exact_boundary = 0
    start_only = 0
    end_only = 0
    no_boundary = 0
    gene_name_match = 0
    product_match = 0

    for i, f1 in enumerate(cds1):
        for j, f2 in enumerate(cds2):
            if j in matched2:
                continue
            if f1.overlaps(f2):
                matched1.add(i)
                matched2.add(j)

                # Check boundaries
                start_match, end_match = f1.boundary_match(f2)
                if start_match and end_match:
                    exact_boundary += 1
                elif start_match:
                    start_only += 1
                elif end_match:
                    end_only += 1
                else:
                    no_boundary += 1

                # Check names
                if f1.gene_name and f2.gene_name and f1.gene_name.lower() == f2.gene_name.lower():
                    gene_name_match += 1
                if f1.product and f2.product and f1.product.lower() == f2.product.lower():
                    product_match += 1

                break

    shared = len(matched1)
    file1_only = len(cds1) - shared
    file2_only = len(cds2) - shared

    # Jaccard similarity
    union = len(cds1) + len(cds2) - shared
    jaccard = shared / union if union > 0 else 0

    return ComparisonMetrics(
        file1=str(file1_path),
        file2=str(file2_path),
        file1_label=file1_label,
        file2_label=file2_label,
        file1_total_features=len(features1),
        file2_total_features=len(features2),
        file1_cds_count=len(cds1),
        file2_cds_count=len(cds2),
        file1_gene_count=len(genes1),
        file2_gene_count=len(genes2),
        shared_cds=shared,
        file1_only_cds=file1_only,
        file2_only_cds=file2_only,
        cds_jaccard_similarity=jaccard,
        exact_boundary_matches=exact_boundary,
        start_only_matches=start_only,
        end_only_matches=end_only,
        no_boundary_matches=no_boundary,
        gene_name_matches=gene_name_match,
        product_matches=product_match
    )


def print_comparison(metrics: ComparisonMetrics):
    """Print comparison results in a readable format."""
    print("\n" + "=" * 60)
    print("ANNOTATION COMPARISON")
    print("=" * 60)

    print(f"\n{metrics.file1_label}: {metrics.file1}")
    print(f"{metrics.file2_label}: {metrics.file2}")

    print("\n--- Feature Counts ---")
    print(f"{'Metric':<30} {metrics.file1_label:>12} {metrics.file2_label:>12}")
    print("-" * 56)
    print(f"{'Total features':<30} {metrics.file1_total_features:>12} {metrics.file2_total_features:>12}")
    print(f"{'CDS features':<30} {metrics.file1_cds_count:>12} {metrics.file2_cds_count:>12}")
    print(f"{'Gene features':<30} {metrics.file1_gene_count:>12} {metrics.file2_gene_count:>12}")

    print("\n--- CDS Overlap ---")
    print(f"Shared CDS: {metrics.shared_cds}")
    print(f"{metrics.file1_label} only: {metrics.file1_only_cds}")
    print(f"{metrics.file2_label} only: {metrics.file2_only_cds}")
    print(f"Jaccard similarity: {metrics.cds_jaccard_similarity:.3f}")

    print("\n--- Boundary Accuracy (for shared CDS) ---")
    total_shared = metrics.shared_cds
    if total_shared > 0:
        print(f"Exact matches (both boundaries): {metrics.exact_boundary_matches} ({100*metrics.exact_boundary_matches/total_shared:.1f}%)")
        print(f"Start only matches: {metrics.start_only_matches} ({100*metrics.start_only_matches/total_shared:.1f}%)")
        print(f"End only matches: {metrics.end_only_matches} ({100*metrics.end_only_matches/total_shared:.1f}%)")
        print(f"No boundary matches: {metrics.no_boundary_matches} ({100*metrics.no_boundary_matches/total_shared:.1f}%)")

    print("\n--- Name Matching (for shared CDS) ---")
    if total_shared > 0:
        print(f"Gene name matches: {metrics.gene_name_matches} ({100*metrics.gene_name_matches/total_shared:.1f}%)")
        print(f"Product matches: {metrics.product_matches} ({100*metrics.product_matches/total_shared:.1f}%)")

    print("\n" + "=" * 60)


def main():
    parser = argparse.ArgumentParser(description="Compare annotation files")
    parser.add_argument("--vicast", help="VICAST GFF3 file")
    parser.add_argument("--vadr", help="VADR output file (.tbl or .gff3)")
    parser.add_argument("--reference", help="Reference annotation file")
    parser.add_argument("--file1", help="First annotation file")
    parser.add_argument("--file2", help="Second annotation file")
    parser.add_argument("--output", help="Output JSON file for metrics")

    args = parser.parse_args()

    # Determine which files to compare
    if args.vicast and args.vadr:
        file1, label1 = Path(args.vicast), "VICAST"
        file2, label2 = Path(args.vadr), "VADR"
    elif args.vicast and args.reference:
        file1, label1 = Path(args.vicast), "VICAST"
        file2, label2 = Path(args.reference), "Reference"
    elif args.file1 and args.file2:
        file1, label1 = Path(args.file1), "File1"
        file2, label2 = Path(args.file2), "File2"
    else:
        parser.error("Must specify --vicast with --vadr or --reference, or --file1 with --file2")

    # Run comparison
    metrics = compare_annotations(file1, file2, label1, label2)

    # Print results
    print_comparison(metrics)

    # Save to JSON if requested
    if args.output:
        with open(args.output, "w") as f:
            json.dump(asdict(metrics), f, indent=2)
        print(f"\nMetrics saved to: {args.output}")


if __name__ == "__main__":
    main()
