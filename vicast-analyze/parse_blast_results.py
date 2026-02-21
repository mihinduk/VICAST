#!/usr/bin/env python3
"""
Parse BLAST results for contamination screening and quality assessment.

For each assembled contig, analyzes BLAST hits to:
1. Confirm the reference virus is present
2. Detect viral contaminants (other viruses / co-infections)
3. Detect non-viral contaminants (bacteria, mycoplasma, fungi, human)
4. Apply query coverage thresholds to filter spurious hits

Outputs:
  - top_hits.tsv:  ALL hits per contig (up to max_target_seqs), with tied-best flag
  - viral_blast.tsv: Viral hits with query coverage
  - blast_report_section.txt: Detailed contamination screening report

Usage:
    python parse_blast_results.py blast_all.tsv contigs.fa --accession NC_001474.2
    python parse_blast_results.py blast_all.tsv contigs.fa --accession NC_001474.2 \
        --top-hits out_top_hits.tsv --viral-hits out_viral.tsv --report out_report.txt \
        --min-coverage 10
"""

import argparse
import json
import os
import re
import sys
from collections import defaultdict, OrderedDict
from pathlib import Path


def resolve_accessions(accession):
    """Resolve a single accession to a set of accessions to match.

    For segmented viruses (e.g., influenza_pr8), the custom accession won't
    appear in BLAST results — the individual segment accessions will. This
    function looks up segment_accessions from the prebuilt manifest.

    Returns a set of accession strings to check against BLAST subject IDs.
    """
    accessions = {accession}
    if not accession:
        return accessions

    # Search for manifest.json in standard locations
    script_dir = Path(__file__).parent
    manifest_paths = [
        script_dir.parent / "prebuilt_databases" / "manifest.json",
        script_dir / "manifest.json",
        Path("manifest.json"),
    ]

    for mpath in manifest_paths:
        if mpath.exists():
            try:
                with open(mpath, 'r') as f:
                    manifest = json.load(f)
                for entry in manifest.get("genomes", []):
                    if entry.get("accession") == accession:
                        segments = entry.get("segment_accessions", [])
                        if segments:
                            accessions.update(segments)
                            print(f"Resolved {accession} → {len(segments)} segment accessions",
                                  file=sys.stderr)
                        break
            except (json.JSONDecodeError, KeyError):
                pass
            break  # Only try first found manifest

    return accessions


def accession_matches(subject_id, accession_set):
    """Check if a BLAST subject_id matches any accession in the set."""
    return any(acc in subject_id for acc in accession_set)


def parse_contig_lengths(fasta_path):
    """Parse contig lengths from FASTA file.

    Handles MEGAHIT headers (len=XXXX) and falls back to counting
    actual sequence length for other FASTA formats.
    """
    lengths = {}
    current_id = None
    current_len = 0

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                # Save previous contig (fallback length)
                if current_id and current_id not in lengths:
                    lengths[current_id] = current_len

                parts = line[1:].split()
                current_id = parts[0]
                current_len = 0

                # Try MEGAHIT header format first
                m = re.search(r"len=(\d+)", line)
                if m:
                    lengths[current_id] = int(m.group(1))
            else:
                current_len += len(line)

        # Save last contig
        if current_id and current_id not in lengths:
            lengths[current_id] = current_len

    return lengths


def classify_kingdom(title):
    """Classify a BLAST subject title into kingdom/type."""
    t = title.lower()
    # Virus detection
    if any(w in t for w in ("virus", "viral", "phage", "viroid")):
        return "Virus"
    # Mycoplasma (check before bacteria - mycoplasma is a distinct category)
    if "ycoplasma" in t:
        return "Mycoplasma"
    # Bacteria
    if any(w in t for w in (
        "acteria", "acillus", "taphylococc", "seudomonas", "scherichia",
        "nterococcus", "lostridium", "treptococc", "almonella",
        "lebsiella", "aemophilus", "eisseria",
    )):
        return "Bacteria"
    # Fungi
    if any(w in t for w in (
        "ungi", "andida", "ryptococc", "accharomyces", "spergillus",
    )):
        return "Fungi"
    # Human
    if any(w in t for w in ("omo sapiens", "uman")):
        return "Human"
    return "Other"


def extract_organism_name(title, kingdom):
    """Extract a short organism name from BLAST subject title.

    For non-viral: returns genus + species (e.g., 'Escherichia coli')
    For viral: returns virus name up to first comma (e.g., 'Dengue virus 2')
    """
    name = title.strip('"')
    # Take everything before common suffixes
    for sep in [", complete genome", ", complete sequence",
                ", partial cds", " non structural"]:
        idx = name.find(sep)
        if idx > 0:
            name = name[:idx]
            break
    # Handle "chromosome X," pattern
    m = re.search(r' chromosome\b', name)
    if m:
        name = name[:m.start()]

    # For non-viral organisms, simplify to genus + species
    if kingdom != "Virus":
        words = name.split()
        if len(words) >= 2:
            name = f"{words[0]} {words[1]}"

    return name.strip()


def parse_blast_tsv(blast_path):
    """Read BLAST outfmt 6 TSV (with header) into list of dicts."""
    rows = []
    with open(blast_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            # Skip header
            if line.startswith("Query_ID") or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            rows.append({
                "query_id": fields[0],
                "subject_id": fields[1],
                "pident": float(fields[2]),
                "align_len": int(fields[3]),
                "mismatches": int(fields[4]),
                "gap_opens": int(fields[5]),
                "qstart": int(fields[6]),
                "qend": int(fields[7]),
                "sstart": int(fields[8]),
                "send": int(fields[9]),
                "evalue": float(fields[10]),
                "bitscore": float(fields[11]),
                # stitle may contain tabs — join remaining fields
                "stitle": "\t".join(fields[12:]) if len(fields) > 12 else "",
            })
    return rows


def get_hits_per_contig(rows, contig_lengths):
    """Get ALL hits per contig with tied-best detection and query coverage.

    Returns list of hit dicts sorted by contig length desc, bitscore desc
    within each contig. Each hit includes is_tied_best flag.

    Tied-best = same Bitscore, E-value, AND %Identity (within 0.001%)
    as the best hit for that contig.
    """
    by_contig = defaultdict(list)
    for r in rows:
        by_contig[r["query_id"]].append(r)

    all_hits = []
    for contig_id, hits in by_contig.items():
        # Sort: bitscore DESC, evalue ASC, pident DESC
        hits.sort(key=lambda h: (-h["bitscore"], h["evalue"], -h["pident"]))
        best = hits[0]

        for h in hits:
            clen = contig_lengths.get(contig_id, 0)
            cov_len = abs(h["qend"] - h["qstart"]) + 1
            qcov = (cov_len / clen * 100) if clen > 0 else 0.0

            # Tied-best: same bitscore, evalue, and pident (with float tolerance)
            is_tied = (h["bitscore"] == best["bitscore"]
                       and h["evalue"] == best["evalue"]
                       and abs(h["pident"] - best["pident"]) < 0.001)

            kingdom = classify_kingdom(h["stitle"])
            all_hits.append({
                "contig_id": contig_id,
                "contig_length": clen,
                "subject_id": h["subject_id"],
                "pident": h["pident"],
                "align_len": h["align_len"],
                "qcov": qcov,
                "evalue": h["evalue"],
                "bitscore": h["bitscore"],
                "kingdom": kingdom,
                "stitle": h["stitle"],
                "is_tied_best": is_tied,
            })

    # Sort: contig length desc, then contig_id for stability, bitscore desc within
    all_hits.sort(key=lambda h: (-h["contig_length"], h["contig_id"], -h["bitscore"]))
    return all_hits


def classify_contigs(all_hits, accession, min_coverage, accession_set=None):
    """Classify contigs into reference, other-viral, contaminant categories.

    A contig is classified by its BEST hit (highest bitscore).
    Hits below min_coverage are excluded from classification.

    accession_set: set of accession strings to match (includes segment accessions).
    Falls back to {accession} if not provided.

    Returns dict with:
      reference:    list of (contig_id, contig_length, [best_hits])
      other_viral:  list of (contig_id, contig_length, [best_hits])
      non_viral:    OrderedDict of kingdom -> list of (contig_id, contig_length, best_hit)
    """
    if accession_set is None:
        accession_set = {accession} if accession else set()
    # Group hits by contig
    by_contig = defaultdict(list)
    for h in all_hits:
        by_contig[h["contig_id"]].append(h)

    reference_contigs = []
    other_viral = []
    non_viral = defaultdict(list)

    for contig_id, hits in by_contig.items():
        # Filter to hits above coverage threshold
        valid_hits = [h for h in hits if h["qcov"] >= min_coverage]
        if not valid_hits:
            continue

        # Best hits = tied-best among valid hits
        # Re-sort valid hits and recalculate tied-best in case filtering changed order
        valid_hits.sort(key=lambda h: (-h["bitscore"], h["evalue"], -h["pident"]))
        vbest = valid_hits[0]
        best_hits = []
        for h in valid_hits:
            if (h["bitscore"] == vbest["bitscore"]
                    and h["evalue"] == vbest["evalue"]
                    and abs(h["pident"] - vbest["pident"]) < 0.001):
                best_hits.append(h)
            else:
                break

        # Check if reference accession (or segment accession) appears in ANY valid hit
        has_ref = any(accession_matches(h["subject_id"], accession_set)
                      for h in valid_hits) if accession_set else False

        if has_ref:
            reference_contigs.append((contig_id, vbest["contig_length"], best_hits))
        elif vbest["kingdom"] == "Virus":
            other_viral.append((contig_id, vbest["contig_length"], best_hits))
        else:
            kingdom = vbest["kingdom"]
            non_viral[kingdom].append((contig_id, vbest["contig_length"], vbest))

    # Sort by contig length desc
    reference_contigs.sort(key=lambda x: -x[1])
    other_viral.sort(key=lambda x: -x[1])
    for k in non_viral:
        non_viral[k].sort(key=lambda x: -x[1])

    return {
        "reference": reference_contigs,
        "other_viral": other_viral,
        "non_viral": dict(non_viral),
    }


def write_top_hits_tsv(all_hits, output_path):
    """Write only tied-best hits per contig (no Tied_Best column)."""
    with open(output_path, "w") as f:
        f.write("Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\t"
                "Alignment_Length\tQuery_Coverage\tE_value\tBit_Score\t"
                "Kingdom/Type\tSubject_Title\n")
        for h in all_hits:
            if not h["is_tied_best"]:
                continue
            f.write(
                f"{h['contig_id']}\t"
                f"{h['contig_length']}\t"
                f"{h['subject_id']}\t"
                f"{h['pident']:.2f}%\t"
                f"{h['align_len']}\t"
                f"{h['qcov']:.1f}%\t"
                f"{h['evalue']}\t"
                f"{h['bitscore']}\t"
                f"{h['kingdom']}\t"
                f"{h['stitle']}\n"
            )


def write_viral_hits_tsv(rows, contig_lengths, output_path):
    """Write viral hits with query coverage information."""
    with open(output_path, "w") as f:
        f.write("Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\t"
                "Alignment_Length\tQuery_Coverage\tE_value\tBit_Score\t"
                "Subject_Title\n")
        count = 0
        for r in rows:
            if classify_kingdom(r["stitle"]) == "Virus":
                clen = contig_lengths.get(r["query_id"], 0)
                cov_len = abs(r["qend"] - r["qstart"]) + 1
                qcov = (cov_len / clen * 100) if clen > 0 else 0.0
                f.write(
                    f"{r['query_id']}\t{clen}\t{r['subject_id']}\t"
                    f"{r['pident']:.2f}%\t{r['align_len']}\t{qcov:.1f}%\t"
                    f"{r['evalue']}\t{r['bitscore']}\t{r['stitle']}\n"
                )
                count += 1
    return count


def write_report_section(all_hits, contig_lengths, accession, min_coverage, output_path,
                         accession_set=None):
    """Write detailed contamination screening report."""
    classification = classify_contigs(all_hits, accession, min_coverage,
                                      accession_set=accession_set)

    # Get all contig IDs (including those with no hits)
    all_contig_ids = set(contig_lengths.keys())
    contigs_with_hits = set(h["contig_id"] for h in all_hits)
    no_hit_contigs = all_contig_ids - contigs_with_hits

    # Contigs whose ONLY hits are below coverage threshold
    by_contig = defaultdict(list)
    for h in all_hits:
        by_contig[h["contig_id"]].append(h)
    contigs_below_threshold = set()
    for cid, hits in by_contig.items():
        if all(h["qcov"] < min_coverage for h in hits):
            contigs_below_threshold.add(cid)

    with open(output_path, "w") as f:
        # Header
        f.write("BLAST CONTAMINATION SCREENING REPORT\n")
        f.write("=" * 70 + "\n")
        f.write(f"Reference accession: {accession}\n")
        f.write(f"Minimum query coverage threshold: {min_coverage:.0f}%\n")
        f.write(f"Total contigs: {len(all_contig_ids)}\n")
        f.write(f"Contigs with BLAST hits: {len(contigs_with_hits)}\n\n")

        # --- Reference virus ---
        ref_contigs = classification["reference"]
        f.write(f"REFERENCE VIRUS ({accession})\n")
        f.write("-" * 70 + "\n")
        if ref_contigs:
            total_bp = sum(cl for _, cl, _ in ref_contigs)
            f.write(f"  CONFIRMED: {len(ref_contigs)} contig(s), {total_bp:,} total bp\n")
            for cid, clen, best_hits in ref_contigs:
                bh = best_hits[0]
                f.write(f"    {cid} ({clen:,} bp): "
                        f"{bh['pident']:.2f}% identity, {bh['qcov']:.1f}% coverage\n")
        else:
            f.write(f"  WARNING: Reference {accession} NOT found in BLAST results\n")
        f.write("\n")

        # --- Other viral hits ---
        other_viral = classification["other_viral"]
        f.write("OTHER VIRAL HITS (potential co-infection/contamination)\n")
        f.write("-" * 70 + "\n")
        if other_viral:
            # Group by organism name for summary
            virus_groups = defaultdict(list)
            for cid, clen, best_hits in other_viral:
                bh = best_hits[0]
                vname = extract_organism_name(bh["stitle"], "Virus")
                virus_groups[vname].append((cid, clen, bh))

            # Summary line
            f.write(f"  {len(other_viral)} contig(s) from "
                    f"{len(virus_groups)} virus(es):\n\n")

            # Per-virus detail
            for vname, contigs in sorted(virus_groups.items(),
                                         key=lambda x: -sum(c[1] for c in x[1])):
                total_bp = sum(c[1] for c in contigs)
                f.write(f"  {vname}: {len(contigs)} contig(s), {total_bp:,} bp\n")
                f.write(f"    {'Contig_ID':<15} {'Length':>8}  {'%ID':>7}  {'%Cov':>7}\n")
                for cid, clen, bh in contigs:
                    f.write(f"    {cid:<15} {clen:>7,}  "
                            f"{bh['pident']:>6.2f}%  {bh['qcov']:>6.1f}%\n")
                f.write("\n")
        else:
            f.write(f"  None detected (>={min_coverage:.0f}% coverage threshold)\n\n")

        # --- Non-viral contaminants ---
        non_viral = classification["non_viral"]
        f.write("NON-VIRAL CONTAMINANTS\n")
        f.write("-" * 70 + "\n")

        # Display in order: Bacteria, Mycoplasma, Fungi, Human, Other
        display_order = ["Bacteria", "Mycoplasma", "Fungi", "Human", "Other"]
        has_any = False
        for kingdom in display_order:
            contigs_list = non_viral.get(kingdom, [])
            if not contigs_list:
                continue
            has_any = True

            # Sub-group by organism species
            org_groups = defaultdict(list)
            for cid, clen, bh in contigs_list:
                oname = extract_organism_name(bh["stitle"], kingdom)
                org_groups[oname].append((cid, clen, bh))

            total_contigs = len(contigs_list)
            f.write(f"\n  {kingdom}: {total_contigs} contig(s)\n")

            for oname, contigs in sorted(org_groups.items(),
                                         key=lambda x: -len(x[1])):
                cids = [c[0] for c in contigs]
                f.write(f"    {oname}: {len(contigs)} contig(s) "
                        f"({', '.join(cids)})\n")

        if not has_any:
            f.write(f"  None detected (>={min_coverage:.0f}% coverage threshold)\n")
        f.write("\n")

        # --- No hits / below threshold ---
        no_hit_or_low = no_hit_contigs | contigs_below_threshold
        f.write("CONTIGS WITH NO SIGNIFICANT HITS\n")
        f.write("-" * 70 + "\n")
        if no_hit_or_low:
            sorted_ids = sorted(no_hit_or_low,
                                key=lambda c: -contig_lengths.get(c, 0))
            f.write(f"  {len(no_hit_or_low)} contig(s):\n")
            for cid in sorted_ids:
                clen = contig_lengths.get(cid, 0)
                if cid in contigs_below_threshold:
                    # Show what the hit was (below threshold)
                    best = by_contig[cid][0]
                    f.write(f"    {cid} ({clen:,} bp) — "
                            f"hits below {min_coverage:.0f}% cov threshold "
                            f"(best: {best['qcov']:.1f}% cov to "
                            f"{extract_organism_name(best['stitle'], best['kingdom'])})\n")
                else:
                    f.write(f"    {cid} ({clen:,} bp) — no BLAST hits\n")
        else:
            f.write("  All contigs have significant BLAST hits\n")


def main():
    parser = argparse.ArgumentParser(
        description="Parse BLAST results with contamination screening and tied-hit support"
    )
    parser.add_argument("blast_tsv", help="BLAST results TSV (with header)")
    parser.add_argument("contigs_fasta", help="Contigs FASTA file")
    parser.add_argument("--accession", default="",
                        help="Reference accession to check for")
    parser.add_argument("--top-hits", default="",
                        help="Output all hits per contig TSV path")
    parser.add_argument("--viral-hits", default="",
                        help="Output viral hits TSV path")
    parser.add_argument("--report", default="",
                        help="Output contamination report text path")
    parser.add_argument("--min-coverage", type=float, default=10.0,
                        help="Minimum query coverage %% to count as contamination (default: 10)")
    args = parser.parse_args()

    # Parse inputs
    contig_lengths = parse_contig_lengths(args.contigs_fasta)
    rows = parse_blast_tsv(args.blast_tsv)

    print(f"Parsed {len(rows)} BLAST hits for {len(contig_lengths)} contigs",
          file=sys.stderr)

    # Get ALL hits per contig (not just tied-best)
    all_hits = get_hits_per_contig(rows, contig_lengths)

    # Count stats
    contigs_with_hits = len(set(h["contig_id"] for h in all_hits))
    tied_best_count = sum(1 for h in all_hits if h["is_tied_best"])
    print(f"Hits per contig: {len(all_hits)} total entries for "
          f"{contigs_with_hits} contigs "
          f"({tied_best_count} tied-best)", file=sys.stderr)

    # Resolve accession to include segment accessions for segmented viruses
    accession_set = resolve_accessions(args.accession) if args.accession else set()

    # Accession check
    if args.accession:
        has_ref = any(accession_matches(h["subject_id"], accession_set) for h in all_hits)
        if has_ref:
            if len(accession_set) > 1:
                print(f"REFERENCE MATCH: {args.accession} confirmed via segment accessions "
                      f"({', '.join(sorted(accession_set - {args.accession}))})",
                      file=sys.stderr)
            else:
                print(f"REFERENCE MATCH: {args.accession} confirmed in BLAST hits",
                      file=sys.stderr)
        else:
            print(f"WARNING: Reference {args.accession} not found in BLAST hits",
                  file=sys.stderr)

    # Write outputs
    if args.top_hits:
        write_top_hits_tsv(all_hits, args.top_hits)
        print(f"Wrote top hits: {args.top_hits}", file=sys.stderr)

    if args.viral_hits:
        viral_count = write_viral_hits_tsv(rows, contig_lengths, args.viral_hits)
        print(f"Wrote {viral_count} viral hits: {args.viral_hits}", file=sys.stderr)

    if args.report:
        write_report_section(all_hits, contig_lengths, args.accession,
                             args.min_coverage, args.report,
                             accession_set=accession_set)
        print(f"Wrote report: {args.report}", file=sys.stderr)

    # Also print top hits to stdout for quick viewing
    if not args.top_hits and not args.report:
        write_top_hits_tsv(all_hits, "/dev/stdout")


if __name__ == "__main__":
    main()
