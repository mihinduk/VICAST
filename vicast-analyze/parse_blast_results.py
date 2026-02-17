#!/usr/bin/env python3
"""
Parse BLAST results into top hits summary with tied-hit support.

Replaces the awk-based top hits processing in viral_diagnostic.sh.
Handles multiple statistically indistinguishable hits (same %ID, E-value,
Bitscore) and flags whether the user's reference accession appears.

Usage:
    python parse_blast_results.py blast_all.tsv contigs.fa --accession NC_001474.2
    python parse_blast_results.py blast_all.tsv contigs.fa --accession NC_001474.2 \
        --top-hits out_top_hits.tsv --viral-hits out_viral.tsv --report out_report.txt
"""

import argparse
import re
import sys
from collections import defaultdict


def parse_contig_lengths(fasta_path):
    """Parse contig lengths from MEGAHIT FASTA headers (len=XXXX)."""
    lengths = {}
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split()
                contig_id = parts[0]
                m = re.search(r"len=(\d+)", line)
                if m:
                    lengths[contig_id] = int(m.group(1))
    return lengths


def classify_kingdom(title):
    """Classify a BLAST subject title into kingdom/type."""
    t = title.lower()
    if any(w in t for w in ("virus", "viral", "phage", "viroid")):
        return "Virus"
    if "ycoplasma" in t:
        return "Mycoplasma"
    if "acteria" in t or "acillus" in t or "taphylococc" in t or "seudomonas" in t or "scherichia" in t:
        return "Bacteria"
    if "ungi" in t or "andida" in t or "ryptococc" in t or "accharomyces" in t:
        return "Fungi"
    if "omo sapiens" in t or "uman" in t:
        return "Human"
    return "Unknown"


def parse_blast_tsv(blast_path):
    """Read BLAST outfmt 6 TSV (with header) into list of dicts."""
    rows = []
    with open(blast_path) as f:
        header = None
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            # Skip header
            if line.startswith("Query_ID") or line.startswith("#"):
                header = line.split("\t")
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
                "stitle": fields[12] if len(fields) > 12 else "",
            })
    return rows


def get_top_hits(rows, contig_lengths):
    """Get top hits per contig, including ALL tied hits.

    Tied = same Bitscore, E-value, AND %Identity as the best hit for that contig.
    """
    # Group by query_id
    by_contig = defaultdict(list)
    for r in rows:
        by_contig[r["query_id"]].append(r)

    top_hits = []
    for contig_id, hits in by_contig.items():
        # Sort: bitscore DESC, evalue ASC, pident DESC
        hits.sort(key=lambda h: (-h["bitscore"], h["evalue"], -h["pident"]))
        best = hits[0]

        # Collect all tied hits
        for h in hits:
            if (h["bitscore"] == best["bitscore"]
                    and h["evalue"] == best["evalue"]
                    and h["pident"] == best["pident"]):
                clen = contig_lengths.get(contig_id, 0)
                cov_len = abs(h["qend"] - h["qstart"]) + 1
                qcov = (cov_len / clen * 100) if clen > 0 else 0.0
                kingdom = classify_kingdom(h["stitle"])
                top_hits.append({
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
                    "is_tied": True,
                })
            else:
                break  # No more ties

    # Sort by contig length descending (largest contigs first)
    top_hits.sort(key=lambda h: -h["contig_length"])
    return top_hits


def check_accession_match(top_hits, accession):
    """Check if user's accession appears in top hits and generate messages."""
    messages = []
    accession_found = False
    same_virus_found = False

    # Get the virus name from the accession's stitle if it appears anywhere
    accession_titles = set()
    for h in top_hits:
        if accession in h["subject_id"]:
            accession_found = True
            accession_titles.add(h["stitle"])

    if accession_found:
        messages.append(f"REFERENCE MATCH: {accession} confirmed among top BLAST hits")
    else:
        # Check if any top hit is the same virus species
        # Look for the virus name pattern from the accession
        for h in top_hits:
            title = h["stitle"].lower()
            # Check common virus name patterns
            if accession.startswith("NC_001474"):  # DENV-2
                if "dengue" in title and ("type 2" in title or "serotype 2" in title):
                    same_virus_found = True
                    break
            elif accession.startswith("NC_045512"):  # SARS-CoV-2
                if "sars-cov-2" in title or "severe acute respiratory syndrome coronavirus 2" in title:
                    same_virus_found = True
                    break

        if same_virus_found:
            messages.append(
                f"NOTE: Reference {accession} not the top hit, but same virus species detected in results"
            )
        elif top_hits:
            messages.append(
                f"WARNING: Reference {accession} not found in BLAST results"
            )

    return messages


def write_top_hits_tsv(top_hits, output_path):
    """Write top hits TSV with tied hits included."""
    with open(output_path, "w") as f:
        f.write("Contig_ID\tContig_Length\tSubject_ID\tPercent_Identity\t"
                "Alignment_Length\tQuery_Coverage\tE_value\tBit_Score\t"
                "Kingdom/Type\tSubject_Title\n")
        for h in top_hits:
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


def write_viral_hits_tsv(rows, output_path):
    """Filter BLAST rows for viral hits and write TSV."""
    with open(output_path, "w") as f:
        f.write("Query_ID\tSubject_ID\tPercent_Identity\tAlignment_Length\t"
                "Mismatches\tGap_Opens\tQuery_Start\tQuery_End\t"
                "Subject_Start\tSubject_End\tE_value\tBit_Score\tSubject_Title\n")
        count = 0
        for r in rows:
            if classify_kingdom(r["stitle"]) == "Virus":
                f.write(
                    f"{r['query_id']}\t{r['subject_id']}\t{r['pident']}\t"
                    f"{r['align_len']}\t{r['mismatches']}\t{r['gap_opens']}\t"
                    f"{r['qstart']}\t{r['qend']}\t{r['sstart']}\t{r['send']}\t"
                    f"{r['evalue']}\t{r['bitscore']}\t{r['stitle']}\n"
                )
                count += 1
    return count


def write_report_section(top_hits, accession_messages, output_path):
    """Write the BLAST results section for the diagnostic report."""
    with open(output_path, "w") as f:
        # Accession matching messages
        for msg in accession_messages:
            f.write(f"{msg}\n")
        f.write("\n")

        if not top_hits:
            f.write("No significant BLAST hits found.\n")
            return

        # Top hits table
        f.write("TOP HITS BY CONTIG (largest contigs first, tied hits grouped):\n")
        f.write("=" * 80 + "\n\n")

        current_contig = None
        for h in top_hits:
            if h["contig_id"] != current_contig:
                if current_contig is not None:
                    f.write("\n")
                current_contig = h["contig_id"]
                f.write(f"{h['contig_id']} ({h['contig_length']} bp):\n")

            f.write(
                f"  {h['pident']:.2f}%  {h['qcov']:5.1f}% cov  "
                f"E={h['evalue']}  Score={h['bitscore']}  "
                f"[{h['kingdom']}]  {h['stitle'][:70]}\n"
            )

        f.write("\n")

        # Kingdom summary
        kingdom_counts = defaultdict(int)
        seen_contigs = set()
        for h in top_hits:
            if h["contig_id"] not in seen_contigs:
                kingdom_counts[h["kingdom"]] += 1
                seen_contigs.add(h["contig_id"])

        f.write("CONTAMINATION SUMMARY (unique contigs):\n")
        for k, count in sorted(kingdom_counts.items(), key=lambda x: -x[1]):
            f.write(f"  {k}: {count} contig{'s' if count > 1 else ''}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Parse BLAST results with tied-hit support and accession matching"
    )
    parser.add_argument("blast_tsv", help="BLAST results TSV (with header)")
    parser.add_argument("contigs_fasta", help="Contigs FASTA (MEGAHIT format with len= headers)")
    parser.add_argument("--accession", default="", help="User's reference accession to check for")
    parser.add_argument("--top-hits", default="", help="Output top hits TSV path")
    parser.add_argument("--viral-hits", default="", help="Output viral hits TSV path")
    parser.add_argument("--report", default="", help="Output report section text path")
    args = parser.parse_args()

    # Parse inputs
    contig_lengths = parse_contig_lengths(args.contigs_fasta)
    rows = parse_blast_tsv(args.blast_tsv)

    print(f"Parsed {len(rows)} BLAST hits for {len(contig_lengths)} contigs", file=sys.stderr)

    # Get top hits with ties
    top_hits = get_top_hits(rows, contig_lengths)

    # Count unique contigs with hits
    contigs_with_hits = len(set(h["contig_id"] for h in top_hits))
    tied_groups = sum(1 for i, h in enumerate(top_hits)
                      if i == 0 or h["contig_id"] != top_hits[i-1]["contig_id"])
    print(f"Top hits: {len(top_hits)} entries for {contigs_with_hits} contigs "
          f"({len(top_hits) - tied_groups} tied hits)", file=sys.stderr)

    # Check accession match
    accession_messages = []
    if args.accession:
        accession_messages = check_accession_match(top_hits, args.accession)
        for msg in accession_messages:
            print(msg, file=sys.stderr)

    # Write outputs
    if args.top_hits:
        write_top_hits_tsv(top_hits, args.top_hits)
        print(f"Wrote top hits: {args.top_hits}", file=sys.stderr)

    if args.viral_hits:
        viral_count = write_viral_hits_tsv(rows, args.viral_hits)
        print(f"Wrote {viral_count} viral hits: {args.viral_hits}", file=sys.stderr)

    if args.report:
        write_report_section(top_hits, accession_messages, args.report)
        print(f"Wrote report section: {args.report}", file=sys.stderr)

    # Also print top hits to stdout for quick viewing
    if not args.top_hits and not args.report:
        write_top_hits_tsv(top_hits, "/dev/stdout")


if __name__ == "__main__":
    main()
