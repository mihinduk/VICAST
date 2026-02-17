#!/usr/bin/env python3
"""
Convert a VCF file to TSV with INFO fields expanded into individual columns.

Parses the INFO column (semicolon-separated key=value pairs) and creates
one column per unique INFO key. Flag fields (no value) are marked as "Yes".

Designed for lofreq _vars.filt.vcf output but works with any standard VCF.

Usage:
    python vcf_to_tsv.py input_vars.filt.vcf
    python vcf_to_tsv.py input_vars.filt.vcf -o output.tsv
    python vcf_to_tsv.py input_vars.filt.vcf --split-dp4
"""

import argparse
import os
import sys


def parse_info(info_str):
    """Parse a VCF INFO string into a dict of key->value pairs.
    Flags (no =) get value 'Yes'."""
    fields = {}
    if info_str == "." or not info_str:
        return fields
    for entry in info_str.split(";"):
        if "=" in entry:
            key, val = entry.split("=", 1)
            fields[key] = val
        else:
            fields[entry] = "Yes"
    return fields


def main():
    parser = argparse.ArgumentParser(
        description="Convert VCF to TSV with INFO fields as separate columns"
    )
    parser.add_argument("vcf", help="Input VCF file (e.g. sample_vars.filt.vcf)")
    parser.add_argument("-o", "--output", help="Output TSV file (default: stdout)")
    parser.add_argument(
        "--split-dp4",
        action="store_true",
        help="Split DP4 into REF_FWD, REF_REV, ALT_FWD, ALT_REV columns",
    )
    args = parser.parse_args()

    # First pass: collect all INFO keys in order of first appearance
    info_keys = []
    seen_keys = set()
    header_line = None
    records = []

    with open(args.vcf) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM") or line.startswith("#"):
                header_line = line
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                continue
            info = parse_info(cols[7])
            for key in info:
                if key not in seen_keys:
                    info_keys.append(key)
                    seen_keys.add(key)
            records.append(cols)

    if not records:
        print("No variant records found in VCF.", file=sys.stderr)
        sys.exit(1)

    # Build output column order
    base_cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]

    # Expand DP4 if requested
    dp4_idx = None
    expanded_keys = []
    for key in info_keys:
        if key == "DP4" and args.split_dp4:
            dp4_idx = len(expanded_keys)
            expanded_keys.extend(["REF_FWD", "REF_REV", "ALT_FWD", "ALT_REV"])
        else:
            expanded_keys.append(key)

    out_header = base_cols + expanded_keys

    out = open(args.output, "w") if args.output else sys.stdout

    try:
        out.write("\t".join(out_header) + "\n")

        for cols in records:
            row = cols[:7]  # CHROM through FILTER
            info = parse_info(cols[7])

            for key in info_keys:
                val = info.get(key, "")
                if key == "DP4" and args.split_dp4:
                    if val:
                        parts = val.split(",")
                        row.extend(parts[:4] if len(parts) >= 4 else parts + [""] * (4 - len(parts)))
                    else:
                        row.extend(["", "", "", ""])
                else:
                    row.append(val)

            out.write("\t".join(row) + "\n")
    finally:
        if args.output:
            out.close()

    if args.output:
        dest = os.path.abspath(args.output)
        # Map Docker container path back to host path if HOST_PWD is set
        host_pwd = os.environ.get("HOST_PWD", "")
        if host_pwd and dest.startswith("/data/"):
            dest = host_pwd + dest[5:]  # replace /data with host_pwd
    else:
        dest = "stdout"
    n = len(records)
    print(f"Converted {n} variants with {len(expanded_keys)} INFO columns -> {dest}", file=sys.stderr)


if __name__ == "__main__":
    main()
