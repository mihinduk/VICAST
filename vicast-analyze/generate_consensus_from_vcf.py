#!/usr/bin/env python3
import sys
import pandas as pd
from Bio import SeqIO
import re

# Read annotated VCF
vcf_file = sys.argv[1]
ref_file = sys.argv[2]
output_prefix = sys.argv[3]

# Parse VCF
variants = []
with open(vcf_file, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        chrom, pos, id_, ref, alt, qual, filt, info = parts[:8]
        
        # Extract depth and frequency
        dp = int(re.search(r'DP=(\d+)', info).group(1))
        af = float(re.search(r'AF=([\d\.]+)', info).group(1))
        qual_score = float(qual)
        
        # Filter: Q>=1000, depth>=20, freq>=0.50 for consensus
        if qual_score >= 1000 and dp >= 20 and af >= 0.50:
            variants.append({
                'POS': int(pos),
                'REF': ref,
                'ALT': alt,
                'QUAL': qual_score,
                'DP': dp,
                'AF': af,
                'INFO': info
            })

print(f"Found {len(variants)} high-quality consensus variants")

# Read reference
ref_seq = str(next(SeqIO.parse(ref_file, 'fasta')).seq)
ref_list = list(ref_seq)

# Apply variants
for var in sorted(variants, key=lambda x: x['POS']):
    pos_idx = var['POS'] - 1  # 0-based
    if ref_list[pos_idx] == var['REF']:
        ref_list[pos_idx] = var['ALT']
        print(f"  Applied: {var['POS']} {var['REF']}>{var['ALT']} (AF={var['AF']:.2%})")

consensus = ''.join(ref_list)

# Write consensus
with open(f"{output_prefix}.fasta", 'w') as f:
    f.write(f">{output_prefix}_consensus\n{consensus}\n")

print(f"\nWrote consensus: {output_prefix}.fasta ({len(consensus)} bp)")
