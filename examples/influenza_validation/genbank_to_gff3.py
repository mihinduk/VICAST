#!/usr/bin/env python3
"""
Convert GenBank files to GFF3 format for VICAST/SnpEff.

Creates properly formatted GFF3 with gene/mRNA/CDS hierarchy.
"""

import sys
from pathlib import Path
from Bio import SeqIO


def genbank_to_gff3(gb_files, output_gff):
    """Convert GenBank files to GFF3."""

    with open(output_gff, 'w') as out:
        out.write("##gff-version 3\n")

        for gb_file in sorted(gb_files):
            print(f"Processing: {gb_file}")

            for record in SeqIO.parse(gb_file, "genbank"):
                seqid = record.id
                seq_len = len(record.seq)

                # Write sequence region
                out.write(f"##sequence-region {seqid} 1 {seq_len}\n")

                # Track gene IDs to avoid duplicates
                gene_count = {}

                for feature in record.features:
                    if feature.type == "CDS":
                        # Extract info
                        start = int(feature.location.start) + 1  # 1-based
                        end = int(feature.location.end)
                        strand = '+' if feature.location.strand == 1 else '-'

                        # Get gene name and product
                        gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
                        product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                        protein_id = feature.qualifiers.get('protein_id', [''])[0]

                        # Handle duplicates (e.g., M1/M2 on same segment)
                        if gene_name in gene_count:
                            gene_count[gene_name] += 1
                            unique_id = f"{gene_name}_{gene_count[gene_name]}"
                        else:
                            gene_count[gene_name] = 1
                            unique_id = gene_name

                        # Handle spliced features (join locations)
                        if hasattr(feature.location, 'parts'):
                            # Spliced gene - write each exon
                            parts = list(feature.location.parts)

                            # Gene feature spans full range
                            gene_start = min(int(p.start) + 1 for p in parts)
                            gene_end = max(int(p.end) for p in parts)

                            out.write(f"{seqid}\tNCBI\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
                            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")

                            out.write(f"{seqid}\tNCBI\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
                            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")

                            # Write CDS parts
                            for i, part in enumerate(parts):
                                p_start = int(part.start) + 1
                                p_end = int(part.end)
                                phase = 0 if i == 0 else (3 - ((p_start - parts[0].start) % 3)) % 3

                                out.write(f"{seqid}\tNCBI\tCDS\t{p_start}\t{p_end}\t.\t{strand}\t{phase}\t")
                                out.write(f"ID=cds_{unique_id}_{i+1};Parent=mRNA_{unique_id};")
                                out.write(f"Name={gene_name};gene={gene_name};product={product}")
                                if protein_id:
                                    out.write(f";protein_id={protein_id}")
                                out.write("\n")
                        else:
                            # Simple unspliced gene
                            out.write(f"{seqid}\tNCBI\tgene\t{start}\t{end}\t.\t{strand}\t.\t")
                            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")

                            out.write(f"{seqid}\tNCBI\tmRNA\t{start}\t{end}\t.\t{strand}\t.\t")
                            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")

                            out.write(f"{seqid}\tNCBI\tCDS\t{start}\t{end}\t.\t{strand}\t0\t")
                            out.write(f"ID=cds_{unique_id};Parent=mRNA_{unique_id};")
                            out.write(f"Name={gene_name};gene={gene_name};product={product}")
                            if protein_id:
                                out.write(f";protein_id={protein_id}")
                            out.write("\n")

    print(f"Created: {output_gff}")


if __name__ == "__main__":
    gb_files = sorted(Path(".").glob("NC_*.gb"))
    genbank_to_gff3(gb_files, "influenza_A_California_2009.gff3")
