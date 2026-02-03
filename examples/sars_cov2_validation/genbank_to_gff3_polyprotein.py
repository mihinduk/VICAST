#!/usr/bin/env python3
"""
Convert SARS-CoV-2 GenBank to GFF3 with mature peptide annotations.

This script creates a GFF3 that includes:
- Standard genes (S, E, M, N, ORF3a, etc.)
- ORF1ab polyprotein with ribosomal frameshift
- All 16 mature peptides (nsp1-16) as child features

This showcases VICAST's polyprotein annotation capability.
"""

import sys
from pathlib import Path
from Bio import SeqIO


def convert_sars_cov2_to_gff3(gb_file, output_gff, output_fasta=None):
    """Convert SARS-CoV-2 GenBank to GFF3 with polyprotein annotation."""

    record = list(SeqIO.parse(gb_file, "genbank"))[0]
    seqid = record.id
    seq_len = len(record.seq)

    print(f"Processing: {seqid} ({seq_len} bp)")

    with open(output_gff, 'w') as out:
        out.write("##gff-version 3\n")
        out.write(f"##sequence-region {seqid} 1 {seq_len}\n")

        # Track what we've written
        genes_written = set()
        cds_count = 0
        mat_peptide_count = 0

        # First pass: collect all CDS and mat_peptide features
        cds_features = []
        mat_peptides = []

        for feature in record.features:
            if feature.type == "CDS":
                cds_features.append(feature)
            elif feature.type == "mat_peptide":
                mat_peptides.append(feature)

        print(f"  Found {len(cds_features)} CDS features")
        print(f"  Found {len(mat_peptides)} mature peptide features")

        # Process CDS features
        for feature in cds_features:
            gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
            product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
            protein_id = feature.qualifiers.get('protein_id', [''])[0]

            # Escape special characters
            product = product.replace(';', '%3B').replace('=', '%3D')

            strand = '+' if feature.location.strand == 1 else '-'

            # Handle compound locations (ribosomal frameshift)
            if hasattr(feature.location, 'parts'):
                parts = list(feature.location.parts)
                gene_start = min(int(p.start) + 1 for p in parts)
                gene_end = max(int(p.end) for p in parts)
            else:
                gene_start = int(feature.location.start) + 1
                gene_end = int(feature.location.end)

            # Skip ORF1a polyprotein (pp1a - we only want pp1ab with frameshift)
            # ORF1a is the non-frameshifted product, ORF1ab includes the frameshift
            if gene_name == "ORF1ab":
                if "orf1a polyprotein" in product.lower() and "orf1ab" not in product.lower():
                    continue
                if "pp1a" in product.lower() and "pp1ab" not in product.lower():
                    continue

            # Create unique ID
            if gene_name in genes_written:
                unique_id = f"{gene_name}_{len([g for g in genes_written if g.startswith(gene_name)]) + 1}"
            else:
                unique_id = gene_name

            genes_written.add(gene_name)
            cds_count += 1

            # Gene feature
            out.write(f"{seqid}\tNCBI\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
            out.write(f"ID=gene_{unique_id};Name={gene_name}\n")

            # mRNA feature
            out.write(f"{seqid}\tNCBI\tmRNA\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t")
            out.write(f"ID=mRNA_{unique_id};Parent=gene_{unique_id};Name={gene_name};gene={gene_name}\n")

            # CDS feature(s)
            if hasattr(feature.location, 'parts'):
                parts = list(feature.location.parts)
                cumulative_len = 0
                for i, part in enumerate(parts):
                    p_start = int(part.start) + 1
                    p_end = int(part.end)
                    phase = (3 - (cumulative_len % 3)) % 3
                    cumulative_len += (p_end - p_start + 1)

                    out.write(f"{seqid}\tNCBI\tCDS\t{p_start}\t{p_end}\t.\t{strand}\t{phase}\t")
                    out.write(f"ID=cds_{unique_id}_{i+1};Parent=mRNA_{unique_id};")
                    out.write(f"Name={gene_name};gene={gene_name};product={product}")
                    if protein_id:
                        out.write(f";protein_id={protein_id}")
                    out.write("\n")
            else:
                out.write(f"{seqid}\tNCBI\tCDS\t{gene_start}\t{gene_end}\t.\t{strand}\t0\t")
                out.write(f"ID=cds_{unique_id};Parent=mRNA_{unique_id};")
                out.write(f"Name={gene_name};gene={gene_name};product={product}")
                if protein_id:
                    out.write(f";protein_id={protein_id}")
                out.write("\n")

        # Process mature peptides (nsp1-16)
        # These are children of ORF1ab
        # Track seen coordinates to avoid duplicates (pp1a and pp1ab have same peptides)
        seen_coords = set()
        import re

        for feature in mat_peptides:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            # Skip duplicates (same coordinates)
            coord_key = (start, end)
            if coord_key in seen_coords:
                continue
            seen_coords.add(coord_key)

            mat_peptide_count += 1

            gene_name = feature.qualifiers.get('gene', ['ORF1ab'])[0]
            product = feature.qualifiers.get('product', ['mature peptide'])[0]
            protein_id = feature.qualifiers.get('protein_id', [''])[0]
            note = feature.qualifiers.get('note', [''])[0]

            # Extract nsp name - prioritize product over note
            # (note can have misleading text like "former nsp1" for nsp3)
            nsp_name = None

            # First check if product is exactly "nspN" or starts with "nspN"
            if product:
                nsp_match = re.match(r'^(nsp\d+)\b', product, re.IGNORECASE)
                if nsp_match:
                    nsp_name = nsp_match.group(1).lower()

            # If not found in product, check note but only for primary nsp designation
            # Look for patterns like "nsp5A_3CLpro" at the start of note
            if not nsp_name and note:
                # Match nsp at word boundary but not after "former"
                nsp_match = re.search(r'(?<!former\s)(nsp\d+)', note, re.IGNORECASE)
                if nsp_match:
                    nsp_name = nsp_match.group(1).lower()

            # Map well-known products to nsp numbers
            if not nsp_name:
                product_to_nsp = {
                    'leader protein': 'nsp1',
                    '3c-like proteinase': 'nsp5',
                    'rna-dependent rna polymerase': 'nsp12',
                    'helicase': 'nsp13',
                    "3'-to-5' exonuclease": 'nsp14',
                    'endornase': 'nsp15',
                    "2'-o-ribose methyltransferase": 'nsp16',
                }
                nsp_name = product_to_nsp.get(product.lower(), None)

            # Fall back to cleaned product name
            if not nsp_name:
                nsp_name = product.split()[0] if product else 'unknown'

            # Clean up for ID
            safe_id = nsp_name.replace(' ', '_').replace('-', '_').replace("'", "")
            safe_id = ''.join(c for c in safe_id if c.isalnum() or c == '_')

            product_escaped = product.replace(';', '%3B').replace('=', '%3D')

            strand = '+' if feature.location.strand == 1 else '-'

            # Write mature peptide as child of ORF1ab mRNA
            out.write(f"{seqid}\tNCBI\tmature_peptide\t{start}\t{end}\t.\t{strand}\t.\t")
            out.write(f"ID=mat_{safe_id};Parent=mRNA_ORF1ab;")
            out.write(f"Name={nsp_name};gene={gene_name};product={product_escaped}")
            if protein_id:
                out.write(f";protein_id={protein_id}")
            out.write("\n")

    print(f"  Wrote {cds_count} CDS features")
    print(f"  Wrote {mat_peptide_count} mature peptide features")
    print(f"Created: {output_gff}")

    # Write FASTA if requested
    if output_fasta:
        SeqIO.write([record], output_fasta, "fasta")
        print(f"Created: {output_fasta}")

    return cds_count, mat_peptide_count


if __name__ == "__main__":
    gb_file = "NC_045512.gb"
    gff_file = "NC_045512.gff3"
    fasta_file = "NC_045512.fasta"

    if not Path(gb_file).exists():
        print(f"Error: {gb_file} not found")
        sys.exit(1)

    convert_sars_cov2_to_gff3(gb_file, gff_file, fasta_file)
