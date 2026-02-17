#!/usr/bin/env python3
import pandas as pd
import re
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description='Parse SnpEff TSV with quality filtering')
    parser.add_argument('input_file', help='Input SnpEff TSV file')
    parser.add_argument('output_file', help='Output filtered TSV file')
    parser.add_argument('-q', '--quality', type=float, default=1000, help='Quality score cutoff (default: 1000)')
    parser.add_argument('-d', '--depth', type=int, default=200, help='Minimum depth cutoff (default: 200)')
    parser.add_argument('-f', '--freq', type=float, default=0.01, help='Minimum allele frequency cutoff (default: 0.01)')
    
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    quality_cutoff = args.quality
    depth_cutoff = args.depth
    freq_cutoff = args.freq

    print(f"Processing SnpEff TSV file: {input_file}")
    print(f"Quality cutoff: {quality_cutoff}")
    print(f"Depth cutoff: {depth_cutoff}")
    print(f"Allele frequency cutoff: {freq_cutoff}")

    # Read the SnpEff TSV file - handle the #CHROM header line
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Find the header line that starts with #CHROM and remove the leading #
    header_line = None
    data_lines = []
    for line in lines:
        if line.startswith('#CHROM'):
            header_line = line[1:].strip()  # Remove leading # and whitespace
        elif not line.startswith('##'):  # Skip SnpEff metadata lines
            data_lines.append(line.strip())

    if not header_line:
        print("❌ Error: Could not find #CHROM header line")
        sys.exit(1)

    # Reconstruct the TSV content
    tsv_content = header_line + '\n' + '\n'.join(data_lines)

    import io
    df = pd.read_csv(io.StringIO(tsv_content), sep='\t')

    print(f"Found {len(df)} total variants")

    # Extract AF from INFO field
    def extract_af(info):
        match = re.search(r'AF=([0-9.]+)', str(info))
        return float(match.group(1)) if match else 0.0

    # Extract DP from INFO field  
    def extract_dp(info):
        match = re.search(r'DP=([0-9]+)', str(info))
        return int(match.group(1)) if match else 0

    # Extract DP4 from INFO field for strand bias information
    def extract_dp4(info):
        match = re.search(r'DP4=([0-9,]+)', str(info))
        return match.group(1) if match else '0,0,0,0'

    # Extract SB (strand bias) from INFO field
    def extract_sb(info):
        match = re.search(r'SB=([0-9]+)', str(info))
        return match.group(1) if match else '0'

    # Parse the data
    df['Allele_Frequency'] = df['INFO'].apply(extract_af)
    df['Total_Depth'] = df['INFO'].apply(extract_dp)
    df['DP4'] = df['INFO'].apply(extract_dp4) 
    df['strand_bias'] = df['INFO'].apply(extract_sb)

    print("Applying quality, depth, and frequency filtering...")
    
    # Apply quality, depth, and frequency filtering
    filtered_df = df[
        (df['QUAL'] >= quality_cutoff) & 
        (df['Total_Depth'] >= depth_cutoff) &
        (df['Allele_Frequency'] >= freq_cutoff)
    ].copy()
    
    print(f"After filtering (QUAL >= {quality_cutoff}, depth >= {depth_cutoff}, freq >= {freq_cutoff}): {len(filtered_df)} variants")

    if len(filtered_df) == 0:
        print("⚠️  Warning: No variants passed filtering. Consider lowering cutoffs.")
        
    # Create output format using the existing annotation columns
    output_df = pd.DataFrame({
        'CHROM': filtered_df['CHROM'],
        'POS': filtered_df['POS'],
        'ID': filtered_df['ID'],
        'REF': filtered_df['REF'],
        'ALT': filtered_df['ALT'],
        'QUAL': filtered_df['QUAL'],
        'INFO': filtered_df['INFO'],
        'Total_Depth': filtered_df['Total_Depth'],
        'Allele_Frequency': filtered_df['Allele_Frequency'],
        'strand_bias': filtered_df['strand_bias'],
        'DP4': filtered_df['DP4'],
        'EFFECT': filtered_df['EFFECT'],
        'PUTATIVE_IMPACT': filtered_df['PUTATIVE_IMPACT'],
        'GENE_NAME': filtered_df['GENE_NAME'],
        'GENE_ID': filtered_df['GENE_ID'],
        'FEATURE_TYPE': filtered_df['FEATURE_TYPE'],
        'FEATURE_ID': filtered_df['FEATURE_ID'],
        'TRANSCRIPT_TYPE': filtered_df['TRANSCRIPT_TYPE'],
        'HGVSc': filtered_df['HGVSc'],
        'HGVSp': filtered_df['HGVSp'],
        'cDNA_POSITION_AND_LENGTH': filtered_df['cDNA_POSITION_AND_LENGTH'],
        'CDS_POSITION_AND_LENGTH': filtered_df['CDS_POSITION_AND_LENGTH'],
        'PROTEIN_POSITION_AND_LENGTH': filtered_df['PROTEIN_POSITION_AND_LENGTH'],
        'ERROR': filtered_df['ERROR']
    })

    # Save the parsed file
    output_df.to_csv(output_file, sep='\t', index=False)
    print(f"✅ Successfully parsed and filtered {len(output_df)} mutations to {output_file}")
    print("   - Applied quality and depth filtering")
    print("   - Preserved all SnpEff annotation columns")
    print("   - Ready for visualization with proper mutation type coloring")

if __name__ == "__main__":
    main()