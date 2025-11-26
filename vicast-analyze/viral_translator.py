#\!/usr/bin/env python3
"""
Custom viral translator for polyprotein cleavage products and spliced genes
Handles the clever tricks viruses use that BioPython doesn't handle well
"""

def get_codon_table():
    """Standard genetic code table"""
    return {
        # Phenylalanine
        'TTT': 'F', 'TTC': 'F',
        # Leucine  
        'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        # Serine
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'AGT': 'S', 'AGC': 'S',
        # Tyrosine
        'TAT': 'Y', 'TAC': 'Y',
        # Stop codons
        'TAA': '*', 'TAG': '*', 'TGA': '*',
        # Cysteine
        'TGT': 'C', 'TGC': 'C',
        # Tryptophan
        'TGG': 'W',
        # Proline
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        # Histidine
        'CAT': 'H', 'CAC': 'H',
        # Glutamine
        'CAA': 'Q', 'CAG': 'Q',
        # Arginine
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
        # Isoleucine
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
        # Methionine (start)
        'ATG': 'M',
        # Threonine
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        # Asparagine
        'AAT': 'N', 'AAC': 'N',
        # Lysine
        'AAA': 'K', 'AAG': 'K',
        # Valine
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        # Alanine
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        # Aspartic acid
        'GAT': 'D', 'GAC': 'D',
        # Glutamic acid
        'GAA': 'E', 'GAG': 'E',
        # Glycine
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

def parse_coordinates(coord_string):
    """
    Parse GenBank coordinate strings like '97..436' or '97..254,290..400'
    Returns list of (start, end) tuples (0-based for Python)
    """
    coordinates = []
    
    # Handle spliced coordinates (comma-separated)
    parts = coord_string.split(',')
    
    for part in parts:
        part = part.strip()
        if '..' in part:
            start, end = part.split('..')
            # Convert from 1-based GenBank to 0-based Python
            coordinates.append((int(start) - 1, int(end)))
        else:
            # Single position
            pos = int(part) - 1
            coordinates.append((pos, pos + 1))
    
    return coordinates

def extract_sequence_from_coordinates(dna_sequence, coordinates):
    """
    Extract DNA sequence from potentially spliced coordinates
    coordinates: list of (start, end) tuples (0-based)
    """
    extracted = ""
    
    for start, end in coordinates:
        if start < len(dna_sequence) and end <= len(dna_sequence):
            extracted += dna_sequence[start:end]
        else:
            print(f"Warning: Coordinates {start}-{end} exceed sequence length {len(dna_sequence)}")
    
    return extracted

def viral_translate(dna_sequence, coordinates=None, stop_at_stop_codon=False):
    """
    Custom translator for viral genes - no start codon detection\!
    
    Args:
        dna_sequence: DNA sequence string
        coordinates: List of (start, end) tuples for spliced genes, or None for whole sequence
        stop_at_stop_codon: If True, stop at first stop codon (default: False for mature peptides)
    
    Returns:
        Protein sequence string
    """
    codon_table = get_codon_table()
    
    # Extract sequence from coordinates if provided
    if coordinates:
        dna = extract_sequence_from_coordinates(dna_sequence, coordinates)
    else:
        dna = dna_sequence.upper()
    
    # Remove any gaps or ambiguous characters
    dna = ''.join([c for c in dna if c in 'ATCG'])
    
    # Translate 3-to-1
    protein = ""
    
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        
        if len(codon) == 3:
            amino_acid = codon_table.get(codon, 'X')  # X for unknown codons
            
            # Handle stop codons
            if amino_acid == '*':
                if stop_at_stop_codon:
                    break  # Stop translation
                else:
                    # For mature peptides, we might want to continue past internal stops
                    # This is virus-dependent behavior
                    amino_acid = '*'  # Keep the stop codon in the sequence
            
            protein += amino_acid
    
    return protein

def translate_genbank_feature(dna_sequence, feature_start, feature_end, spliced_coords=None):
    """
    Translate a GenBank feature with exact coordinates
    
    Args:
        dna_sequence: Full genome DNA sequence
        feature_start: Start position (1-based GenBank)
        feature_end: End position (1-based GenBank)  
        spliced_coords: Optional list of coordinate ranges for spliced genes
    
    Returns:
        Protein sequence
    """
    if spliced_coords:
        # Parse spliced coordinates like "97..254,290..400"
        coordinates = parse_coordinates(spliced_coords)
    else:
        # Simple coordinates
        coordinates = [(feature_start - 1, feature_end)]  # Convert to 0-based
    
    return viral_translate(dna_sequence, coordinates, stop_at_stop_codon=False)

def reverse_complement(dna_sequence):
    """Generate reverse complement for negative strand genes"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                  'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
    
    return ''.join([complement.get(c, c) for c in dna_sequence[::-1]])

# Test function
if __name__ == '__main__':
    # Test with a simple sequence
    test_dna = "TTTTTGAGGCTGAAATAA"  # Should translate to "FLR...*"
    result = viral_translate(test_dna)
    print(f"Test DNA: {test_dna}")
    print(f"Translation: {result}")
    
    # Test spliced coordinates
    test_coords = [(0, 6), (9, 15)]  # Skip middle portion
    spliced_result = viral_translate(test_dna, test_coords)
    print(f"Spliced translation: {spliced_result}")
