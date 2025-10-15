#!/usr/bin/env python3
"""
VICAST Validation Module
Provides validation and reporting functions for viral genome curation pipeline.
"""

import sys
import os
import re
import pandas as pd
from collections import defaultdict, Counter
from datetime import datetime
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation


def validate_gff_for_snpeff(gff_file, fasta_file=None):
    """
    Validate GFF3 file for snpEff compatibility.
    
    Args:
        gff_file: Path to GFF3 file
        fasta_file: Optional path to FASTA file for sequence validation
        
    Returns:
        tuple: (is_valid, errors, warnings)
            - is_valid: Boolean indicating if GFF passes critical validation
            - errors: List of error messages
            - warnings: List of warning messages
    """
    errors = []
    warnings = []
    
    # Check file exists
    if not os.path.exists(gff_file):
        errors.append(f"GFF file not found: {gff_file}")
        return False, errors, warnings
    
    # Load FASTA sequence if provided
    seq_lengths = {}
    if fasta_file and os.path.exists(fasta_file):
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                seq_lengths[record.id] = len(record.seq)
        except Exception as e:
            warnings.append(f"Could not parse FASTA file: {e}")
    
    # Parse GFF file
    features = []
    seen_ids = set()
    cds_by_gene = defaultdict(list)
    line_num = 0
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                line_num += 1
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith('#'):
                    continue
                
                # Parse GFF line
                parts = line.split('\t')
                if len(parts) != 9:
                    errors.append(f"Line {line_num}: Invalid GFF format (expected 9 columns, got {len(parts)})")
                    continue
                
                seqid, source, feat_type, start, end, score, strand, phase, attributes = parts
                
                # Validate coordinates
                try:
                    start = int(start)
                    end = int(end)
                except ValueError:
                    errors.append(f"Line {line_num}: Invalid coordinates (start={start}, end={end})")
                    continue
                
                if start > end:
                    errors.append(f"Line {line_num}: Start position ({start}) > end position ({end})")
                
                if start < 1:
                    errors.append(f"Line {line_num}: Start position ({start}) < 1")
                
                # Check against sequence length if available
                if seqid in seq_lengths:
                    if end > seq_lengths[seqid]:
                        errors.append(f"Line {line_num}: End position ({end}) exceeds sequence length ({seq_lengths[seqid]})")
                
                # Validate strand
                if strand not in ['+', '-', '.']:
                    errors.append(f"Line {line_num}: Invalid strand '{strand}' (must be +, -, or .)")
                
                # Parse attributes
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                # Check for required ID in genes and CDS
                if feat_type in ['gene', 'CDS', 'mRNA']:
                    if 'ID' not in attr_dict:
                        errors.append(f"Line {line_num}: {feat_type} feature missing required ID attribute")
                    else:
                        feat_id = attr_dict['ID']
                        if feat_id in seen_ids:
                            errors.append(f"Line {line_num}: Duplicate ID '{feat_id}'")
                        seen_ids.add(feat_id)
                
                # Track CDS features for overlap checking
                if feat_type == 'CDS':
                    gene_id = attr_dict.get('gene', attr_dict.get('Parent', 'unknown'))
                    cds_by_gene[gene_id].append((start, end, strand, line_num))
                
                # Store feature for additional checks
                features.append({
                    'line': line_num,
                    'seqid': seqid,
                    'type': feat_type,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'attributes': attr_dict
                })
                
    except Exception as e:
        errors.append(f"Error parsing GFF file: {e}")
        return False, errors, warnings
    
    # Check for overlapping CDS features within same gene
    for gene_id, cds_list in cds_by_gene.items():
        if len(cds_list) > 1:
            # Sort by start position
            cds_list.sort(key=lambda x: x[0])
            for i in range(len(cds_list) - 1):
                if cds_list[i][1] >= cds_list[i+1][0]:  # End of current >= start of next
                    warnings.append(
                        f"Overlapping CDS in gene '{gene_id}': "
                        f"Lines {cds_list[i][3]} ({cds_list[i][0]}-{cds_list[i][1]}) and "
                        f"{cds_list[i+1][3]} ({cds_list[i+1][0]}-{cds_list[i+1][1]})"
                    )
    
    # Check feature hierarchy
    gene_features = [f for f in features if f['type'] == 'gene']
    cds_features = [f for f in features if f['type'] == 'CDS']
    
    if len(gene_features) == 0 and len(cds_features) > 0:
        warnings.append("GFF contains CDS features but no gene features")
    
    # Check for features with same coordinates but different types
    coord_map = defaultdict(list)
    for f in features:
        coord_key = (f['seqid'], f['start'], f['end'])
        coord_map[coord_key].append(f['type'])
    
    for coord, types in coord_map.items():
        if len(set(types)) > 1:
            type_counts = Counter(types)
            warnings.append(
                f"Multiple feature types at position {coord[0]}:{coord[1]}-{coord[2]}: "
                f"{dict(type_counts)}"
            )
    
    # Check for essential attributes in CDS features
    for f in cds_features:
        if 'gene' not in f['attributes'] and 'Parent' not in f['attributes']:
            warnings.append(
                f"Line {f['line']}: CDS feature missing 'gene' or 'Parent' attribute"
            )
        if 'product' not in f['attributes']:
            warnings.append(
                f"Line {f['line']}: CDS feature missing 'product' attribute (recommended for snpEff)"
            )
    
    # Determine overall validity
    is_valid = len(errors) == 0
    
    return is_valid, errors, warnings


def compare_features(original_features, final_features):
    """
    Compare original and final features to identify changes.
    
    Args:
        original_features: List of features from original file
        final_features: List of features from final file
        
    Returns:
        dict: Statistics about changes
    """
    stats = {
        'original_count': len(original_features),
        'final_count': len(final_features),
        'added': [],
        'deleted': [],
        'modified': [],
        'unchanged': []
    }
    
    # Create feature maps for comparison
    orig_map = {f.get('id', f.get('key', str(i))): f for i, f in enumerate(original_features)}
    final_map = {f.get('id', f.get('key', str(i))): f for i, f in enumerate(final_features)}
    
    # Find deletions
    for feat_id, feat in orig_map.items():
        if feat_id not in final_map:
            stats['deleted'].append(feat)
    
    # Find additions
    for feat_id, feat in final_map.items():
        if feat_id not in orig_map:
            stats['added'].append(feat)
    
    # Find modifications
    for feat_id in set(orig_map.keys()) & set(final_map.keys()):
        orig_feat = orig_map[feat_id]
        final_feat = final_map[feat_id]
        
        # Compare key attributes
        if (orig_feat.get('start') != final_feat.get('start') or
            orig_feat.get('end') != final_feat.get('end') or
            orig_feat.get('gene') != final_feat.get('gene') or
            orig_feat.get('product') != final_feat.get('product')):
            stats['modified'].append({
                'id': feat_id,
                'original': orig_feat,
                'final': final_feat
            })
        else:
            stats['unchanged'].append(feat_id)
    
    return stats


def generate_curation_report(gb_file, tsv_file, gff_file, output_file):
    """
    Generate a comprehensive curation report.
    
    Args:
        gb_file: Path to original GenBank file
        tsv_file: Path to edited TSV file
        gff_file: Path to final GFF3 file
        output_file: Path to output report file
    """
    report_lines = []
    report_lines.append("="*80)
    report_lines.append("VICAST GENOME CURATION REPORT")
    report_lines.append("="*80)
    report_lines.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report_lines.append("")
    
    # Parse original GenBank
    original_features = []
    genome_info = {}
    
    if os.path.exists(gb_file):
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
                genome_info['id'] = record.id
                genome_info['description'] = record.description
                genome_info['length'] = len(record.seq)
                
                for feature in record.features:
                    if feature.type in ['CDS', 'gene', 'mat_peptide']:
                        feat_dict = {
                            'type': feature.type,
                            'start': int(feature.location.start) + 1,  # Convert to 1-based
                            'end': int(feature.location.end),
                            'strand': '+' if feature.location.strand == 1 else '-',
                            'gene': feature.qualifiers.get('gene', [''])[0],
                            'product': feature.qualifiers.get('product', [''])[0],
                            'id': feature.qualifiers.get('locus_tag', [''])[0]
                        }
                        original_features.append(feat_dict)
        except Exception as e:
            report_lines.append(f"Error parsing GenBank file: {e}")
    
    # Parse edited TSV
    edited_features = []
    if os.path.exists(tsv_file):
        try:
            df = pd.read_csv(tsv_file, sep='\t')
            for _, row in df.iterrows():
                if 'action' in df.columns and row.get('action') == 'DELETE':
                    continue
                feat_dict = {
                    'type': row.get('type', ''),
                    'start': row.get('start', 0),
                    'end': row.get('end', 0),
                    'strand': row.get('strand', '.'),
                    'gene': row.get('gene_name', '') or row.get('gene', ''),
                    'product': row.get('product', ''),
                    'id': row.get('ID', '')
                }
                edited_features.append(feat_dict)
        except Exception as e:
            report_lines.append(f"Error parsing TSV file: {e}")
    
    # Parse final GFF
    final_features = []
    if os.path.exists(gff_file):
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) == 9:
                    attrs = {}
                    for attr in parts[8].split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attrs[key] = value
                    
                    feat_dict = {
                        'type': parts[2],
                        'start': int(parts[3]),
                        'end': int(parts[4]),
                        'strand': parts[6],
                        'gene': attrs.get('gene', ''),
                        'product': attrs.get('product', ''),
                        'id': attrs.get('ID', '')
                    }
                    final_features.append(feat_dict)
    
    # Generate genome information section
    if genome_info:
        report_lines.append("GENOME INFORMATION")
        report_lines.append("-"*40)
        report_lines.append(f"Genome ID: {genome_info.get('id', 'N/A')}")
        report_lines.append(f"Description: {genome_info.get('description', 'N/A')}")
        report_lines.append(f"Length: {genome_info.get('length', 0):,} bp")
        report_lines.append("")
    
    # Compare features
    if original_features and final_features:
        stats = compare_features(original_features, final_features)
        
        report_lines.append("CURATION SUMMARY")
        report_lines.append("-"*40)
        report_lines.append(f"Original features: {stats['original_count']}")
        report_lines.append(f"Final features: {stats['final_count']}")
        report_lines.append(f"Net change: {stats['final_count'] - stats['original_count']:+d}")
        report_lines.append("")
        report_lines.append(f"Features added: {len(stats['added'])}")
        report_lines.append(f"Features deleted: {len(stats['deleted'])}")
        report_lines.append(f"Features modified: {len(stats['modified'])}")
        report_lines.append(f"Features unchanged: {len(stats['unchanged'])}")
        report_lines.append("")
        
        # Detail changes
        if stats['deleted']:
            report_lines.append("DELETED FEATURES")
            report_lines.append("-"*40)
            for feat in stats['deleted'][:20]:  # Show first 20
                report_lines.append(
                    f"  - {feat['type']} at {feat['start']}-{feat['end']}: "
                    f"{feat.get('gene', 'unknown')} ({feat.get('product', 'no product')})"
                )
            if len(stats['deleted']) > 20:
                report_lines.append(f"  ... and {len(stats['deleted']) - 20} more")
            report_lines.append("")
        
        if stats['added']:
            report_lines.append("ADDED FEATURES")
            report_lines.append("-"*40)
            for feat in stats['added'][:20]:  # Show first 20
                report_lines.append(
                    f"  + {feat['type']} at {feat['start']}-{feat['end']}: "
                    f"{feat.get('gene', 'unknown')} ({feat.get('product', 'no product')})"
                )
            if len(stats['added']) > 20:
                report_lines.append(f"  ... and {len(stats['added']) - 20} more")
            report_lines.append("")
        
        if stats['modified']:
            report_lines.append("MODIFIED FEATURES")
            report_lines.append("-"*40)
            for item in stats['modified'][:20]:  # Show first 20
                orig = item['original']
                final = item['final']
                changes = []
                if orig['start'] != final['start'] or orig['end'] != final['end']:
                    changes.append(f"coordinates: {orig['start']}-{orig['end']} -> {final['start']}-{final['end']}")
                if orig.get('gene') != final.get('gene'):
                    changes.append(f"gene: '{orig.get('gene')}' -> '{final.get('gene')}'")
                if orig.get('product') != final.get('product'):
                    changes.append(f"product changed")
                
                report_lines.append(f"  * {item['id']}: {', '.join(changes)}")
            if len(stats['modified']) > 20:
                report_lines.append(f"  ... and {len(stats['modified']) - 20} more")
            report_lines.append("")
    
    # Feature type summary
    if final_features:
        type_counts = Counter(f['type'] for f in final_features)
        report_lines.append("FINAL FEATURE TYPES")
        report_lines.append("-"*40)
        for feat_type, count in sorted(type_counts.items()):
            report_lines.append(f"  {feat_type}: {count}")
        report_lines.append("")
    
    # Validation results
    report_lines.append("VALIDATION RESULTS")
    report_lines.append("-"*40)
    
    if os.path.exists(gff_file):
        is_valid, errors, warnings = validate_gff_for_snpeff(gff_file)
        
        if is_valid:
            report_lines.append("[SUCCESS] GFF validation passed")
        else:
            report_lines.append("[ERROR] GFF validation failed")
        
        if errors:
            report_lines.append(f"\nErrors ({len(errors)}):")
            for error in errors[:10]:
                report_lines.append(f"  - {error}")
            if len(errors) > 10:
                report_lines.append(f"  ... and {len(errors) - 10} more")
        
        if warnings:
            report_lines.append(f"\nWarnings ({len(warnings)}):")
            for warning in warnings[:10]:
                report_lines.append(f"  - {warning}")
            if len(warnings) > 10:
                report_lines.append(f"  ... and {len(warnings) - 10} more")
    else:
        report_lines.append("GFF file not found for validation")
    
    report_lines.append("")
    report_lines.append("="*80)
    report_lines.append("END OF REPORT")
    report_lines.append("="*80)
    
    # Write report to file
    with open(output_file, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Curation report generated: {output_file}")
    
    # Also print summary to console
    print("\nCuration Summary:")
    if original_features and final_features:
        print(f"  Original features: {len(original_features)}")
        print(f"  Final features: {len(final_features)}")
        print(f"  Net change: {len(final_features) - len(original_features):+d}")


def validate_tsv_for_conversion(tsv_file):
    """
    Quick validation of TSV file before conversion to GFF.
    
    Args:
        tsv_file: Path to TSV file
        
    Returns:
        tuple: (is_valid, errors)
    """
    errors = []
    
    try:
        df = pd.read_csv(tsv_file, sep='\t')
    except Exception as e:
        errors.append(f"Cannot read TSV file: {e}")
        return False, errors
    
    # Check required columns
    required_cols = ['seqid', 'source', 'type', 'start', 'end', 'strand']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        errors.append(f"Missing required columns: {missing}")
    
    # Check for valid coordinates
    if 'start' in df.columns and 'end' in df.columns:
        invalid_coords = df[(df['start'] > df['end']) | (df['start'] < 1)]
        if not invalid_coords.empty:
            errors.append(f"{len(invalid_coords)} rows with invalid coordinates")
    
    # Check for valid strand values
    if 'strand' in df.columns:
        invalid_strand = df[~df['strand'].isin(['+', '-', '.'])]
        if not invalid_strand.empty:
            errors.append(f"{len(invalid_strand)} rows with invalid strand values")
    
    is_valid = len(errors) == 0
    return is_valid, errors


if __name__ == '__main__':
    # Quick test if run directly
    import sys
    if len(sys.argv) > 1:
        gff_file = sys.argv[1]
        fasta_file = sys.argv[2] if len(sys.argv) > 2 else None
        
        print("Validating GFF file...")
        is_valid, errors, warnings = validate_gff_for_snpeff(gff_file, fasta_file)
        
        print(f"\nValidation result: {'PASSED' if is_valid else 'FAILED'}")
        if errors:
            print(f"\nErrors ({len(errors)}):")
            for e in errors[:10]:
                print(f"  - {e}")
        if warnings:
            print(f"\nWarnings ({len(warnings)}):")
            for w in warnings[:10]:
                print(f"  - {w}")
    else:
        print("Usage: python3 vicast_validation.py <gff_file> [fasta_file]")