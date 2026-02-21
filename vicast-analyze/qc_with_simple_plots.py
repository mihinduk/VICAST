#!/usr/bin/env python3
"""
QC Report with Simple Plots (minimal dependencies)
Creates text report + basic ASCII plots for presentations
"""

import os
import sys
import argparse
from pathlib import Path
import urllib.request
import urllib.parse
import json
import re

def lookup_virus_name(accession):
    """Look up virus name from NCBI accession number"""
    if not accession:
        return None
        
    try:
        # Use NCBI E-utilities to get sequence information
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            'db': 'nucleotide',
            'id': accession,
            'rettype': 'fasta',
            'retmode': 'text'
        }
        
        url = f"{base_url}?{urllib.parse.urlencode(params)}"
        
        with urllib.request.urlopen(url, timeout=10) as response:
            fasta_content = response.read().decode('utf-8')
            
        # Extract organism name from FASTA header
        if fasta_content.startswith('>'):
            header = fasta_content.split('\n')[0]
            
            # Look for virus name patterns
            virus_patterns = [
                r'(\w+\s+virus)',  # "Zika virus", "West Nile virus"
                r'(\w+\s+\w+\s+virus)',  # "Chikungunya virus", "Dengue virus"
                r'(\w+virus)',  # "Poliovirus"
            ]
            
            for pattern in virus_patterns:
                match = re.search(pattern, header, re.IGNORECASE)
                if match:
                    return match.group(1).lower()
                    
    except Exception as e:
        print(f"Warning: Could not lookup accession {accession}: {e}")
        
    return None

def parse_diagnostic_report(report_file):
    """Parse viral diagnostic report to extract key metrics"""
    data = {
        'sample': '',
        'total_reads': 0,
        'duplicate_reads': 0,
        'unique_reads': 0,
        'raw_mapping_reads': 0,
        'raw_mapping_percent': 0,
        'dedup_mapping_reads': 0,
        'dedup_mapping_percent': 0,
        'duplication_rate': 0,
        'total_contigs': 0,
        'contigs_gt1000': 0
    }
    
    try:
        with open(report_file, 'r') as f:
            content = f.read()
            
        # Extract sample name
        for line in content.split('\n'):
            if line.startswith('Sample:'):
                data['sample'] = line.split(':', 1)[1].strip()
                break
                
        # Extract mapping statistics
        lines = content.split('\n')
        for i, line in enumerate(lines):
            try:
                if 'Total Reads:' in line:
                    value = line.split(':')[1].strip().replace(',', '')
                    data['total_reads'] = int(value) if value else 0
                elif 'Duplicate Reads:' in line:
                    parts = line.split(':')[1].strip().split()
                    if parts:
                        data['duplicate_reads'] = int(parts[0].replace(',', ''))
                    if '(' in line:
                        data['duplication_rate'] = float(line.split('(')[1].split('%')[0])
                elif 'Unique Reads:' in line:
                    value = line.split(':')[1].strip().replace(',', '')
                    data['unique_reads'] = int(value) if value else 0
                elif 'Raw Mapping:' in line:
                    parts = line.split(':')[1].strip().split()
                    if parts and parts[0]:
                        data['raw_mapping_reads'] = int(parts[0].replace(',', ''))
                    if '(' in line:
                        data['raw_mapping_percent'] = float(line.split('(')[1].split('%')[0])
                elif 'Deduplicated Mapping:' in line:
                    parts = line.split(':')[1].strip().split()
                    if parts and parts[0]:
                        data['dedup_mapping_reads'] = int(parts[0].replace(',', ''))
                    if '(' in line:
                        data['dedup_mapping_percent'] = float(line.split('(')[1].split('%')[0])
                elif 'Total Contigs:' in line:
                    value = line.split(':')[1].strip()
                    data['total_contigs'] = int(value) if value else 0
                elif 'Contigs >1000bp:' in line:
                    value = line.split(':')[1].strip()
                    data['contigs_gt1000'] = int(value) if value else 0
            except (ValueError, IndexError):
                # Skip lines that can't be parsed
                continue
                
    except Exception as e:
        print(f"Warning: Could not parse {report_file}: {e}")
        
    return data

def resolve_segment_accessions(accession):
    """Resolve a custom accession to its segment accessions from the manifest.

    For segmented viruses (e.g., influenza_pr8), the custom accession won't
    appear in BLAST results. This returns a set including individual segment
    accessions so they can be matched.
    """
    accessions = {accession.lower()} if accession else set()
    if not accession:
        return accessions

    script_dir = Path(__file__).parent
    manifest_paths = [
        script_dir.parent / "prebuilt_databases" / "manifest.json",
        script_dir / "manifest.json",
        Path("/opt/vicast/prebuilt_databases/manifest.json"),
    ]

    for mpath in manifest_paths:
        if mpath.exists():
            try:
                with open(mpath, 'r') as f:
                    manifest = json.load(f)
                for entry in manifest.get("databases", []):
                    if entry.get("accession") == accession:
                        segments = entry.get("segment_accessions", [])
                        if segments:
                            accessions.update(s.lower() for s in segments)
                        break
            except (json.JSONDecodeError, KeyError):
                pass
            break
    return accessions


def parse_blast_results(blast_file, target_accession=None, target_virus_name=None):
    """Parse BLAST top_hits.tsv to extract contamination and target genome information.

    Header-aware: looks up columns by name to handle format changes.
    top_hits.tsv contains only best hits (tied-best when multiple have identical stats).
    Deduplicates by contig_id: 1 entry per contig unless multiple tied-best hits exist.
    """
    # Resolve segment accessions for segmented viruses
    target_accession_set = resolve_segment_accessions(target_accession)
    contamination_data = []
    target_data = []

    try:
        if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
            with open(blast_file, 'r') as f:
                lines = f.readlines()

            if len(lines) > 1:
                # Parse header to find column indices by name
                header = lines[0].strip().split('\t')
                col = {name: i for i, name in enumerate(header)}

                # Required columns (with fallback indices for old formats)
                idx_contig = col.get('Contig_ID', 0)
                idx_length = col.get('Contig_Length', 1)
                idx_subject = col.get('Subject_ID', 2)
                idx_pident = col.get('Percent_Identity', 3)
                idx_qcov = col.get('Query_Coverage', None)
                idx_kingdom = col.get('Kingdom/Type', -2)
                idx_title = col.get('Subject_Title', -1)

                for line in lines[1:]:
                    parts = line.strip().split('\t')
                    if len(parts) < 4:
                        continue

                    # Parse percent identity (remove % symbol)
                    try:
                        identity = float(parts[idx_pident].replace('%', ''))
                    except (ValueError, IndexError):
                        identity = 0

                    # Parse query coverage (remove % symbol)
                    coverage = 0
                    if idx_qcov is not None and idx_qcov < len(parts):
                        try:
                            coverage = float(parts[idx_qcov].replace('%', ''))
                        except (ValueError, IndexError):
                            coverage = 0

                    contig_id = parts[idx_contig]
                    contig_len = parts[idx_length] if idx_length < len(parts) else 'Unknown'
                    subject_id = parts[idx_subject] if idx_subject < len(parts) else ''
                    title = parts[idx_title].strip('"') if abs(idx_title) <= len(parts) else ''
                    kingdom = parts[idx_kingdom].strip('"') if abs(idx_kingdom) <= len(parts) else ''
                    title_lower = title.lower()

                    # Check if this matches target genome
                    is_target = False
                    sid_lower = subject_id.lower()
                    if target_accession_set and any(
                            acc in sid_lower or acc in title_lower
                            for acc in target_accession_set):
                        is_target = True
                    elif target_virus_name and target_virus_name in title_lower:
                        is_target = True

                    if is_target:
                        target_data.append({
                            'contig_id': contig_id,
                            'accession': subject_id,
                            'organism': title[:80] if title else 'Unknown',
                            'identity': identity,
                            'coverage': coverage,
                            'length': contig_len,
                        })

                    if not is_target:
                        # Categorize organisms using kingdom column or title
                        if kingdom == 'Virus' or 'virus' in title_lower:
                            category = 'Virus'
                        elif kingdom == 'Mycoplasma' or any(x in title_lower for x in ['mycoplasma', 'mesomycoplasma']):
                            category = 'Mycoplasma'
                        elif any(x in title_lower for x in ['escherichia', 'e. coli']):
                            category = 'E. coli'
                        elif any(x in title_lower for x in ['staphylococcus', 'staph']):
                            category = 'Staphylococcus'
                        elif any(x in title_lower for x in ['pseudomonas']):
                            category = 'Pseudomonas'
                        elif any(x in title_lower for x in ['candida']):
                            category = 'Candida'
                        elif any(x in title_lower for x in ['saccharomyces']):
                            category = 'Yeast'
                        elif any(x in title_lower for x in ['cryptococcus']):
                            category = 'Cryptococcus'
                        else:
                            category = 'Other'

                        contamination_data.append({
                            'contig_id': contig_id,
                            'accession': subject_id,
                            'category': category,
                            'organism': title[:60] + '...' if len(title) > 60 else title if title else 'Unknown',
                            'identity': identity,
                            'coverage': coverage,
                            'length': contig_len,
                        })

        # Deduplicate target_data by contig_id (keep first/best entry per contig)
        seen_target_contigs = set()
        deduped_target = []
        for entry in target_data:
            if entry['contig_id'] not in seen_target_contigs:
                seen_target_contigs.add(entry['contig_id'])
                deduped_target.append(entry)
        target_data = deduped_target

        # Deduplicate contamination_data by contig_id (keep first/best entry per contig)
        seen_contam_contigs = set()
        deduped_contam = []
        for entry in contamination_data:
            if entry['contig_id'] not in seen_contam_contigs:
                seen_contam_contigs.add(entry['contig_id'])
                deduped_contam.append(entry)
        contamination_data = deduped_contam

    except Exception as e:
        print(f"Warning: Could not parse {blast_file}: {e}")

    return contamination_data, target_data

def create_simple_html_charts(samples_data, all_contamination, all_target_data, target_accession, target_virus_name, output_dir):
    """Create simple HTML file with basic charts using only CSS"""
    
    # Create unique filename based on output directory
    output_basename = os.path.basename(output_dir.rstrip('/'))
    html_file = os.path.join(output_dir, f'{output_basename}_presentation_ready_report.html')
    
    # Calculate contamination counts
    contamination_counts = {}
    for item in all_contamination:
        cat = item['category']
        contamination_counts[cat] = contamination_counts.get(cat, 0) + 1
    
    total_contamination = sum(contamination_counts.values())
    
    with open(html_file, 'w') as f:
        f.write("""<!DOCTYPE html>
<html>
<head>
    <title>Viral Culture QC Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
        .container { max-width: 1000px; margin: 0 auto; background: white; padding: 30px; border-radius: 10px; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
        h1 { color: #2c3e50; text-align: center; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; border-left: 4px solid #3498db; padding-left: 15px; }
        .metric-card { background: #ecf0f1; padding: 15px; margin: 10px 0; border-radius: 8px; border-left: 5px solid #3498db; }
        .excellent { border-left-color: #27ae60 !important; }
        .moderate { border-left-color: #f39c12 !important; }
        .poor { border-left-color: #e74c3c !important; }
        .bar-chart { margin: 20px 0; }
        .bar { height: 30px; margin: 5px 0; position: relative; background: #ecf0f1; border-radius: 5px; overflow: hidden; }
        .bar-fill { height: 100%; border-radius: 5px; position: relative; }
        .bar-label { position: absolute; left: 10px; top: 50%; transform: translateY(-50%); font-weight: bold; z-index: 2; color: white; }
        .bar-value { position: absolute; right: 10px; top: 50%; transform: translateY(-50%); font-weight: bold; z-index: 2; color: white; }
        .excellent-bar { background: linear-gradient(90deg, #27ae60, #2ecc71); }
        .moderate-bar { background: linear-gradient(90deg, #f39c12, #e67e22); }
        .poor-bar { background: linear-gradient(90deg, #e74c3c, #c0392b); }
        .contamination-item { padding: 10px; margin: 5px 0; background: #fff; border-radius: 5px; border-left: 4px solid #e74c3c; }
        .alert { padding: 15px; margin: 15px 0; border-radius: 8px; font-weight: bold; }
        .alert-critical { background: #f8d7da; color: #721c24; border: 1px solid #f5c6cb; }
        .alert-warning { background: #fff3cd; color: #856404; border: 1px solid #ffeaa7; }
        .alert-success { background: #d4edda; color: #155724; border: 1px solid #c3e6cb; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px; margin: 20px 0; }
        .stat-box { text-align: center; padding: 20px; background: #3498db; color: white; border-radius: 8px; }
        .stat-number { font-size: 2em; font-weight: bold; display: block; }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Viral Culture Quality Control Report</h1>
        
        <div class="stats-grid">""")
        
        # Summary statistics
        total_samples = len(samples_data)
        total_contigs = sum(s['contigs_gt1000'] for s in samples_data)
        avg_mapping = sum(s['dedup_mapping_percent'] for s in samples_data) / total_samples if total_samples > 0 else 0
        
        f.write(f"""
            <div class="stat-box">
                <span class="stat-number">{total_samples}</span>
                Samples Analyzed
            </div>
            <div class="stat-box">
                <span class="stat-number">{total_contigs}</span>
                Contigs >1kb Processed
            </div>
            <div class="stat-box">
                <span class="stat-number">{avg_mapping:.1f}%</span>
                Average Mapping
            </div>
            <div class="stat-box">
                <span class="stat-number">{total_contamination}</span>
                Contamination Contigs
            </div>
        </div>
        
        <h2>📊 Sample Quality Overview</h2>
        <div class="bar-chart">""")
        
        # Sample quality bars
        for sample in samples_data:
            mapping_pct = sample['dedup_mapping_percent']
            
            if mapping_pct >= 70:
                bar_class = "excellent-bar"
                card_class = "excellent"
                status = "🟢 EXCELLENT"
            elif mapping_pct >= 30:
                bar_class = "moderate-bar" 
                card_class = "moderate"
                status = "🟡 MODERATE"
            else:
                bar_class = "poor-bar"
                card_class = "poor"
                status = "🔴 POOR"
            
            f.write(f"""
            <div class="metric-card {card_class}">
                <strong>{sample['sample']}</strong> - {status}
                <div class="bar">
                    <div class="bar-fill {bar_class}" style="width: {mapping_pct}%">
                        <span class="bar-label">Mapping: {mapping_pct:.1f}%</span>
                    </div>
                </div>
                <div style="margin-top: 5px; font-size: 0.9em; color: #555;">{sample['contigs_gt1000']} contigs analyzed</div>
            </div>""")
        
        f.write("""
        </div>
        
        <h2>🦠 Contamination Analysis</h2>""")
        
        # Target genome analysis
        target_contigs = len(all_target_data)
        target_name = "Unknown"
        if all_target_data:
            target_name = all_target_data[0]['organism'].split(',')[0]  # Get main organism name
        elif target_virus_name:
            target_name = target_virus_name.title()  # Capitalize properly
        elif target_accession:
            target_name = f"Target genome ({target_accession})"
            
        # Contamination alerts - ANY organism that is not target genome is contamination
        mycoplasma_found = any(item['category'] == 'Mycoplasma' for item in all_contamination)
        bacterial_found = any(item['category'] in ['E. coli', 'Staphylococcus', 'Pseudomonas'] for item in all_contamination)
        viral_contamination = any(item['category'] == 'Virus' for item in all_contamination)
        any_contamination = len(all_contamination) > 0
        
        if mycoplasma_found:
            f.write('<div class="alert alert-critical">🔴 CRITICAL: Mycoplasma contamination detected! Immediate culture cleaning required.</div>')
        if viral_contamination:
            f.write(f'<div class="alert alert-critical">🔴 CRITICAL: Viral contamination detected! Sample is NOT pure {target_name}.</div>')
        if bacterial_found:
            f.write('<div class="alert alert-warning">🟡 WARNING: Bacterial contamination detected. Culture cleaning recommended.</div>')
        if any_contamination and not mycoplasma_found and not viral_contamination and not bacterial_found:
            f.write('<div class="alert alert-warning">🟡 WARNING: Other contamination detected. Review BLAST results.</div>')
        if not any_contamination:
            f.write('<div class="alert alert-success">✅ SUCCESS: No contamination detected in BLAST analysis.</div>')
        else:
            f.write(f'<div class="alert alert-critical">🔴 CONTAMINATION DETECTED: Sample contains non-target organisms and is NOT pure {target_name}.</div>')
            
        # Target genome information
        if target_contigs > 0:
            f.write(f'<div class="alert alert-success">✅ TARGET GENOME: {target_contigs} contigs match {target_name} ({target_accession})</div>')
            
            # Detailed target genome table
            f.write('<h3>🎯 Target Genome Contigs</h3>')
            f.write('<table style="width: 100%; border-collapse: collapse; margin: 20px 0; border: 1px solid #ddd;">')
            f.write('<tr style="background: #27ae60; color: white;"><th style="padding: 12px; text-align: left; border: 1px solid #ddd;">Contig ID</th><th style="padding: 12px; border: 1px solid #ddd;">Accession</th><th style="padding: 12px; border: 1px solid #ddd;">Identity %</th><th style="padding: 12px; border: 1px solid #ddd;">Length (bp)</th><th style="padding: 12px; border: 1px solid #ddd;">% Coverage</th><th style="padding: 12px; border: 1px solid #ddd;">Organism</th></tr>')

            for i, target in enumerate(all_target_data):
                bg_color = "#d5f4e6" if i % 2 == 0 else "#ffffff"
                length_display = f"{target['length']:,}" if isinstance(target['length'], (int, float)) else target['length']
                accession_display = target.get('accession', '')
                coverage_display = f"{target.get('coverage', 0):.1f}%"
                f.write(f'<tr style="background: {bg_color};"><td style="padding: 12px; border: 1px solid #ddd;">{target["contig_id"]}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{accession_display}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{target["identity"]:.1f}%</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{length_display}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{coverage_display}</td><td style="padding: 12px; border: 1px solid #ddd;">{target["organism"][:80] + "..." if len(target["organism"]) > 80 else target["organism"]}</td></tr>')
            
            f.write('</table>')
        else:
            f.write(f'<div class="alert alert-critical">🔴 NO TARGET CONTIGS: No contigs match expected {target_name} ({target_accession})</div>')
        
        # Contamination breakdown
        if contamination_counts:
            # Summary bar chart for all categories
            f.write('<h3>Contamination Summary</h3>')
            f.write('<div class="bar-chart">')
            for category, count in sorted(contamination_counts.items()):
                percentage = (count / total_contamination) * 100
                f.write(f"""
                <div class="contamination-item">
                    <strong>{category}</strong>: {count} contigs ({percentage:.1f}%)
                    <div class="bar">
                        <div class="bar-fill poor-bar" style="width: {percentage}%"></div>
                    </div>
                </div>""")
            f.write('</div>')

            # Detailed table for viral contaminants
            viral_contam = [item for item in all_contamination if item['category'] == 'Virus']
            if viral_contam:
                f.write('<h3>Non-Target Viral Contigs</h3>')
                f.write('<table style="width: 100%; border-collapse: collapse; margin: 20px 0; border: 1px solid #ddd;">')
                f.write('<tr style="background: #c0392b; color: white;"><th style="padding: 12px; text-align: left; border: 1px solid #ddd;">Contig ID</th><th style="padding: 12px; border: 1px solid #ddd;">Accession</th><th style="padding: 12px; border: 1px solid #ddd;">Identity %</th><th style="padding: 12px; border: 1px solid #ddd;">Length (bp)</th><th style="padding: 12px; border: 1px solid #ddd;">% Coverage</th><th style="padding: 12px; border: 1px solid #ddd;">Organism</th></tr>')
                for i, item in enumerate(viral_contam):
                    bg_color = "#fde8e8" if i % 2 == 0 else "#ffffff"
                    length_display = f"{item['length']:,}" if isinstance(item['length'], (int, float)) else item['length']
                    f.write(f'<tr style="background: {bg_color};"><td style="padding: 12px; border: 1px solid #ddd;">{item["contig_id"]}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{item.get("accession", "")}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{item["identity"]:.1f}%</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{length_display}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{item.get("coverage", 0):.1f}%</td><td style="padding: 12px; border: 1px solid #ddd;">{item["organism"]}</td></tr>')
                f.write('</table>')
        
        f.write("""
        <h2>📋 Detailed Sample Information</h2>""")
        
        # Detailed sample table
        f.write('<table style="width: 100%; border-collapse: collapse; margin: 20px 0; border: 1px solid #ddd;">')
        f.write('<tr style="background: #34495e; color: white;"><th style="padding: 12px; text-align: left; border: 1px solid #ddd;">Sample</th><th style="padding: 12px; border: 1px solid #ddd;">Mapping %</th><th style="padding: 12px; border: 1px solid #ddd;">Duplication %</th><th style="padding: 12px; border: 1px solid #ddd;">Total Reads</th><th style="padding: 12px; border: 1px solid #ddd;">Contigs >1kb</th></tr>')
        
        for i, sample in enumerate(samples_data):
            bg_color = "#f8f9fa" if i % 2 == 0 else "#ffffff"
            f.write(f'<tr style="background: {bg_color};"><td style="padding: 12px; border: 1px solid #ddd;">{sample["sample"]}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{sample["dedup_mapping_percent"]:.1f}%</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{sample["duplication_rate"]:.1f}%</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{sample["total_reads"]:,}</td><td style="text-align: center; padding: 12px; border: 1px solid #ddd;">{sample["contigs_gt1000"]}</td></tr>')
        
        f.write("""</table>
        
        <div style="text-align: center; margin-top: 40px; padding: 20px; background: #ecf0f1; border-radius: 8px;">
            <strong>Report Generated:</strong> """ + __import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S') + """<br>
            <em>Ready for presentation! 🚀</em>
        </div>
        
    </div>
</body>
</html>""")

def extract_target_accession(diagnostic_dir):
    """Extract target accession from diagnostic report"""
    sample_name = os.path.basename(diagnostic_dir)
    if sample_name.startswith('diagnostic_'):
        sample_name = sample_name[11:]
    
    report_file = os.path.join(diagnostic_dir, f"{sample_name}_diagnostic_report.txt")
    try:
        with open(report_file, 'r') as f:
            content = f.read()
        for line in content.split('\n'):
            if 'Reference:' in line:
                return line.split(':', 1)[1].strip()
    except:
        pass
    return None

def main():
    parser = argparse.ArgumentParser(description='Create QC report with simple visualizations')
    parser.add_argument('diagnostic_dirs', nargs='+', help='Diagnostic output directories')
    parser.add_argument('-o', '--output', help='Output directory for report (defaults to first diagnostic_dir)')
    
    args = parser.parse_args()
    
    # Default output to first diagnostic directory if not specified
    if not args.output:
        args.output = args.diagnostic_dirs[0]
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Extract target accession from first diagnostic directory
    target_accession = extract_target_accession(args.diagnostic_dirs[0])
    
    # Look up the virus name for the target accession
    target_virus_name = lookup_virus_name(target_accession)
    if target_virus_name:
        print(f"✅ Target virus identified: {target_virus_name.title()} ({target_accession})")
    else:
        print(f"⚠️  Could not identify virus type for {target_accession}, using accession-only matching")
    
    # Process all samples
    samples_data = []
    all_contamination = []
    all_target_data = []
    
    for diagnostic_dir in args.diagnostic_dirs:
        if not os.path.exists(diagnostic_dir):
            print(f"Warning: Directory not found: {diagnostic_dir}")
            continue
            
        # Find sample name
        sample_name = os.path.basename(diagnostic_dir)
        if sample_name.startswith('diagnostic_'):
            sample_name = sample_name[11:]  # Remove 'diagnostic_' prefix
            
        # Parse diagnostic report
        report_file = os.path.join(diagnostic_dir, f"{sample_name}_diagnostic_report.txt")
        sample_data = parse_diagnostic_report(report_file)
        
        if sample_data['sample']:
            samples_data.append(sample_data)
            
            # Parse contamination and target data
            blast_file = os.path.join(diagnostic_dir, f"{sample_name}_top_hits.tsv")
            contamination_data, target_data = parse_blast_results(blast_file, target_accession, target_virus_name)
            
            # Add sample name to contamination data
            for item in contamination_data:
                item['sample'] = sample_data['sample']
            all_contamination.extend(contamination_data)
            
            # Add sample name to target data
            for item in target_data:
                item['sample'] = sample_data['sample']
            all_target_data.extend(target_data)
            
            print(f"✅ Processed sample: {sample_data['sample']}")
    
    # Create reports
    if samples_data:
        create_simple_html_charts(samples_data, all_contamination, all_target_data, target_accession, target_virus_name, args.output)
        
        output_basename = os.path.basename(args.output.rstrip('/'))
        print(f"\n🎉 QC report saved to: {args.output}/")
        print(f"🌐 {output_basename}_presentation_ready_report.html - HTML report with visual charts")
        print("\nReady for your presentation tomorrow! 🚀")
    else:
        print("❌ No valid diagnostic data found")

if __name__ == "__main__":
    main()