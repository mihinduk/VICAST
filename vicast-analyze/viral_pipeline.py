#!/usr/bin/env python3
"""
Shotgun Viral Genomics Pipeline

A comprehensive pipeline for processing shotgun sequencing data from viral samples,
performing quality control, mapping to reference genomes, variant calling, and annotation.

Copyright (c) 2024
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Dict, Union, Any, Tuple
import glob  # Added for glob file pattern support
from Bio import Entrez  # For NCBI genome downloads

# Set email for NCBI Entrez
Entrez.email = "vicast@example.com"

__version__ = "0.1.0"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

# Progress tracking helpers
def print_step_header(step_num: int, total_steps: int, step_name: str) -> None:
    """Print a formatted step header with progress indicator."""
    separator = "=" * 70
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logger.info(separator)
    logger.info(f"⏳ STEP {step_num}/{total_steps}: {step_name}")
    logger.info(f"Started at: {timestamp}")
    logger.info(separator)

def print_step_complete(step_name: str, output_files: Optional[List[str]] = None) -> None:
    """Print a completion message with output files if provided."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    logger.info(f"✓ COMPLETE: {step_name}")
    logger.info(f"Finished at: {timestamp}")
    if output_files:
        logger.info("Output files:")
        for f in output_files:
            if os.path.exists(f):
                logger.info(f"  → {f}")
    logger.info("")

def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Shotgun Viral Genomics Pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    
    # Input options
    parser.add_argument("--r1", required=True, help="Forward reads (R1) FASTQ files (accepts wildcards like *_R1.fastq.gz or *_R1_001.fastq.gz)")
    parser.add_argument("--r2", required=True, help="Reverse reads (R2) FASTQ files (accepts wildcards like *_R2.fastq.gz or *_R2_001.fastq.gz)")
    
    # Reference genome options
    ref_group = parser.add_argument_group("Reference Genome")
    ref_group.add_argument("--accession", help="Reference genome accession number (e.g., NC_045512.2)")
    ref_group.add_argument("--reference", help="Path to local reference genome (when not using accession)")
    ref_group.add_argument("--skip-download", action="store_true", help="Skip genome download (use local reference)")
    ref_group.add_argument("--force-download", action="store_true", help="Force download even if reference exists locally")
    
    # Pipeline control
    pipeline_group = parser.add_argument_group("Pipeline Control")
    pipeline_group.add_argument("--skip-qc", action="store_true", help="Skip quality control step")
    pipeline_group.add_argument("--skip-stats", action="store_true", help="Skip read statistics step")
    pipeline_group.add_argument("--skip-mapping", action="store_true", help="Skip read mapping step")
    pipeline_group.add_argument("--skip-variants", action="store_true", help="Skip variant calling step")
    pipeline_group.add_argument("--skip-annotation", action="store_true", help="Skip variant annotation step")
    pipeline_group.add_argument("--skip-vector-filter", action="store_true", help="Skip cloning vector read filtering (default: filter using UniVec)")
    pipeline_group.add_argument("--resume-from-vcf", action="store_true", help="Resume from existing VCF files (skip Steps 1-6, run only annotation Steps 7-9)")
    pipeline_group.add_argument("--create-genbank", action="store_true", help="Create GenBank file from FASTA using BLAST annotation")
    
    # Output options
    output_group = parser.add_argument_group("Output")
    output_group.add_argument("--outdir", default=".", help="Output directory")
    output_group.add_argument("--min-depth", type=int, default=200, help="Minimum read depth for lofreq filter (default: 200)")
    output_group.add_argument("--min-qual", type=int, default=1000, help="Minimum variant quality (phred) for lofreq filter (default: 1000)")
    output_group.add_argument("--keep-tmp", action="store_true", help="Keep temporary files")
    
    # Performance options
    perf_group = parser.add_argument_group("Performance")
    perf_group.add_argument("--threads", type=int, default=1, help="Number of CPU threads to use")
    perf_group.add_argument("--large-files", action="store_true", help="Enable high-memory mode for large files (1-5GB). Increases Java heap to 32GB and uses 8GB per thread for sorting.")
    perf_group.add_argument("--extremely-large-files", action="store_true", help="Enable extreme high-memory mode for massive files (>5GB). Increases Java heap to 64GB and uses 32GB per thread for sorting. Requires high-memory compute nodes.")
    
    # SnpEff options
    snpeff_group = parser.add_argument_group("SnpEff")
    snpeff_group.add_argument("--add-to-snpeff", action="store_true", help="Attempt to add genome to snpEff if not present")
    snpeff_group.add_argument("--snpeff-jar", default="snpEff.jar", help="Path to snpEff.jar")
    snpeff_group.add_argument("--java-path", default="java", help="Path to Java executable")
    
    # BLAST options for GenBank creation
    blast_group = parser.add_argument_group("BLAST Options (for GenBank creation)")
    blast_group.add_argument("--blast-db", default="nt", help="BLAST database to use for GenBank creation")
    blast_group.add_argument("--email", help="Email for NCBI Entrez queries (required for GenBank creation)")
    
    return parser.parse_args()

def download_reference_genome(accession: str, output_dir: str, force: bool = False) -> str:
    """
    Get a reference genome: try NCBI download first, then fall back to
    SnpEff database sequences.fa for custom genomes (e.g. segmented viruses).

    Args:
        accession: The genome accession number (e.g., NC_045512.2) or custom name
        output_dir: Directory to save the downloaded genome
        force: Force download even if the file exists

    Returns:
        Path to the downloaded genome FASTA file
    """
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"{accession}.fasta")

    # Check if file already exists
    if os.path.exists(output_file) and not force:
        logger.info(f"Reference genome file already exists: {output_file}")
        return output_file

    logger.info(f"Downloading reference genome: {accession}")

    # Use BioPython Entrez (more reliable than command-line tools)
    ncbi_success = False
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession,
                               rettype="fasta", retmode="text")
        with open(output_file, 'w') as f:
            f.write(handle.read())
        handle.close()
        # Verify it's a real FASTA (not an error page)
        with open(output_file) as f:
            first_line = f.readline().strip()
        if first_line.startswith(">"):
            ncbi_success = True
            logger.info(f"Reference genome downloaded successfully via BioPython")
        else:
            logger.warning(f"NCBI returned non-FASTA content for {accession}")
    except Exception as e:
        logger.warning(f"BioPython Entrez failed: {e}")

    if not ncbi_success:
        # Fall back to wget
        logger.info("Falling back to NCBI E-utilities via wget")
        try:
            ncbi_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id={accession}&rettype=fasta&retmode=text"
            cmd = f"wget -q -O {output_file} \"{ncbi_url}\""
            subprocess.run(cmd, shell=True, check=True)
            with open(output_file) as f:
                first_line = f.readline().strip()
            if first_line.startswith(">"):
                ncbi_success = True
        except subprocess.SubprocessError:
            pass

    if not ncbi_success:
        # Fall back to SnpEff database sequences.fa (for custom genomes)
        logger.info(f"NCBI download failed for '{accession}', checking SnpEff database")
        snpeff_data_dirs = [
            os.environ.get("SNPEFF_DATA_CUSTOM", "/opt/vicast/snpeff_data_custom"),
            os.environ.get("SNPEFF_DATA", ""),
            os.environ.get("SNPEFF_DATA_BUILTIN", ""),
        ]
        for data_dir in snpeff_data_dirs:
            if not data_dir:
                continue
            seq_path = os.path.join(data_dir, accession, "sequences.fa")
            if os.path.exists(seq_path) and os.path.getsize(seq_path) > 0:
                logger.info(f"Found reference in SnpEff database: {seq_path}")
                shutil.copy2(seq_path, output_file)
                ncbi_success = True
                break

    if not ncbi_success:
        raise RuntimeError(
            f"Failed to obtain reference genome for '{accession}'. "
            f"Not a valid NCBI accession and not found in SnpEff database. "
            f"Use --reference to provide a local FASTA file."
        )

    # Verify the downloaded file
    if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
        raise RuntimeError(f"Reference genome file is empty or missing: {output_file}")

    logger.info(f"Reference genome ready: {output_file}")
    return output_file

def run_command(cmd: Union[str, List[str]], shell: bool = False, check: bool = True) -> subprocess.CompletedProcess:
    """
    Run a shell command and log the output.
    
    Args:
        cmd: Command to run (string or list of arguments)
        shell: Whether to run the command in a shell
        check: Whether to check the return code
        
    Returns:
        CompletedProcess instance with stdout and stderr
    """
    if isinstance(cmd, list):
        cmd_str = " ".join(cmd)
    else:
        cmd_str = cmd
        
    logger.debug(f"Running command: {cmd_str}")
    
    result = subprocess.run(
        cmd, 
        shell=shell, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        universal_newlines=True,
        check=False
    )
    
    if result.returncode != 0:
        logger.error(f"Command failed with exit code {result.returncode}")
        logger.error(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
        if check:
            raise subprocess.SubprocessError(f"Command failed: {cmd_str}")
    
    return result

def check_snpeff_database(accession: str, snpeff_jar: str, java_path: str = "java") -> bool:
    """
    Check if a genome is in the snpEff database.
    Checks custom data directories first (fast, works offline),
    then falls back to snpEff databases command for built-in genomes.

    Args:
        accession: Genome accession number
        snpeff_jar: Path to snpEff.jar
        java_path: Path to Java executable

    Returns:
        True if the genome is in the database, False otherwise
    """
    logger.info(f"Checking if {accession} is in snpEff database")

    # Check data directories directly (works for custom databases)
    snpeff_dir = os.path.dirname(snpeff_jar)
    for data_dir in [os.environ.get("SNPEFF_DATA_CUSTOM", ""),
                     os.environ.get("SNPEFF_DATA", ""),
                     os.path.join(snpeff_dir, "data")]:
        if data_dir and os.path.exists(os.path.join(data_dir, accession, "snpEffectPredictor.bin")):
            return True

    # Fallback: snpEff databases command (built-in genomes)
    cmd = f"{java_path} -jar {snpeff_jar} databases 2>/dev/null | grep {accession}"
    result = run_command(cmd, shell=True, check=False)
    return result.returncode == 0

def check_genome_complexity(accession: str, fasta_path: str) -> bool:
    """
    Check if a genome requires custom GFF3 processing due to complex annotation structure.
    
    Args:
        accession: Genome accession number
        fasta_path: Path to the genome FASTA file
        
    Returns:
        True if genome needs custom GFF3, False if standard GenBank processing is adequate
    """
    logger.info(f"Analyzing annotation complexity for genome: {accession}")
    
    # Download GenBank file to analyze annotation structure
    gb_path = fasta_path.replace(".fasta", ".gb")
    if not os.path.exists(gb_path):
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gbwithparts", retmode="text")
            with open(gb_path, 'w') as f:
                f.write(handle.read())
            handle.close()
        except:
            logger.warning(f"Could not download GenBank for complexity analysis. Using standard processing.")
            return False
    
    try:
        # Check for specific complex annotation patterns
        complex_patterns = [
            # Check for polyproteins
            "codon_start=1.*product=\".*polyprotein.*\"",
            # Check for mat_peptide features
            r"mat_peptide[ ]+[0-9]+\.\.[0-9]+",
            # Check for CDS with gene name in product
            r"CDS[ ]+[0-9]+\.\.[0-9]+.*product=\".*protein.*\"",
            # Check for multiple genes with same product
            "product=\"nonstructural polyprotein.*\".*product=\"nonstructural polyprotein",
            # Check for POLY gene
            "gene=\"POLY\"",
            # Check for Alphavirus-like patterns
            "product=\".*nsP[1-4].*\"",
            # Check for structural proteins common in viruses
            "product=\"capsid protein\"",
            "product=\"envelope glycoprotein\"",
            "product=\"membrane protein\"",
            "product=\"nucleocapsid protein\"",
            # Check for mature peptides by their common names
            "product=\"helicase\"",
            "product=\"protease\"",
            "product=\"RNA-dependent RNA polymerase\"",
            "product=\".*protease.*\"",
            "product=\".*replicase.*\"",
            # Check specifically for common virus families that need custom handling
            "Togaviridae|Picornaviridae|Flaviviridae|Coronaviridae|Hepeviridae|Astroviridae"
        ]
        
        # Check each pattern
        needs_custom_gff = False
        for pattern in complex_patterns:
            cmd = f"grep -E '{pattern}' {gb_path}"
            result = run_command(cmd, shell=True, check=False)
            if result.returncode == 0 and result.stdout.strip():
                logger.info(f"Complex annotation pattern detected: {pattern}")
                needs_custom_gff = True
                break
        
        if needs_custom_gff:
            # Get virus family for additional context
            cmd = f"grep -A2 'ORGANISM' {gb_path}"
            family_result = run_command(cmd, shell=True, check=False)
            if family_result.returncode == 0:
                logger.info(f"Taxonomy info: {family_result.stdout.strip()}")
                
            logger.info(f"Genome {accession} requires custom GFF3 processing due to complex annotation structure")
        else:
            logger.info(f"Genome {accession} has standard annotation structure - using regular GenBank processing")
            
        return needs_custom_gff
    
    except Exception as e:
        logger.warning(f"Error analyzing genome complexity: {str(e)}. Using standard processing.")
        return False
        
def generate_custom_gff3(accession: str, gb_path: str, output_dir: str) -> str:
    """
    Generate custom GFF3 for complex viral genomes.
    
    Args:
        accession: Genome accession number
        gb_path: Path to the GenBank file
        output_dir: Directory to save the GFF3 file
        
    Returns:
        Path to the generated GFF3 file
    """
    logger.info(f"Generating custom GFF3 for genome: {accession}")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Path to the GFF3 converter script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # First try the advanced converter (requires Biopython)
    converter_script = os.path.join(script_dir, "convert_gb_to_gff_advanced.py")
    
    # If the advanced script doesn't exist, use the simple version
    if not os.path.exists(converter_script):
        simple_converter_script = os.path.join(script_dir, "convert_gb_to_gff_simple.py")
        if os.path.exists(simple_converter_script):
            converter_script = simple_converter_script
            logger.info(f"Using simplified converter script: {converter_script}")
        else:
            # Look for any converter in common locations
            potential_paths = [
                "./convert_gb_to_gff_advanced.py",
                "./convert_gb_to_gff_simple.py",
                "../convert_gb_to_gff_advanced.py",
                "../convert_gb_to_gff_simple.py"
            ]
            
            for path in potential_paths:
                if os.path.exists(path):
                    converter_script = path
                    logger.info(f"Found converter script at: {converter_script}")
                    break
            
            # If still not found, create a simple version locally
            if not os.path.exists(converter_script):
                logger.warning("Creating simple GFF3 converter script...")
                try:
                    simple_script = """#!/usr/bin/env python3
import sys
import re
import os

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 convert_gb_to_gff_simple.py input.gb output.gff3")
        sys.exit(1)
    
    gb_file = sys.argv[1]
    gff_file = sys.argv[2]
    
    # Extract accession from GenBank file
    accession = os.path.basename(gb_file).split('.')[0]
    with open(gb_file, 'r') as f:
        for line in f:
            if line.startswith('ACCESSION'):
                accession = line.split()[1].strip()
                break
    
    # Create a simple GFF3 file
    with open(gff_file, 'w') as out:
        out.write('##gff-version 3\\n')
        out.write(f'##sequence-region {accession} 1 10000\\n')
        out.write(f'{accession}\\tGenBank\\tregion\\t1\\t10000\\t.\\t+\\t.\\tID=region_1\\n')
        
        # Add a generic CDS feature for snpEff
        out.write(f'{accession}\\tGenBank\\tCDS\\t1\\t10000\\t.\\t+\\t0\\tID=CDS_1;Name=viral_protein;product=viral protein\\n')

if __name__ == "__main__":
    main()
"""
                    simple_path = os.path.join(script_dir, "convert_gb_to_gff_simple.py")
                    with open(simple_path, 'w') as f:
                        f.write(simple_script)
                    os.chmod(simple_path, 0o755)  # Make executable
                    converter_script = simple_path
                    logger.info(f"Created simple converter script: {converter_script}")
                except:
                    logger.error("Could not create a simple converter script. Using standard processing.")
                    return None
    
    # Generate the GFF3 file
    gff3_path = os.path.join(output_dir, f"{accession}.gff3")
    
    try:
        cmd = f"python3 {converter_script} {gb_path} {gff3_path}"
        run_command(cmd, shell=True)
        
        # Verify the GFF3 file was created and is not empty
        if os.path.exists(gff3_path) and os.path.getsize(gff3_path) > 0:
            logger.info(f"Successfully generated custom GFF3 file: {gff3_path}")
            return gff3_path
        else:
            logger.error(f"GFF3 generation failed or produced empty file")
            return None
    except Exception as e:
        logger.error(f"Error generating custom GFF3: {str(e)}")
        return None

def add_genome_to_snpeff(accession: str, fasta_path: str, snpeff_jar: str, java_path: str = "java") -> bool:
    """
    Add a genome to the snpEff database.
    
    Args:
        accession: Genome accession number
        fasta_path: Path to the genome FASTA file
        snpeff_jar: Path to snpEff.jar
        java_path: Path to Java executable
        
    Returns:
        True if successful, False otherwise
    """
    logger.info(f"Adding genome {accession} to snpEff database")
    
    # Extract snpEff directory
    snpeff_dir = os.path.dirname(os.path.abspath(snpeff_jar))
    
    # Check for genome complexity
    needs_custom_gff = check_genome_complexity(accession, fasta_path)
    
    if needs_custom_gff:
        logger.info(f"Using custom GFF3 approach for genome {accession}")
        
        # Generate GenBank file path
        gb_path = fasta_path.replace(".fasta", ".gb")
        if not os.path.exists(gb_path):
            handle = Entrez.efetch(db="nucleotide", id=accession,
                                   rettype="gbwithparts", retmode="text")
            with open(gb_path, 'w') as f:
                f.write(handle.read())
            handle.close()
        
        # Generate custom GFF3
        temp_dir = os.path.dirname(fasta_path)
        gff3_path = generate_custom_gff3(accession, gb_path, temp_dir)
        
        if not gff3_path:
            logger.warning("Custom GFF3 generation failed, falling back to standard GenBank approach")
            cmd = f"{java_path} -jar {snpeff_jar} build -genbank -v -noCheckProtein {accession}"
            result = run_command(cmd, shell=True, check=False)
        else:
            # Set up snpEff with GFF3 and FASTA
            data_dir = os.path.join(snpeff_dir, "data", accession)
            os.makedirs(data_dir, exist_ok=True)
            
            # Copy FASTA file to snpEff data directory
            sequences_path = os.path.join(data_dir, "sequences.fa")
            shutil.copy2(fasta_path, sequences_path)
            
            # Copy GFF3 file to snpEff data directory
            genes_path = os.path.join(data_dir, "genes.gff")
            shutil.copy2(gff3_path, genes_path)
            
            # Modify snpEff.config
            config_path = os.path.join(snpeff_dir, "snpEff.config")
            
            # Ensure the config file has the genome entry
            # First check if accession already exists in config
            check_cmd = f"grep -E '^{accession}\\.genome:' {config_path}"
            check_result = run_command(check_cmd, shell=True, check=False)
            
            if check_result.returncode != 0:
                # Extract organism and strain from GenBank
                organism = accession
                try:
                    with open(gb_path, 'r') as f:
                        in_organism = False
                        for line in f:
                            if line.startswith('ORGANISM'):
                                in_organism = True
                                continue
                            elif in_organism and not line.startswith(' '):
                                in_organism = False
                                
                            if in_organism:
                                organism = line.strip()
                                break
                except:
                    logger.warning(f"Could not extract organism from GenBank file. Using accession as genome name.")
                
                # Try to extract strain
                strain = ""
                try:
                    strain_cmd = f"grep -E 'strain=\"[^\"]+\"' {gb_path} | head -1"
                    strain_result = run_command(strain_cmd, shell=True, check=False)
                    if strain_result.returncode == 0 and strain_result.stdout.strip():
                        strain_match = re.search(r'strain="([^"]+)"', strain_result.stdout)
                        if strain_match:
                            strain = strain_match.group(1)
                except:
                    logger.warning("Could not extract strain from GenBank file.")
                
                # Create abbreviated name
                abbr_name = accession
                if organism != accession:
                    # Extract first letters of each word in organism
                    try:
                        organism_abbr = ''.join([word[0].upper() for word in organism.split() if word and word[0].isalpha()])
                        if strain:
                            abbr_name = f"{organism_abbr}-{strain}"
                        else:
                            abbr_name = organism_abbr
                    except:
                        logger.warning("Error creating abbreviated name. Using accession.")
                
                # Add to the local config first (not the global snpEff config)
                local_config = os.path.join(os.path.dirname(fasta_path), "snpEff.config")
                with open(local_config, 'w') as f:
                    f.write(f"# {accession} genome configuration\n")
                    f.write(f"{accession}.genome: {abbr_name}\n")
                    f.write(f"{accession}.chromosomes: {accession}\n")
                    f.write(f"{accession}.codonTable: Standard\n")
                
                # Copy the global config to a backup in case we mess it up
                backup_config = f"{config_path}.bak"
                if not os.path.exists(backup_config):
                    shutil.copy2(config_path, backup_config)
                
                # Now add to the global config too
                try:
                    with open(config_path, 'a') as f:
                        f.write(f"\n# {accession}\n")
                        f.write(f"{accession}.genome: {abbr_name}\n")
                        f.write(f"{accession}.chromosomes: {accession}\n")
                        f.write(f"{accession}.codonTable: Standard\n")
                except Exception as e:
                    logger.error(f"Error updating snpEff config: {str(e)}")
                    logger.info(f"Using local config: {local_config}")
                    
                # For safety, also create a standalone config file for this run
                os.environ["SNPEFF_CONFIG"] = local_config
                logger.info(f"Set SNPEFF_CONFIG environment variable to: {local_config}")
            
            # Build the database with more forgiving parameters for viral genomes
            cmd = f"{java_path} -jar {snpeff_jar} build -gff3 -v -noCheckProtein -noCheckCds -noLog -treatAllAsProteinCoding {accession}"
            result = run_command(cmd, shell=True, check=False)
    else:
        # Use standard GenBank approach with forgiving parameters
        logger.info(f"Using standard GenBank approach for genome {accession}")
        cmd = f"{java_path} -jar {snpeff_jar} build -genbank -v -noCheckProtein -noCheckCds -noLog -treatAllAsProteinCoding {accession}"
        result = run_command(cmd, shell=True, check=False)
    
    # Check if build was successful
    success = result.returncode == 0
    if success:
        logger.info(f"Successfully added {accession} to snpEff database")
    else:
        logger.error(f"Failed to add {accession} to snpEff database")
        logger.error(f"STDOUT: {result.stdout}")
        logger.error(f"STDERR: {result.stderr}")
    
    return success

def calculate_read_stats(r1_pattern: str, output_file: str, threads: int = 1, r2_pattern: str = None) -> None:
    """
    Calculate read statistics using seqkit.
    
    Args:
        r1_pattern: Pattern to match R1 FASTQ files
        output_file: Output file for statistics
        threads: Number of threads to use
        r2_pattern: Optional pattern to match R2 FASTQ files
    """
    logger.info("Calculating read statistics")
    
    # Create patterns for both R1 and R2 files if needed
    patterns = [r1_pattern]
    if r2_pattern:
        patterns.append(r2_pattern)
    
    # Combine patterns for seqkit
    combined_pattern = " ".join(patterns)
    
    # Need to handle the case where the patterns might include shell special characters
    # If there are any potential glob patterns, use shell=True to let the shell expand them
    if any(char in combined_pattern for char in '*?[]{}'):
        cmd = f"seqkit stats {combined_pattern} -T -j {threads} > {output_file}"
        run_command(cmd, shell=True)
    else:
        # For non-glob patterns, we can use a more direct approach
        matching_files = []
        for pattern in patterns:
            if '*' in pattern:
                # Expand the glob pattern
                matching_files.extend(glob.glob(pattern))
            else:
                # Add the file directly if it exists
                if os.path.exists(pattern):
                    matching_files.append(pattern)
        
        # If we found files, run seqkit on them
        if matching_files:
            files_str = " ".join(matching_files)
            cmd = f"seqkit stats {files_str} -T -j {threads} > {output_file}"
            run_command(cmd, shell=True)
        else:
            logger.warning(f"No FASTQ files found matching patterns: {combined_pattern}")
            # Create an empty output file to prevent errors
            with open(output_file, 'w') as f:
                f.write("No matching FASTQ files found\n")
    
    logger.info(f"Read statistics saved to {output_file}")

def clean_reads(output_dir: str, r1_pattern: str, r2_pattern: str, threads: int = 1) -> Dict[str, Dict[str, str]]:
    """
    Clean reads using fastp.
    
    Args:
        output_dir: Output directory for cleaned reads
        r1_pattern: Pattern to match R1 FASTQ files
        r2_pattern: Pattern to match R2 FASTQ files
        threads: Number of threads to use
        
    Returns:
        Dictionary mapping sample names to input and output files
    """
    logger.info("Cleaning reads with fastp")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get R1 fastq files based on provided pattern
    # Handle both glob patterns and direct file paths
    if '*' in r1_pattern:
        # It's a glob pattern
        r1_files = glob.glob(r1_pattern)
    else:
        # It might be a specific file
        r1_files = [r1_pattern] if os.path.exists(r1_pattern) else []
    
    # Sort the files to ensure consistent processing
    r1_files.sort()
    
    if not r1_files:
        raise FileNotFoundError(f"No R1 FASTQ files found matching pattern: {r1_pattern}")
    
    logger.info(f"Found {len(r1_files)} R1 FASTQ files")
    for i, f in enumerate(r1_files):
        logger.info(f"  {i+1}. {f}")
    
    cleaned_files = {}
    
    for r1_file in r1_files:
        # Extract sample name from the filename
        r1_basename = os.path.basename(r1_file)

        # Match the sample name part before _R1/_R1_001 (Illumina) or _1 (SRA)
        # Try Illumina pattern first (_R1 or _R1_001)
        match = re.search(r'(.+?)_R1(?:_001)?\.fastq\.gz$', r1_basename)
        r1_suffix = '_R1'
        r2_suffix = '_R2'

        # If no match, try SRA pattern (_1 or _2)
        if not match:
            match = re.search(r'(.+?)_1\.fastq\.gz$', r1_basename)
            if match:
                r1_suffix = '_1'
                r2_suffix = '_2'

        if not match:
            logger.warning(f"Skipping file with unusual naming pattern: {r1_file}")
            logger.warning(f"  Expected patterns: *_R1.fastq.gz, *_R1_001.fastq.gz (Illumina), or *_1.fastq.gz (SRA)")
            continue

        sample_name = match.group(1)

        # Find matching R2 file - handle both patterns and paths
        r1_dir = os.path.dirname(r1_file)
        r2_basename = r1_basename.replace(r1_suffix, r2_suffix)
        r2_file = os.path.join(r1_dir, r2_basename)

        # If exact file doesn't exist and r2_pattern contains wildcards
        if not os.path.exists(r2_file) and '*' in r2_pattern:
            # Try to find a matching R2 file using the glob pattern
            r2_candidates = glob.glob(r2_pattern)
            for candidate in r2_candidates:
                candidate_basename = os.path.basename(candidate)
                # Check for matching sample name and R2 suffix
                if candidate_basename.startswith(sample_name) and r2_suffix in candidate_basename:
                    r2_file = candidate
                    break

        if not os.path.exists(r2_file):
            logger.warning(f"No matching R2 file found for {r1_file}")
            logger.warning(f"  Expected R2 file: {r2_file}")
            continue
        
        logger.info(f"Processing sample: {sample_name}")
        logger.info(f"  R1: {r1_file}")
        logger.info(f"  R2: {r2_file}")
        
        # Output file paths - standardized naming convention
        r1_out = os.path.join(output_dir, f"{sample_name}_R1.qc.fastq.gz")
        r2_out = os.path.join(output_dir, f"{sample_name}_R2.qc.fastq.gz")
        html_report = os.path.join(output_dir, f"{sample_name}_fastp_report.html")
        json_report = os.path.join(output_dir, f"{sample_name}_fastp_report.json")
        qc_report = os.path.join(output_dir, f"{sample_name}_fastp_QC.txt")
        
        # Run fastp
        cmd = (
            f"fastp -w {threads} -q 30 "
            f"-i {r1_file} -I {r2_file} "
            f"-o {r1_out} -O {r2_out} "
            f"-h {html_report} -j {json_report}"
        )
        
        # Capture fastp output for QC report
        result = run_command(cmd, shell=True)
        
        # Extract and save QC information
        qc_info = []
        qc_info.append(f"{sample_name}")
        
        output_lines = result.stderr.split('\n')
        capture = False
        for line in output_lines:
            if "Filtering result:" in line:
                capture = True
                qc_info.append(line)
            elif capture and line.strip():
                qc_info.append(line)
            elif capture and "Insert size peak" in line:
                qc_info.append(line)
                break
        
        with open(qc_report, 'w') as f:
            f.write('\n'.join(qc_info))
        
        # Add files to result dictionary
        cleaned_files[sample_name] = {
            'r1_in': r1_file,
            'r2_in': r2_file,
            'r1_out': r1_out,
            'r2_out': r2_out,
            'qc_report': qc_report
        }
    
    logger.info(f"Cleaned {len(cleaned_files)} samples")
    return cleaned_files

def map_and_call_variants(
    reference: str,
    output_dir: str,
    threads: int = 1,
    cleaned_files: Dict[str, Tuple[str, str]] = None,
    large_files: bool = False,
    extremely_large_files: bool = False,
    skip_vector_filter: bool = False
) -> Dict[str, Dict[str, str]]:
    """
    Map reads to reference and call variants.

    Args:
        reference: Path to reference genome
        output_dir: Base output directory
        threads: Number of threads to use
        cleaned_files: Dictionary mapping sample names to (R1, R2) cleaned file paths
        skip_vector_filter: If True, skip cloning vector read filtering

    Returns:
        Dictionary mapping sample names to output files
    """
    logger.info(f"Mapping reads to reference and calling variants: {reference}")
    
    # Create output directories
    mapping_dir = os.path.join(output_dir, "mapping")
    variants_dir = os.path.join(output_dir, "variants")
    os.makedirs(mapping_dir, exist_ok=True)
    os.makedirs(variants_dir, exist_ok=True)
    
    # Index reference genome
    logger.info(f"Indexing reference genome: {reference}")
    run_command(f"bwa index {reference}", shell=True)
    
    # If cleaned_files provided, use those; otherwise find all in directory
    if cleaned_files:
        # Use the specific cleaned files provided
        files_to_process = cleaned_files
    else:
        # Legacy behavior: find all cleaned files in directory
        r1_files = [f for f in os.listdir(output_dir) if f.endswith('_R1.qc.fastq.gz')]
        
        if not r1_files:
            raise FileNotFoundError(f"No cleaned R1 files found in {output_dir}")
        
        files_to_process = {}
        for r1_file in r1_files:
            sample_name = r1_file.replace('_R1.qc.fastq.gz', '')
            r2_file = r1_file.replace('_R1', '_R2')
            r1_path = os.path.join(output_dir, r1_file)
            r2_path = os.path.join(output_dir, r2_file)
            files_to_process[sample_name] = (r1_path, r2_path)
    
    result_files = {}
    
    # Process each sample
    for sample_name, file_info in files_to_process.items():
        # Handle both tuple and dict formats
        if isinstance(file_info, tuple):
            r1_path, r2_path = file_info
        else:
            r1_path = file_info['r1_out']
            r2_path = file_info['r2_out']
        
        if not os.path.exists(r2_path):
            logger.warning(f"No matching R2 file found for {r1_path}")
            continue
        
        logger.info(f"Processing sample: {sample_name}")

        # ── Vector filtering: remove reads matching cloning vectors ──
        novector_r1 = None
        novector_r2 = None
        if not skip_vector_filter:
            vector_db = os.path.join(
                os.environ.get('BLAST_DB_DIR', '/opt/vicast/blast_db'),
                'cloning_vectors.fasta')
            if os.path.exists(vector_db):
                logger.info(f"Filtering cloning vector reads for {sample_name}")
                novector_r1 = os.path.join(mapping_dir, f"{sample_name}_novector_R1.fastq.gz")
                novector_r2 = os.path.join(mapping_dir, f"{sample_name}_novector_R2.fastq.gz")

                # Index vector DB if not already indexed
                if not os.path.exists(vector_db + '.bwt'):
                    run_command(f"bwa index {vector_db}", shell=True)

                # Map to vectors, keep only unmapped read pairs (-f 12)
                run_command(
                    f"bwa mem -t {threads} {vector_db} {r1_path} {r2_path} | "
                    f"samtools fastq -f 12 -F 256 "
                    f"-1 {novector_r1} -2 {novector_r2} -s /dev/null -",
                    shell=True
                )

                # Log filtering stats
                try:
                    import subprocess as _sp
                    orig_count = _sp.run(
                        f"zcat {r1_path} | wc -l",
                        shell=True, capture_output=True, text=True
                    ).stdout.strip()
                    filt_count = _sp.run(
                        f"zcat {novector_r1} | wc -l",
                        shell=True, capture_output=True, text=True
                    ).stdout.strip()
                    orig_reads = int(orig_count) // 4
                    filt_reads = int(filt_count) // 4
                    removed = orig_reads - filt_reads
                    logger.info(f"Vector filtering: {orig_reads} → {filt_reads} read pairs "
                                f"({removed} removed, {removed/max(orig_reads,1)*100:.2f}%)")
                except Exception:
                    pass

                # Use filtered reads for all downstream steps
                r1_path = novector_r1
                r2_path = novector_r2
            else:
                logger.debug("No vector database found — skipping vector filtering")

        # Prepare output file paths
        sam_file = os.path.join(mapping_dir, f"{sample_name}.sam")
        fixmate_file = os.path.join(mapping_dir, f"{sample_name}.fixmate.bam")
        bam_file = os.path.join(mapping_dir, f"{sample_name}.bam")
        dedupe_file = os.path.join(mapping_dir, f"{sample_name}.dedupe.bam")
        realign_file = os.path.join(mapping_dir, f"{sample_name}.lofreq.realign.bam")
        indel_file = os.path.join(mapping_dir, f"{sample_name}.lofreq.indel.bam")
        final_bam = os.path.join(mapping_dir, f"{sample_name}.lofreq.final.bam")
        final_bai = os.path.join(mapping_dir, f"{sample_name}.lofreq.final.bam.bai")
        vars_file = os.path.join(variants_dir, f"{sample_name}_vars.vcf")
        
        # Alignment steps
        
        # Check if we should use streaming for extremely large files
        if extremely_large_files:
            # Stream everything to avoid intermediate files for extremely large datasets
            logger.info(f"Using streaming pipeline for extremely large files - {sample_name}")
            sort_mem = "-m 32G"  # 32GB per thread for extremely large files
            
            # Create a single piped command: align -> fixmate -> sort -> markdup
            logger.info(f"Running streaming alignment pipeline for {sample_name}")
            streaming_cmd = (
                f"bwa mem -t {threads} {reference} {r1_path} {r2_path} | "
                f"samtools fixmate -m -O bam --threads {threads} - - | "
                f"samtools sort {sort_mem} --threads {threads} -O bam - | "
                f"samtools markdup --threads {threads} -S - {dedupe_file}"
            )
            run_command(streaming_cmd, shell=True)
            
            # For compatibility, create sorted bam reference
            bam_file = dedupe_file
            
        else:
            # Original non-streaming approach for normal files
            # 1. BWA MEM alignment
            logger.info(f"Aligning reads for {sample_name}")
            run_command(
                f"bwa mem -t {threads} {reference} {r1_path} {r2_path} > {sam_file}",
                shell=True
            )
            
            # 2. Fix mate information
            logger.info(f"Fixing mate information for {sample_name}")
            run_command(
                f"samtools fixmate -O bam -m --threads {threads} {sam_file} {fixmate_file}",
                shell=True
            )
            
            # 3. Sort BAM
            logger.info(f"Sorting BAM file for {sample_name}")
            # Use more memory for large files
            if large_files:
                sort_mem = "-m 8G"   # 8GB per thread for large files
            else:
                sort_mem = ""        # Default for normal files
            run_command(
                f"samtools sort {sort_mem} --threads {threads} -O bam {fixmate_file} > {bam_file}",
                shell=True
            )
            
            # 4. Mark duplicates
            logger.info(f"Marking duplicates for {sample_name}")
            run_command(
                f"samtools markdup --threads {threads} -S {bam_file} {dedupe_file}",
                shell=True
            )
        
        # 5. LoFreq Viterbi realignment
        logger.info(f"LoFreq Viterbi realignment for {sample_name}")
        if extremely_large_files:
            sort_mem = "-m 32G"  # Already set above
        elif large_files:
            sort_mem = "-m 8G"
        else:
            sort_mem = ""
        run_command(
            f"lofreq viterbi -f {reference} {dedupe_file} | samtools sort - {sort_mem} --threads {threads} > {realign_file}",
            shell=True
        )
        
        # 6. LoFreq indel quality calibration
        logger.info(f"LoFreq indel quality calibration for {sample_name}")
        run_command(
            f"lofreq indelqual --dindel -f {reference} {realign_file} | samtools sort - {sort_mem} --threads {threads} > {indel_file}",
            shell=True
        )
        
        # 7. LoFreq alignment quality calibration
        logger.info(f"LoFreq alignment quality calibration for {sample_name}")
        run_command(
            f"lofreq alnqual -b {indel_file} {reference} > {final_bam}",
            shell=True
        )
        
        # 8. Index final BAM
        logger.info(f"Indexing final BAM for {sample_name}")
        run_command(f"samtools index {final_bam}", shell=True)
        
        # 9. Call variants with LoFreq
        logger.info(f"Calling variants for {sample_name}")
        # lofreq filter (invoked internally by call-parallel) has no
        # --force-overwrite flag and will abort if the output file exists.
        # Unconditionally remove any stale VCF from a previous run.
        for stale in [vars_file, vars_file + ".gz"]:
            try:
                os.remove(stale)
                logger.info(f"Removed pre-existing variant file: {stale}")
            except FileNotFoundError:
                pass

        # Call variants with LoFreq
        run_command(
            f"lofreq call-parallel --pp-threads {threads} --force-overwrite --no-default-filter --call-indels -f {reference} -o {vars_file} {final_bam}",
            shell=True
        )
        
        # Clean up vector-filtered intermediate files
        for tmp_file in [novector_r1, novector_r2]:
            if tmp_file and os.path.exists(tmp_file):
                os.remove(tmp_file)
                logger.debug(f"Cleaned up: {tmp_file}")

        # Store output files
        result_files[sample_name] = {
            'sam': sam_file,
            'fixmate': fixmate_file,
            'bam': bam_file,
            'dedupe': dedupe_file,
            'realign': realign_file,
            'indel': indel_file,
            'final_bam': final_bam,
            'final_bai': final_bai,
            'vars': vars_file
        }
    
    logger.info(f"Mapping and variant calling completed for {len(result_files)} samples")
    return result_files

def filter_variants(variants_dir: str, specific_files: Dict[str, Dict[str, str]], min_depth: int = 200, min_qual: int = 1000) -> Dict[str, str]:
    """
    Filter variant calls from LoFreq.

    Args:
        variants_dir: Directory containing variant files
        specific_files: Dictionary of specific files to filter (from variant calling step) - REQUIRED
        min_depth: Minimum read depth for variant filter (lofreq -v)
        min_qual: Minimum variant quality in phred (lofreq -Q for SNVs, -K for indels)

    Returns:
        Dictionary mapping sample names to filtered variant files
    """
    logger.info("Filtering variant calls")
    
    if not specific_files:
        raise ValueError("specific_files parameter is required - cannot process all files in directory")
    
    # Use only the specific files from the current run
    vcf_files = []
    for sample_name, file_info in specific_files.items():
        if 'vars' in file_info:
            vcf_filename = os.path.basename(file_info['vars'])
            vcf_files.append(vcf_filename)
    
    if not vcf_files:
        raise FileNotFoundError(f"No variant VCF files found in specific_files for {variants_dir}")
    
    filtered_files = {}
    
    for vcf_file in vcf_files:
        # Extract sample name
        sample_name = vcf_file.replace('_vars.vcf', '')
        vcf_path = os.path.join(variants_dir, vcf_file)
        
        filtered_path = os.path.join(variants_dir, f"{sample_name}_vars.filt.vcf")
        
        logger.info(f"Filtering variants for {sample_name}")

        # Check if output file exists and remove it
        if os.path.exists(filtered_path):
            logger.info(f"Removing existing filtered file: {filtered_path}")
            os.remove(filtered_path)
      
        # Run LoFreq filter with depth and quality thresholds
        logger.info(f"  Filters: min_depth={min_depth}, min_qual(SNV)={min_qual}, min_qual(indel)={min_qual}")
        run_command(
            f"lofreq filter -i {vcf_path} -o {filtered_path} -v {min_depth} -Q {min_qual} -K {min_qual}",
            shell=True
        )
        
        filtered_files[sample_name] = filtered_path
    
    logger.info(f"Variant filtering completed for {len(filtered_files)} samples")
    return filtered_files

def annotate_variants(variants_dir: str, accession: str, snpeff_jar: str, java_path: str = "java", specific_files: Dict[str, str] = None, large_files: bool = False, extremely_large_files: bool = False) -> Dict[str, Dict[str, str]]:
    """
    Annotate filtered variants using snpEff.
    
    Args:
        variants_dir: Directory containing filtered variant files
        accession: Reference genome accession
        snpeff_jar: Path to snpEff.jar
        java_path: Path to Java executable
        specific_files: Dictionary of specific filtered files to annotate (from filter step) - REQUIRED
        
    Returns:
        Dictionary mapping sample names to annotation files
    """
    logger.info(f"Annotating variants with snpEff using database: {accession}")
    
    if not specific_files:
        raise ValueError("specific_files parameter is required - cannot process all files in directory")
    
    # Use only the specific files from the current run
    filt_files = []
    for sample_name, filtered_path in specific_files.items():
        filt_filename = os.path.basename(filtered_path)
        filt_files.append(filt_filename)
    
    if not filt_files:
        raise FileNotFoundError(f"No filtered variant VCF files found in specific_files for {variants_dir}")
    
    # Verify snpEff database has the genome by checking for database files directly
    # (snpEff databases command only lists built-in DBs, not custom ones in data/)
    snpeff_dir = os.path.dirname(snpeff_jar)
    snpeff_data_dirs = [
        os.environ.get("SNPEFF_DATA_CUSTOM", ""),
        os.environ.get("SNPEFF_DATA", ""),
        os.path.join(snpeff_dir, "data"),
    ]
    db_found = False
    for data_dir in snpeff_data_dirs:
        if not data_dir:
            continue
        predictor = os.path.join(data_dir, accession, "snpEffectPredictor.bin")
        if os.path.exists(predictor):
            logger.info(f"Verified genome {accession} exists in snpEff database: {data_dir}/{accession}/")
            db_found = True
            break

    if not db_found:
        # Fallback: try snpEff databases command (catches built-in genomes)
        verify_cmd = f"{java_path} -jar {snpeff_jar} databases 2>/dev/null | grep -i {accession}"
        verify_result = run_command(verify_cmd, shell=True, check=False)
        if verify_result.returncode == 0:
            logger.info(f"Verified genome {accession} exists in snpEff built-in database.")
            db_found = True

    if not db_found:
        checked = [d for d in snpeff_data_dirs if d]
        logger.error(f"CRITICAL ERROR: Genome {accession} NOT FOUND in snpEff database!")
        logger.error(f"Checked directories: {checked}")
        raise RuntimeError(f"Genome {accession} not found in snpEff database. Use --add-to-snpeff to add it.")
    
    annotation_files = {}
    
    for filt_file in filt_files:
        # Extract sample name
        sample_name = filt_file.replace('_vars.filt.vcf', '')
        filt_path = os.path.join(variants_dir, filt_file)
        
        # Output file paths
        ann_vcf = os.path.join(variants_dir, f"{sample_name}.snpEFF.ann.vcf")
        ann_tsv = os.path.join(variants_dir, f"{sample_name}.snpEFF.ann.tsv")
        summary_html = os.path.join(variants_dir, f"{sample_name}_summary.html")
        summary_genes = os.path.join(variants_dir, f"{sample_name}_summary.genes.txt")
        
        logger.info(f"Annotating variants for {sample_name}")

        # Check if output files exist and remove them
        if os.path.exists(ann_vcf):
            logger.info(f"Removing existing annotation file: {ann_vcf}")
            os.remove(ann_vcf)
      
        # Fix VCF file to remove problematic IUB ambiguity codes that can cause snpEff to fail
        # IUB codes like R, Y, S, W, K, M etc. in ALT field cause snpEff errors
        logger.info(f"Fixing VCF file to remove problematic IUB ambiguity codes: {filt_path}")
        try:
            valid_bases = set("ACGTNacgtn.,*")
            fixed_lines = []
            removed_count = 0
            with open(filt_path, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        fixed_lines.append(line)
                    else:
                        fields = line.rstrip('\n').split('\t')
                        if len(fields) >= 5:
                            alt = fields[4]
                            # Check if ALT contains only valid bases (including multi-allelic commas)
                            if all(c in valid_bases or c == ',' for c in alt):
                                fixed_lines.append(line)
                            else:
                                removed_count += 1
                        else:
                            fixed_lines.append(line)
            if removed_count > 0:
                logger.info(f"Removed {removed_count} variants with IUB ambiguity codes")
                with open(filt_path, 'w') as f:
                    f.writelines(fixed_lines)
            logger.info(f"VCF file fixed successfully: {filt_path}")
        except Exception as e:
            logger.warning(f"Could not fix VCF file due to error: {str(e)}")
            logger.warning("Proceeding with original VCF file - snpEff may fail if problematic variants are present")
        
        # Validate the filtered VCF file
        variant_count_in_input = 0
        with open(filt_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    variant_count_in_input += 1
        
        logger.info(f"Input VCF contains {variant_count_in_input} variants")
        
        if variant_count_in_input == 0:
            logger.warning(f"Input VCF file {filt_path} has no variants! Creating empty output files.")
            # Create empty files
            with open(ann_vcf, 'w') as f:
                f.write("")
            header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tEFFECT\tPUTATIVE_IMPACT\tGENE_NAME\tGENE_ID\tFEATURE_TYPE\tFEATURE_ID\tTRANSCRIPT_TYPE\tEXON_INTRON_RANK\tHGVSc\tHGVSp\tcDNA_POSITION_AND_LENGTH\tCDS_POSITION_AND_LENGTH\tPROTEIN_POSITION_AND_LENGTH\tDISTANCE_TO_FEATURE\tERROR"
            with open(ann_tsv, 'w') as f:
                f.write(header + "\n")
            # Add to annotation_files even for empty results
            annotation_files[sample_name] = {
                'ann_vcf': ann_vcf,
                'ann_tsv': ann_tsv,
                'summary_html': summary_html,
                'summary_genes': summary_genes
            }
            continue

        # Make a copy of the VCF file for safety
        safe_filt_path = f"{filt_path}.safe"
        shutil.copy2(filt_path, safe_filt_path)

        # Run snpEff with more robust error handling and debug output
        # Use higher memory for large files
        if extremely_large_files:
            java_mem = "-Xmx64g"   # 64GB for extremely large files
        elif large_files:
            java_mem = "-Xmx32g"   # 32GB for large files
        else:
            java_mem = "-Xmx4g"    # 4GB default
        # Build snpEff command - use custom config if available (has genome entries for custom DBs)
        snpeff_config_flag = ""
        for candidate_config in [
            os.environ.get("SNPEFF_CONFIG_CUSTOM", ""),
            os.environ.get("SNPEFF_CONFIG", ""),
            os.path.join(os.environ.get("SNPEFF_DATA_CUSTOM", ""), "snpEff.config"),
        ]:
            if candidate_config and os.path.exists(candidate_config):
                snpeff_config_flag = f"-c {candidate_config}"
                logger.info(f"Using snpEff config: {candidate_config}")
                break

        cmd = f"{java_path} -jar {java_mem} {snpeff_jar} {snpeff_config_flag} -v {accession} {safe_filt_path} -s {summary_html} > {ann_vcf}"
        logger.info(f"Running command: {cmd}")
        
        try:
            result = run_command(cmd, shell=True)
            logger.info(f"SnpEff command completed successfully")
            logger.debug(f"STDOUT: {result.stdout}")
            logger.debug(f"STDERR: {result.stderr}")
        except subprocess.SubprocessError as e:
            logger.error(f"SnpEff command failed: {str(e)}")
            
            # Check if this is the IUB code error
            if "Unkown IUB code for SNP" in str(e) or "invalid start byte" in str(e):
                logger.warning("SnpEff failed due to IUB codes - attempting additional filtering")
                
                # Try more aggressive filtering
                try:
                    filtered_tmp = f"{safe_filt_path}.filtered"
                    
                    # Extract header
                    run_command(f"grep '^#' {safe_filt_path} > {filtered_tmp}", shell=True)
                    
                    # Extract and filter content lines more aggressively
                    aggressive_filter = (
                        f"grep -v '^#' {safe_filt_path} | "
                        f"grep -E '^[^[:space:]]+[[:space:]]+[0-9]+[[:space:]]+\\.[[:space:]]+[ACGT][[:space:]]+[ACGT][[:space:]]' >> {filtered_tmp}"
                    )
                    run_command(aggressive_filter, shell=True, check=False)
                    
                    # Check if we still have variants
                    count_cmd = f"grep -v '^#' {filtered_tmp} | wc -l"
                    count_result = run_command(count_cmd, shell=True, check=False)
                    aggressive_count = int(count_result.stdout.strip())
                    
                    if aggressive_count > 0:
                        logger.info(f"Filtered VCF now has {aggressive_count} variants (originally {variant_count_in_input})")
                        
                        # Try snpEff with the more aggressively filtered file
                        retry_cmd = f"{java_path} -jar {java_mem} {snpeff_jar} {snpeff_config_flag} -v {accession} {filtered_tmp} -s {summary_html} > {ann_vcf}"
                        try:
                            logger.info(f"Retrying with aggressively filtered VCF: {retry_cmd}")
                            run_command(retry_cmd, shell=True)
                        except subprocess.SubprocessError as retry_error:
                            logger.error(f"Retry also failed: {str(retry_error)}")
                            # Create empty output file
                            with open(ann_vcf, 'w') as f:
                                f.write("")
                    else:
                        logger.warning("No variants left after aggressive filtering!")
                        # Create empty output file
                        with open(ann_vcf, 'w') as f:
                            f.write("")
                except Exception as filter_error:
                    logger.error(f"Additional filtering failed: {str(filter_error)}")
                    # Create empty output file
                    with open(ann_vcf, 'w') as f:
                        f.write("")
            else:
                # Create empty output file for this case too
                with open(ann_vcf, 'w') as f:
                    f.write("")
                logger.error(f"Error not related to IUB codes: {str(e)}")
        
        # Clean up temporary files
        if os.path.exists(safe_filt_path):
            os.remove(safe_filt_path)
        
        # Process the annotated VCF file into TSV format (pure Python - no shell dependency)
        logger.info(f"Converting VCF to TSV for {sample_name}")

        ann_info_header = ["INFO", "EFFECT", "PUTATIVE_IMPACT", "GENE_NAME", "GENE_ID",
                           "FEATURE_TYPE", "FEATURE_ID", "TRANSCRIPT_TYPE", "EXON_INTRON_RANK",
                           "HGVSc", "HGVSp", "cDNA_POSITION_AND_LENGTH", "CDS_POSITION_AND_LENGTH",
                           "PROTEIN_POSITION_AND_LENGTH", "DISTANCE_TO_FEATURE", "ERROR"]
        base_header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER"]
        full_header = base_header + ann_info_header

        try:
            data_lines = []
            with open(ann_vcf, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.rstrip('\n').split('\t')
                    if len(fields) < 8:
                        continue
                    base_cols = fields[:7]
                    info_field = fields[7]
                    # Split INFO field on '|' to extract snpEff annotation columns
                    info_parts = info_field.split('|')
                    # Take first 16 fields (INFO + 15 annotation columns)
                    if len(info_parts) >= 16:
                        ann_cols = info_parts[:16]
                    else:
                        ann_cols = info_parts + [''] * (16 - len(info_parts))
                    data_lines.append(base_cols + ann_cols)

            logger.info(f"Found {len(data_lines)} annotated variants in output VCF")

            with open(ann_tsv, 'w') as f:
                # Prefix header with # so downstream Perl parser skips it
                f.write('#' + '\t'.join(full_header) + '\n')
                for row in data_lines:
                    f.write('\t'.join(row) + '\n')
        except Exception as e:
            logger.error(f"Error converting VCF to TSV: {str(e)}")
            logger.warning(f"Creating basic TSV file as fallback")
            with open(ann_tsv, 'w') as f:
                f.write('#' + '\t'.join(full_header) + '\n')
        
        # Check for ERROR_CHROMOSOME_NOT_FOUND
        if os.path.exists(ann_tsv) and os.path.getsize(ann_tsv) > 0:
            check_cmd = f"grep -q 'ERROR_CHROMOSOME_NOT_FOUND' {ann_tsv}"
            result = run_command(check_cmd, shell=True, check=False)
            
            if result.returncode == 0:
                logger.error(f"ERROR: Genome {accession} is NOT properly formatted for use with snpEff annotation.")
                logger.error("This often happens when the chromosome name in the VCF doesn't match what's in the snpEff database.")
                
                # Look at chromosome names in the VCF
                check_chrom_cmd = f"grep -v '^#' {filt_path} | head -1 | cut -f1"
                chrom_result = run_command(check_chrom_cmd, shell=True, check=False)
                if chrom_result.returncode == 0 and chrom_result.stdout.strip():
                    logger.error(f"Chromosome name in VCF: {chrom_result.stdout.strip()}")
                    logger.error(f"Make sure this matches what's in the snpEff database for {accession}.")
        
        annotation_files[sample_name] = {
            'ann_vcf': ann_vcf,
            'ann_tsv': ann_tsv,
            'summary_html': summary_html,
            'summary_genes': summary_genes
        }
    
    logger.info(f"Variant annotation completed for {len(annotation_files)} samples")
    return annotation_files

def parse_annotations(variants_dir: str, min_depth: int, specific_files: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    """
    Parse annotated variants using the Perl script.
    
    Args:
        variants_dir: Directory containing annotation files
        min_depth: Minimum read depth for reporting
        specific_files: Dictionary of specific annotation files to parse (from annotation step) - REQUIRED
        
    Returns:
        Dictionary mapping sample names to parsed annotation files
    """
    logger.info(f"Parsing annotated variants with minimum depth {min_depth}")
    
    if not specific_files:
        raise ValueError("specific_files parameter is required - cannot process all files in directory")
    
    # Use only the specific files from the current run
    ann_files = []
    for sample_name, file_info in specific_files.items():
        if 'ann_tsv' in file_info:
            ann_filename = os.path.basename(file_info['ann_tsv'])
            ann_files.append(ann_filename)
    
    if not ann_files:
        raise FileNotFoundError(f"No annotation TSV files found in specific_files for {variants_dir}")
    
    parsed_files = {}
    
    # Find the Perl script by looking in multiple locations
    # Get the directory where this Python script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Potential locations for the Perl script
    perl_script_locations = [
        # Current directory
        "./parse_snpEff_annotated_vcf_for_collaborators.pl",
        # Same directory as this Python script
        os.path.join(script_dir, "parse_snpEff_annotated_vcf_for_collaborators.pl"),
        # Parent directory of this Python script
        os.path.join(os.path.dirname(script_dir), "parse_snpEff_annotated_vcf_for_collaborators.pl"),
    ]
    
    perl_script = None
    for location in perl_script_locations:
        if os.path.exists(location):
            perl_script = location
            logger.info(f"Found Perl script at: {perl_script}")
            break
    
    # Check if Perl script was found
    if not perl_script:
        raise FileNotFoundError(
            "Perl script 'parse_snpEff_annotated_vcf_for_collaborators.pl' not found. "
            "Please make sure it's in the same directory as the viral_pipeline.py script "
            "or specify its path manually."
        )
    
    for ann_file in ann_files:
        # Extract sample name
        sample_name = ann_file.replace('.snpEFF.ann.tsv', '')
        ann_path = os.path.join(variants_dir, ann_file)
        
        # Parse with Perl script
        logger.info(f"Parsing annotations for {sample_name} with minimum depth {min_depth}")
        
        run_command(
            f"perl {perl_script} {ann_path} {min_depth}",
            shell=True
        )
        
        # The output file should be named by the script
        output_file = f"{ann_path.rstrip('.tsv')}_{min_depth}.tsv"
        parsed_files[sample_name] = output_file
    
    logger.info(f"Annotation parsing completed for {len(parsed_files)} samples")
    return parsed_files

def discover_existing_vcf_files(output_dir: str, r1_filename: str) -> Dict[str, Dict[str, str]]:
    """
    Discover existing VCF files for a specific sample from a previous QC run.

    This function finds VCF files for the sample specified by the R1 filename,
    matching the behavior of the normal workflow which processes only the
    specified sample.

    Args:
        output_dir: Base output directory (default is current directory)
        r1_filename: R1 FASTQ filename to determine which sample to process

    Returns:
        Dictionary mapping sample name to file paths, format:
        {
            'sample_name': {
                'vars': '/path/to/sample_vars.vcf',
                'final_bam': '/path/to/sample.lofreq.final.bam',
                'final_bai': '/path/to/sample.lofreq.final.bam.bai'
            }
        }

    Raises:
        FileNotFoundError: If the specified sample's VCF file is not found
    """
    logger.info("Discovering existing VCF files from previous QC run...")

    # Extract sample name from R1 filename (same logic as main pipeline)
    # Support both Illumina (_R1) and SRA (_1) naming conventions
    r1_base = os.path.basename(r1_filename)
    sample_name = re.sub(r'(_R?[12])?(_001)?\.fastq\.gz$', '', r1_base)

    logger.info(f"Looking for VCF files for sample: {sample_name}")

    # Define expected directory structure
    cleaned_dir = os.path.join(output_dir, "cleaned_seqs")
    variants_dir = os.path.join(cleaned_dir, "variants")
    mapping_dir = os.path.join(cleaned_dir, "mapping")

    # Check if variants directory exists
    if not os.path.exists(variants_dir):
        raise FileNotFoundError(
            f"Variants directory not found: {variants_dir}\n"
            "Please run the QC workflow first (Steps 1-6) using:\n"
            f"  ./run_vicast_analyze_qc_only.sh {r1_filename} <R2> <accession>"
        )

    # Look for VCF file for this specific sample
    vcf_file = os.path.join(variants_dir, f"{sample_name}_vars.vcf")

    if not os.path.exists(vcf_file):
        raise FileNotFoundError(
            f"VCF file not found for sample '{sample_name}': {vcf_file}\n"
            "Please run the QC workflow first (Steps 1-6) using:\n"
            f"  ./run_vicast_analyze_qc_only.sh {r1_filename} <R2> <accession>"
        )

    logger.info(f"Found VCF file for sample: {sample_name}")

    # Look for corresponding BAM files in mapping directory
    final_bam = os.path.join(mapping_dir, f"{sample_name}.lofreq.final.bam")
    final_bai = os.path.join(mapping_dir, f"{sample_name}.lofreq.final.bam.bai")

    # Check if BAM files exist (needed for depth calculation)
    if not os.path.exists(final_bam):
        logger.warning(
            f"BAM file not found for sample {sample_name}: {final_bam}\n"
            "Depth file generation may fail for this sample."
        )

    if not os.path.exists(final_bai):
        logger.warning(
            f"BAM index not found for sample {sample_name}: {final_bai}\n"
            f"You may need to regenerate it with: samtools index {final_bam}"
        )

    # Build result structure (matching map_and_call_variants output)
    result_files = {
        sample_name: {
            'vars': vcf_file,
            'final_bam': final_bam if os.path.exists(final_bam) else None,
            'final_bai': final_bai if os.path.exists(final_bai) else None
        }
    }

    logger.info(f"  → Discovered sample: {sample_name}")
    logger.info(f"      VCF: {vcf_file}")
    if os.path.exists(final_bam):
        logger.info(f"      BAM: {final_bam}")

    logger.info(f"Successfully discovered sample '{sample_name}' for resume")
    return result_files

def main():
    """Main function to run the pipeline."""
    args = parse_args()

    # Create output directory structure
    output_dir = args.outdir
    cleaned_dir = os.path.join(output_dir, "cleaned_seqs")

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cleaned_dir, exist_ok=True)

    # Create temp directory
    temp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(temp_dir, exist_ok=True)

    # Calculate total steps
    total_steps = 9
    current_step = 0

    # Print pipeline header
    logger.info("")
    logger.info("=" * 70)
    logger.info("VICAST-ANALYZE: Viral Genomics Analysis Pipeline")
    logger.info(f"Version: {__version__}")
    logger.info(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info("=" * 70)
    logger.info("")

    try:
        # Check for resume mode
        if args.resume_from_vcf:
            logger.info("=" * 70)
            logger.info("RESUME MODE: Skipping Steps 1-6 (QC workflow)")
            logger.info("Running Steps 7-9 only (Annotation workflow)")
            logger.info("=" * 70)
            logger.info("")

            # Discover existing VCF files from previous QC run
            variant_files = discover_existing_vcf_files(output_dir, args.r1)

            # Get accession for annotation
            accession = args.accession
            if not accession:
                raise ValueError("--accession is required when using --resume-from-vcf")

            # Extract sample name from R1 filename (for consistency with full pipeline)
            # Support both Illumina (_R1) and SRA (_1) naming conventions
            r1_base = os.path.basename(args.r1)
            sample_name = re.sub(r'(_R?[12])?(_001)?\.fastq\.gz$', '', r1_base)

            # Mark Steps 1-6 as skipped
            for step_num in range(1, 7):
                current_step += 1
                logger.info(f"⊘ STEP {current_step}/{total_steps}: (SKIPPED - Resume Mode)")

            logger.info("")
            logger.info("Proceeding directly to annotation workflow (Steps 7-9)...")
            logger.info("")

        # Normal workflow: Run Steps 1-6
        elif not args.resume_from_vcf:
            # Step 1: Prepare reference genome
            current_step += 1
            print_step_header(current_step, total_steps, "Prepare Reference Genome")

            reference_path = None
            if args.reference:
                reference_path = args.reference
                logger.info(f"Using provided reference genome: {reference_path}")
            elif args.accession:
                # Download reference genome if needed
                if not args.skip_download:
                    reference_path = download_reference_genome(
                        args.accession,
                        cleaned_dir,
                        args.force_download
                    )
                else:
                    # Use existing reference
                    reference_path = os.path.join(cleaned_dir, f"{args.accession}.fasta")
                    if not os.path.exists(reference_path):
                        raise FileNotFoundError(f"Reference genome not found: {reference_path}")
            else:
                raise ValueError("Either --accession or --reference must be provided")

            print_step_complete("Reference genome prepared", [reference_path] if reference_path else None)

            # Check if reference genome is in snpEff database
            accession = args.accession
            if not accession and reference_path:
                # Try to extract accession from reference file name
                accession = os.path.basename(reference_path).split('.')[0]

            if accession:
                logger.info("Checking snpEff database...")
                in_database = check_snpeff_database(accession, args.snpeff_jar, args.java_path)
                if not in_database:
                    logger.warning(f"Genome {accession} not found in snpEff database")
                    if args.add_to_snpeff:
                        logger.info(f"Attempting to add {accession} to snpEff database")
                        success = add_genome_to_snpeff(accession, reference_path, args.snpeff_jar, args.java_path)
                        if not success:
                            logger.error(f"Failed to add {accession} to snpEff database. Annotation will likely fail.")
                    else:
                        logger.warning("Use --add-to-snpeff to attempt adding the genome to snpEff")

            # Step 2: Calculate read statistics
            current_step += 1
            if not args.skip_stats:
                print_step_header(current_step, total_steps, "Calculate Read Statistics")
                stats_file = os.path.join(output_dir, "input_stats.txt")
                calculate_read_stats(args.r1, stats_file, args.threads, args.r2)
                print_step_complete("Read statistics calculated", [stats_file])
            else:
                logger.info(f"⊘ STEP {current_step}/{total_steps}: Calculate Read Statistics (SKIPPED)")
                logger.info("")

            # Step 3: Clean reads
            current_step += 1
            if not args.skip_qc:
                print_step_header(current_step, total_steps, "Clean Reads (Quality Control)")
                cleaned_files = clean_reads(cleaned_dir, args.r1, args.r2, args.threads)
                # Extract file paths from nested dict structure
                output_files = []
                for file_dict in cleaned_files.values():
                    if isinstance(file_dict, dict):
                        for filepath in file_dict.values():
                            if filepath and isinstance(filepath, str):
                                output_files.append(filepath)
                print_step_complete("Reads cleaned", output_files)
            else:
                logger.info(f"⊘ STEP {current_step}/{total_steps}: Clean Reads (SKIPPED)")
                logger.info("")
                cleaned_files = {}

            # Step 4: Map reads and call variants
            current_step += 1
            if not args.skip_mapping:
                print_step_header(current_step, total_steps, "Map Reads and Call Variants")
                variant_files = map_and_call_variants(
                    reference_path,
                    cleaned_dir,
                    args.threads,
                    cleaned_files,  # Pass the specific cleaned files
                    args.large_files,
                    args.extremely_large_files,
                    skip_vector_filter=getattr(args, 'skip_vector_filter', False)
                )
                # Extract file paths from nested dict structure
                output_files = []
                for file_dict in variant_files.values():
                    if isinstance(file_dict, dict):
                        for filepath in file_dict.values():
                            if filepath and isinstance(filepath, str) and os.path.exists(filepath):
                                output_files.append(filepath)
                print_step_complete("Read mapping and variant calling completed", output_files)
            else:
                logger.info(f"⊘ STEP {current_step}/{total_steps}: Map Reads and Call Variants (SKIPPED)")
                logger.info("")
                variant_files = {}

            # Step 5: Generate Depth File (QC) - only in normal mode
            current_step += 1
            print_step_header(current_step, total_steps, "Generate Depth File (QC)")

            # Extract sample name from R1 filename
            # Support both Illumina (_R1) and SRA (_1) naming conventions
            r1_base = os.path.basename(args.r1)
            sample_name = re.sub(r'(_R?[12])?(_001)?\.fastq\.gz$', '', r1_base)

            # Create results directory
            results_dir = os.path.join(os.getcwd(), f"{sample_name}_results")
            os.makedirs(results_dir, exist_ok=True)

            # Generate depth file from final BAM
            depth_file = os.path.join(results_dir, f"{sample_name}_depth.txt")
            logger.info(f"Generating per-position depth file: {depth_file}")

            # Get final BAM file from variant_files
            final_bam = None
            for sample, file_dict in variant_files.items():
                if isinstance(file_dict, dict) and 'final_bam' in file_dict:
                    final_bam = file_dict['final_bam']
                    break

            if final_bam and os.path.exists(final_bam):
                # Write header
                with open(depth_file, 'w') as f:
                    f.write("chrom\tposition\tdepth\n")

                # Run samtools depth
                depth_cmd = f"samtools depth {final_bam} >> {depth_file}"
                run_command(depth_cmd, shell=True)
                logger.info(f"Depth file created: {depth_file}")
                print_step_complete("Depth file generation completed", [depth_file])
            else:
                logger.warning("No final BAM file found, skipping depth file generation")
                logger.info("")

            # Step 6: Run Diagnostic Report (QC) - only in normal mode
            current_step += 1
            print_step_header(current_step, total_steps, "Generate Diagnostic Report (QC)")

            # Get script directory for diagnostic script
            script_dir = os.path.dirname(os.path.abspath(__file__))
            diagnostic_script = os.path.join(script_dir, "viral_diagnostic.sh")

            if os.path.exists(diagnostic_script):
                logger.info("Running comprehensive viral diagnostic analysis...")
                logger.info("This includes: mapping stats, de novo assembly, viral BLAST")

                diagnostic_cmd = f"bash {diagnostic_script} {args.r1} {args.r2} {accession} {sample_name} {args.threads}"
                if args.extremely_large_files:
                    diagnostic_cmd += " --extremely-large-files"

                try:
                    run_command(diagnostic_cmd, shell=True)

                    # List diagnostic outputs
                    diagnostic_dir = f"./diagnostic_{sample_name}"
                    diagnostic_outputs = [
                        os.path.join(diagnostic_dir, f"{sample_name}_diagnostic_report.txt"),
                        os.path.join(diagnostic_dir, f"diagnostic_{sample_name}_presentation_ready_report.html"),
                        os.path.join(diagnostic_dir, f"{sample_name}_top_hits.tsv"),
                        os.path.join(diagnostic_dir, f"{sample_name}_viral_blast.tsv")
                    ]

                    existing_outputs = [f for f in diagnostic_outputs if os.path.exists(f)]
                    print_step_complete("Diagnostic report generation completed", existing_outputs)
                except Exception as e:
                    logger.warning(f"Diagnostic report failed: {str(e)}")
                    logger.info("Pipeline will continue despite diagnostic failure")
                    logger.info("")
            else:
                logger.warning(f"Diagnostic script not found at: {diagnostic_script}")
                logger.info("Skipping diagnostic report generation")
                logger.info("")


        # Steps 7-9: Annotation workflow (runs in both normal and resume modes)
        # Step 7: Filter variants
        current_step += 1
        if not args.skip_variants:
            print_step_header(current_step, total_steps, "Filter Variants")
            variants_dir = os.path.join(cleaned_dir, "variants")
            filtered_files = filter_variants(variants_dir, variant_files, min_depth=args.min_depth, min_qual=args.min_qual)
            output_files = [os.path.join(variants_dir, v) for v in filtered_files.values() if v]
            print_step_complete("Variant filtering completed", output_files)
        else:
            logger.info(f"⊘ STEP {current_step}/{total_steps}: Filter Variants (SKIPPED)")
            logger.info("")
            filtered_files = {}

        # Step 8: Annotate variants
        current_step += 1
        if not args.skip_annotation and accession:
            print_step_header(current_step, total_steps, "Annotate Variants with snpEff")
            variants_dir = os.path.join(cleaned_dir, "variants")
            annotation_files = annotate_variants(variants_dir, accession, args.snpeff_jar, args.java_path, filtered_files, args.large_files, args.extremely_large_files)
            # Extract file paths from nested dict structure
            output_files = []
            for file_dict in annotation_files.values():
                if isinstance(file_dict, dict):
                    for filepath in file_dict.values():
                        if filepath and isinstance(filepath, str) and os.path.exists(filepath):
                            output_files.append(filepath)
            print_step_complete("Variant annotation completed", output_files)

            # Step 9: Parse annotations
            current_step += 1
            print_step_header(current_step, total_steps, "Parse Annotations")
            parsed_files = parse_annotations(variants_dir, args.min_depth, annotation_files)
            output_files = [v for v in parsed_files.values() if v and os.path.exists(v)]
            print_step_complete("Annotation parsing completed", output_files)
        else:
            logger.info(f"⊘ STEP {current_step}/{total_steps}: Annotate Variants (SKIPPED)")
            logger.info("")
            current_step += 1
            logger.info(f"⊘ STEP {current_step}/{total_steps}: Parse Annotations (SKIPPED)")
            logger.info("")

        # Pipeline completion summary
        logger.info("")
        logger.info("=" * 70)
        logger.info("✓ PIPELINE COMPLETED SUCCESSFULLY")
        logger.info(f"Finished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        logger.info("=" * 70)
        logger.info("")
        
    except Exception as e:
        logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)
    finally:
        # Cleanup
        if not args.keep_tmp and os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

if __name__ == "__main__":
    main()