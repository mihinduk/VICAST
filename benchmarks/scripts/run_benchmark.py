#!/usr/bin/env python3
"""
VICAST Benchmark Runner

Compares VICAST annotation against VADR and/or reference annotations.

Usage:
    python run_benchmark.py --genome NC_001477 --tools vicast,vadr
    python run_benchmark.py --genome-list genomes.txt --output results/
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Optional, Dict, Tuple
from datetime import datetime


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run."""
    genome_id: str
    tool: str
    runtime_seconds: float
    num_genes: int
    num_cds: int
    num_features: int
    success: bool
    error_message: Optional[str] = None
    output_file: Optional[str] = None


@dataclass
class ComparisonResult:
    """Comparison between two annotation sets."""
    genome_id: str
    tool1: str
    tool2: str
    shared_genes: int
    tool1_only: int
    tool2_only: int
    boundary_matches: int
    boundary_mismatches: int
    jaccard_similarity: float


def check_tool_available(tool: str) -> bool:
    """Check if a tool is available on the system."""
    if tool == "vicast":
        # Check for vicast scripts
        vicast_script = Path(__file__).parent.parent.parent / "vicast-annotate" / "vicast_annotate.py"
        return vicast_script.exists()
    elif tool == "vadr":
        # Check for v-annotate.pl
        try:
            result = subprocess.run(
                ["which", "v-annotate.pl"],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except Exception:
            return False
    return False


def run_vicast(genome_id: str, output_dir: Path) -> BenchmarkResult:
    """Run VICAST annotation on a genome."""
    start_time = time.time()

    vicast_script = Path(__file__).parent.parent.parent / "vicast-annotate" / "vicast_annotate.py"

    try:
        result = subprocess.run(
            [sys.executable, str(vicast_script), genome_id, "--auto", "--output-dir", str(output_dir)],
            capture_output=True,
            text=True,
            timeout=300  # 5 minute timeout
        )

        runtime = time.time() - start_time
        success = result.returncode == 0

        # Parse output to get counts
        num_genes, num_cds, num_features = 0, 0, 0
        gff_file = output_dir / f"{genome_id}.gff3"

        if gff_file.exists():
            with open(gff_file) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        feature_type = parts[2]
                        num_features += 1
                        if feature_type == "gene":
                            num_genes += 1
                        elif feature_type == "CDS":
                            num_cds += 1

        return BenchmarkResult(
            genome_id=genome_id,
            tool="vicast",
            runtime_seconds=runtime,
            num_genes=num_genes,
            num_cds=num_cds,
            num_features=num_features,
            success=success,
            error_message=result.stderr if not success else None,
            output_file=str(gff_file) if gff_file.exists() else None
        )

    except subprocess.TimeoutExpired:
        return BenchmarkResult(
            genome_id=genome_id,
            tool="vicast",
            runtime_seconds=300,
            num_genes=0,
            num_cds=0,
            num_features=0,
            success=False,
            error_message="Timeout after 300 seconds"
        )
    except Exception as e:
        return BenchmarkResult(
            genome_id=genome_id,
            tool="vicast",
            runtime_seconds=time.time() - start_time,
            num_genes=0,
            num_cds=0,
            num_features=0,
            success=False,
            error_message=str(e)
        )


def run_vadr(genome_id: str, fasta_file: Path, output_dir: Path, model_dir: Optional[Path] = None) -> BenchmarkResult:
    """Run VADR annotation on a genome."""
    start_time = time.time()

    try:
        cmd = ["v-annotate.pl"]
        if model_dir:
            cmd.extend(["--mdir", str(model_dir)])
        cmd.extend([str(fasta_file), str(output_dir)])

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300
        )

        runtime = time.time() - start_time
        success = result.returncode == 0

        # Parse VADR output
        num_genes, num_cds, num_features = 0, 0, 0
        gff_files = list(output_dir.glob("*.vadr.pass.tbl"))

        for tbl_file in gff_files:
            with open(tbl_file) as f:
                for line in f:
                    if line.startswith("\t\t\tgene"):
                        num_genes += 1
                    elif line.startswith("\t\t\tCDS"):
                        num_cds += 1
                    num_features += 1

        return BenchmarkResult(
            genome_id=genome_id,
            tool="vadr",
            runtime_seconds=runtime,
            num_genes=num_genes,
            num_cds=num_cds,
            num_features=num_features,
            success=success,
            error_message=result.stderr if not success else None,
            output_file=str(gff_files[0]) if gff_files else None
        )

    except subprocess.TimeoutExpired:
        return BenchmarkResult(
            genome_id=genome_id,
            tool="vadr",
            runtime_seconds=300,
            num_genes=0,
            num_cds=0,
            num_features=0,
            success=False,
            error_message="Timeout after 300 seconds"
        )
    except FileNotFoundError:
        return BenchmarkResult(
            genome_id=genome_id,
            tool="vadr",
            runtime_seconds=0,
            num_genes=0,
            num_cds=0,
            num_features=0,
            success=False,
            error_message="VADR (v-annotate.pl) not found in PATH"
        )
    except Exception as e:
        return BenchmarkResult(
            genome_id=genome_id,
            tool="vadr",
            runtime_seconds=time.time() - start_time,
            num_genes=0,
            num_cds=0,
            num_features=0,
            success=False,
            error_message=str(e)
        )


def download_genome(genome_id: str, output_dir: Path) -> Optional[Path]:
    """Download genome FASTA from NCBI."""
    from Bio import Entrez, SeqIO

    Entrez.email = os.environ.get("NCBI_EMAIL", "benchmark@example.com")

    try:
        handle = Entrez.efetch(db="nucleotide", id=genome_id, rettype="fasta", retmode="text")
        fasta_file = output_dir / f"{genome_id}.fasta"
        with open(fasta_file, "w") as f:
            f.write(handle.read())
        handle.close()
        return fasta_file
    except Exception as e:
        print(f"Failed to download {genome_id}: {e}")
        return None


def generate_report(results: List[BenchmarkResult], output_file: Path):
    """Generate benchmark report."""
    report = {
        "timestamp": datetime.now().isoformat(),
        "num_genomes": len(set(r.genome_id for r in results)),
        "tools_tested": list(set(r.tool for r in results)),
        "results": [asdict(r) for r in results],
        "summary": {
            "vicast": {
                "success_rate": sum(1 for r in results if r.tool == "vicast" and r.success) / max(1, sum(1 for r in results if r.tool == "vicast")),
                "avg_runtime": sum(r.runtime_seconds for r in results if r.tool == "vicast") / max(1, sum(1 for r in results if r.tool == "vicast")),
                "avg_cds": sum(r.num_cds for r in results if r.tool == "vicast" and r.success) / max(1, sum(1 for r in results if r.tool == "vicast" and r.success))
            },
            "vadr": {
                "success_rate": sum(1 for r in results if r.tool == "vadr" and r.success) / max(1, sum(1 for r in results if r.tool == "vadr")),
                "avg_runtime": sum(r.runtime_seconds for r in results if r.tool == "vadr") / max(1, sum(1 for r in results if r.tool == "vadr")),
                "avg_cds": sum(r.num_cds for r in results if r.tool == "vadr" and r.success) / max(1, sum(1 for r in results if r.tool == "vadr" and r.success))
            }
        }
    }

    with open(output_file, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nBenchmark Report: {output_file}")
    print(f"  Genomes tested: {report['num_genomes']}")
    print(f"  Tools: {', '.join(report['tools_tested'])}")

    for tool in ["vicast", "vadr"]:
        if tool in report["summary"]:
            s = report["summary"][tool]
            print(f"\n  {tool.upper()}:")
            print(f"    Success rate: {s['success_rate']*100:.1f}%")
            print(f"    Avg runtime: {s['avg_runtime']:.2f}s")
            print(f"    Avg CDS count: {s['avg_cds']:.1f}")


def main():
    parser = argparse.ArgumentParser(description="VICAST Benchmark Runner")
    parser.add_argument("--genome", help="Single genome ID to benchmark")
    parser.add_argument("--genome-list", help="File with genome IDs (one per line)")
    parser.add_argument("--tools", default="vicast", help="Comma-separated tools to test (vicast,vadr)")
    parser.add_argument("--output", default="benchmarks/results", help="Output directory")
    parser.add_argument("--vadr-models", help="VADR model directory")

    args = parser.parse_args()

    # Parse tools
    tools = [t.strip().lower() for t in args.tools.split(",")]

    # Check tool availability
    for tool in tools:
        if not check_tool_available(tool):
            print(f"Warning: {tool} not available on this system")

    # Get genome list
    genomes = []
    if args.genome:
        genomes = [args.genome]
    elif args.genome_list:
        with open(args.genome_list) as f:
            genomes = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    else:
        parser.error("Must specify --genome or --genome-list")

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run benchmarks
    results = []

    for genome_id in genomes:
        print(f"\nBenchmarking {genome_id}...")

        # Create temp directory for this genome
        genome_dir = output_dir / genome_id
        genome_dir.mkdir(exist_ok=True)

        # Download genome if needed
        fasta_file = genome_dir / f"{genome_id}.fasta"
        if not fasta_file.exists():
            print(f"  Downloading {genome_id}...")
            fasta_file = download_genome(genome_id, genome_dir)
            if not fasta_file:
                continue

        # Run each tool
        if "vicast" in tools and check_tool_available("vicast"):
            print(f"  Running VICAST...")
            result = run_vicast(genome_id, genome_dir)
            results.append(result)
            print(f"    {'✓' if result.success else '✗'} {result.runtime_seconds:.2f}s, {result.num_cds} CDS")

        if "vadr" in tools and check_tool_available("vadr"):
            print(f"  Running VADR...")
            vadr_dir = genome_dir / "vadr_output"
            vadr_dir.mkdir(exist_ok=True)
            result = run_vadr(
                genome_id,
                fasta_file,
                vadr_dir,
                Path(args.vadr_models) if args.vadr_models else None
            )
            results.append(result)
            print(f"    {'✓' if result.success else '✗'} {result.runtime_seconds:.2f}s, {result.num_cds} CDS")

    # Generate report
    report_file = output_dir / "benchmark_report.json"
    generate_report(results, report_file)


if __name__ == "__main__":
    main()
