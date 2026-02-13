#!/usr/bin/env python3
"""
VICAST Database Manager (Admin Only)

Compare local SnpEff databases to the published manifest, package new genomes,
and update the manifest for distribution.

Usage:
    # Compare local databases to manifest
    python manage_prebuilt_databases.py --compare --db-dir /path/to/snpeff_db

    # Package a new genome and add to manifest
    python manage_prebuilt_databases.py --package NC_002549.1 --db-dir /path/to/snpeff_db

    # Package with metadata
    python manage_prebuilt_databases.py --package NC_002549.1 --db-dir /path/to/snpeff_db \
        --name "Ebola-Zaire" --description "Zaire ebolavirus" \
        --tags filovirus,hemorrhagic,BSL4

    # List what's in the manifest
    python manage_prebuilt_databases.py --list
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tarfile
from datetime import date
from pathlib import Path


# Default paths
SCRIPT_DIR = Path(__file__).parent.resolve()
MANIFEST_FILE = SCRIPT_DIR / "manifest.json"
GITHUB_REPO = "shandley/VICAST"
GITHUB_BRANCH = "main"


def load_manifest():
    """Load the manifest.json file."""
    if not MANIFEST_FILE.exists():
        return {"version": "1.0", "last_updated": str(date.today()), "databases": []}
    with open(MANIFEST_FILE) as f:
        return json.load(f)


def save_manifest(manifest):
    """Save manifest.json with updated date."""
    manifest["last_updated"] = str(date.today())
    with open(MANIFEST_FILE, "w") as f:
        json.dump(manifest, f, indent=2)
        f.write("\n")
    print(f"  Updated: {MANIFEST_FILE}")


def get_manifest_accessions(manifest):
    """Return set of accessions in the manifest."""
    return {db["accession"] for db in manifest["databases"]}


def scan_local_databases(db_dir):
    """Scan a local SnpEff database directory for built genomes."""
    db_path = Path(db_dir)
    if not db_path.exists():
        print(f"Error: Database directory not found: {db_dir}")
        sys.exit(1)

    databases = {}
    for item in sorted(db_path.iterdir()):
        if item.is_dir():
            predictor = item / "snpEffectPredictor.bin"
            genes_gff = item / "genes.gff"
            sequences = item / "sequences.fa"

            if predictor.exists():
                # Count features in genes.gff
                feature_count = 0
                if genes_gff.exists():
                    with open(genes_gff) as f:
                        for line in f:
                            if not line.startswith("#") and line.strip():
                                feature_count += 1

                databases[item.name] = {
                    "path": str(item),
                    "has_predictor": True,
                    "has_genes": genes_gff.exists(),
                    "has_sequences": sequences.exists(),
                    "features": feature_count,
                    "size_bytes": sum(f.stat().st_size for f in item.iterdir() if f.is_file()),
                }

    return databases


def md5_file(filepath):
    """Calculate MD5 hash of a file."""
    md5 = hashlib.md5()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            md5.update(chunk)
    return md5.hexdigest()


def do_compare(db_dir):
    """Compare local databases to manifest."""
    manifest = load_manifest()
    manifest_accessions = get_manifest_accessions(manifest)
    local_dbs = scan_local_databases(db_dir)
    local_accessions = set(local_dbs.keys())

    # In manifest but not local
    missing_local = manifest_accessions - local_accessions
    # Local but not in manifest
    new_local = local_accessions - manifest_accessions
    # In both
    shared = manifest_accessions & local_accessions

    print("=" * 60)
    print("VICAST Database Comparison")
    print("=" * 60)
    print(f"\nManifest: {len(manifest_accessions)} databases")
    print(f"Local:    {len(local_accessions)} databases in {db_dir}")
    print()

    if new_local:
        print(f"NEW (local only, not in manifest) [{len(new_local)}]:")
        for acc in sorted(new_local):
            db = local_dbs[acc]
            size_mb = db["size_bytes"] / 1048576
            print(f"  + {acc:<25} {db['features']:>3} features  {size_mb:.2f} MB")
        print()
        print("  To publish these, run:")
        for acc in sorted(new_local):
            print(f"    python {Path(__file__).name} --package {acc} --db-dir {db_dir}")
        print()

    if missing_local:
        print(f"AVAILABLE (in manifest, not local) [{len(missing_local)}]:")
        manifest_lookup = {db["accession"]: db for db in manifest["databases"]}
        for acc in sorted(missing_local):
            db = manifest_lookup[acc]
            print(f"  - {acc:<25} {db.get('name', ''):<20} {db.get('description', '')}")
        print()

    if shared:
        print(f"INSTALLED (in both) [{len(shared)}]:")
        for acc in sorted(shared):
            print(f"  = {acc}")
        print()

    if not new_local:
        print("All local databases are already in the manifest.")


def do_package(accession, db_dir, name=None, description=None, tags=None):
    """Package a local database and add to manifest."""
    local_dbs = scan_local_databases(db_dir)

    if accession not in local_dbs:
        print(f"Error: {accession} not found in {db_dir}")
        print(f"Available: {', '.join(sorted(local_dbs.keys()))}")
        sys.exit(1)

    db_info = local_dbs[accession]
    output_dir = SCRIPT_DIR
    tar_file = output_dir / f"{accession}.tar.gz"
    md5_file_path = output_dir / f"{accession}.md5"

    print(f"Packaging {accession}...")
    print(f"  Source: {db_info['path']}")
    print(f"  Features: {db_info['features']}")

    # Create tar.gz
    with tarfile.open(tar_file, "w:gz") as tar:
        tar.add(db_info["path"], arcname=accession)

    # Calculate MD5
    md5 = md5_file(tar_file)
    with open(md5_file_path, "w") as f:
        f.write(f"{md5}  {accession}.tar.gz\n")

    size_bytes = tar_file.stat().st_size
    size_mb = round(size_bytes / 1048576, 2)

    print(f"  Created: {tar_file}")
    print(f"  Size: {size_mb} MB")
    print(f"  MD5: {md5}")

    # Update manifest
    manifest = load_manifest()

    # Check if already in manifest
    existing_idx = None
    for i, db in enumerate(manifest["databases"]):
        if db["accession"] == accession:
            existing_idx = i
            break

    entry = {
        "accession": accession,
        "name": name or accession,
        "description": description or f"Viral genome {accession}",
        "size_mb": size_mb,
        "features": db_info["features"],
        "tags": tags.split(",") if tags else [],
        "url": f"https://github.com/{GITHUB_REPO}/raw/{GITHUB_BRANCH}/prebuilt_databases/{accession}.tar.gz",
        "md5": md5,
    }

    if existing_idx is not None:
        # Preserve existing name/description/tags if not provided
        old = manifest["databases"][existing_idx]
        if not name:
            entry["name"] = old.get("name", accession)
        if not description:
            entry["description"] = old.get("description", "")
        if not tags:
            entry["tags"] = old.get("tags", [])
        manifest["databases"][existing_idx] = entry
        print(f"\n  Updated existing entry in manifest")
    else:
        manifest["databases"].append(entry)
        print(f"\n  Added new entry to manifest")

    save_manifest(manifest)

    print(f"\n{'='*60}")
    print(f"NEXT STEPS:")
    print(f"{'='*60}")
    if not name or name == accession:
        print(f"\n  1. Edit manifest.json to set proper name/description/tags:")
        print(f"     nano {MANIFEST_FILE}")
    print(f"\n  2. Commit and push:")
    print(f"     cd {SCRIPT_DIR.parent}")
    print(f"     git add prebuilt_databases/{accession}.tar.gz prebuilt_databases/{accession}.md5 prebuilt_databases/manifest.json")
    print(f"     git commit -m 'Add pre-built database: {accession}'")
    print(f"     git push")


def do_list():
    """List all databases in the manifest."""
    manifest = load_manifest()

    print("=" * 60)
    print(f"VICAST Pre-built Databases ({len(manifest['databases'])} total)")
    print(f"Last updated: {manifest.get('last_updated', 'unknown')}")
    print("=" * 60)

    # Group by tags
    print(f"\n  {'Accession':<25} {'Name':<20} {'Features':>8}  Tags")
    print(f"  {'-'*23}  {'-'*18}  {'-'*8}  {'-'*25}")

    for db in manifest["databases"]:
        tags_str = ", ".join(db.get("tags", []))
        print(f"  {db['accession']:<25} {db.get('name', ''):<20} {db.get('features', '?'):>8}  {tags_str}")


def main():
    parser = argparse.ArgumentParser(
        description="VICAST Database Manager (Admin Only)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Compare local databases to manifest
  %(prog)s --compare --db-dir /mnt/pathogen2/kathie/vicast_annotate/snpeff_db

  # Package and add a new genome to manifest
  %(prog)s --package NC_002549.1 --db-dir /mnt/pathogen2/kathie/vicast_annotate/snpeff_db \\
      --name "Ebola-Zaire" --description "Zaire ebolavirus" --tags filovirus,hemorrhagic

  # List all databases in manifest
  %(prog)s --list
        """,
    )

    parser.add_argument("--compare", action="store_true", help="Compare local databases to manifest")
    parser.add_argument("--package", metavar="ACCESSION", help="Package a genome and add to manifest")
    parser.add_argument("--list", action="store_true", help="List all databases in manifest")
    parser.add_argument("--db-dir", help="Path to local SnpEff database directory")
    parser.add_argument("--name", help="Short virus name (e.g., Ebola-Zaire)")
    parser.add_argument("--description", help="Full description")
    parser.add_argument("--tags", help="Comma-separated tags (e.g., filovirus,hemorrhagic)")

    args = parser.parse_args()

    if not any([args.compare, args.package, args.list]):
        parser.print_help()
        sys.exit(1)

    if args.list:
        do_list()

    if args.compare:
        if not args.db_dir:
            print("Error: --db-dir required for --compare")
            sys.exit(1)
        do_compare(args.db_dir)

    if args.package:
        if not args.db_dir:
            print("Error: --db-dir required for --package")
            sys.exit(1)
        do_package(args.package, args.db_dir, args.name, args.description, args.tags)


if __name__ == "__main__":
    main()
