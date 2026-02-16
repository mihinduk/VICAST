#!/bin/bash

# =============================================================================
# VICAST Pre-built Database Installer
# Install pre-built SnpEff viral genome databases
# =============================================================================

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default manifest URL (can be overridden)
MANIFEST_URL="${VICAST_MANIFEST_URL:-https://raw.githubusercontent.com/mihinduk/VICAST/main/prebuilt_databases/manifest.json}"
SNPEFF_DATA="${SNPEFF_DATA:-${SNPEFF_HOME}/data}"

# Usage message
usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] [ACCESSION]

Install pre-built SnpEff viral genome databases for VICAST.

Commands:
    --list              List all available pre-built databases
    --installed         Show currently installed databases
    --install ACCESSION Install specific database (e.g., NC_001474.2)
    --remove ACCESSION  Remove installed database
    --help              Show this help message

Options:
    --force             Overwrite existing database
    --tag TAG           Install all databases with specific tag (e.g., dengue)

Examples:
    # List available databases
    $(basename "$0") --list

    # Install DENV-2 database
    $(basename "$0") --install NC_001474.2

    # Install all dengue databases
    $(basename "$0") --tag dengue

    # Force reinstall (overwrite custom annotations)
    $(basename "$0") --install NC_001474.2 --force

Inspecting Databases:
    # After installation, inspect the annotation
    snpeff dump NC_001474.2 | head
    cat \$SNPEFF_DATA/NC_001474.2/genes.gff

    # Don't like it? Build your own and overwrite
    step1_parse_viral_genome.py NC_001474.2
    step2_add_to_snpeff.py NC_001474.2 NC_001474.2_features.tsv --force

Environment Variables:
    SNPEFF_DATA         SnpEff data directory (default: \$SNPEFF_HOME/data)
    VICAST_MANIFEST_URL Custom manifest URL

EOF
    exit 0
}

# Check if SnpEff is configured
check_snpeff() {
    if [[ -z "${SNPEFF_HOME:-}" ]]; then
        echo -e "${RED}Error: SNPEFF_HOME not set${NC}" >&2
        echo "Please set SNPEFF_HOME environment variable" >&2
        exit 1
    fi

    if [[ ! -d "$SNPEFF_DATA" ]]; then
        echo -e "${YELLOW}Creating SnpEff data directory: $SNPEFF_DATA${NC}"
        mkdir -p "$SNPEFF_DATA"
    fi
}

# Download manifest
get_manifest() {
    local temp_manifest="/tmp/vicast_manifest_$$.json"

    if command -v curl &> /dev/null; then
        curl -fsSL "$MANIFEST_URL" -o "$temp_manifest" 2>/dev/null || {
            echo -e "${RED}Error: Failed to download manifest${NC}" >&2
            echo "URL: $MANIFEST_URL" >&2
            return 1
        }
    elif command -v wget &> /dev/null; then
        wget -q "$MANIFEST_URL" -O "$temp_manifest" 2>/dev/null || {
            echo -e "${RED}Error: Failed to download manifest${NC}" >&2
            return 1
        }
    else
        echo -e "${RED}Error: curl or wget required${NC}" >&2
        return 1
    fi

    echo "$temp_manifest"
}

# List available databases
list_databases() {
    echo -e "${BLUE}=== VICAST Pre-built Databases ===${NC}\n"

    local manifest
    manifest=$(get_manifest) || exit 1

    # Use python for rich table display with new metadata fields
    python3 -c "
import json, sys

with open('$manifest') as f:
    data = json.load(f)

total = len(data['databases'])
print(f'  {total} genomes available (manifest v{data.get(\"version\", \"1.0\")})')
print(f'  Last updated: {data.get(\"last_updated\", \"unknown\")}')
print()

# Column headers
fmt = '  {:<18} {:<40} {:<14} {:<8} {:<22} {}'
print(fmt.format('ACCESSION', 'VIRUS NAME', 'ABBREVIATION', 'BALT.', 'GENOME TYPE', 'NOTES'))
print('  ' + '-' * 120)

# Group by Baltimore classification for organized display
from collections import OrderedDict
groups = OrderedDict()
for db in data['databases']:
    balt = db.get('baltimore', '?')
    key = {'II': 'Group II  - ssDNA',
           'III': 'Group III - dsRNA',
           'IV': 'Group IV  - (+)ssRNA',
           'V': 'Group V   - (-)ssRNA'}.get(balt, f'Group {balt}')
    groups.setdefault(key, []).append(db)

for group_name, dbs in groups.items():
    print(f'\n  {group_name}')
    print('  ' + '~' * 120)
    for db in dbs:
        acc = db['accession']
        virus = db.get('virus_name', db.get('description', ''))[:40]
        abbr = db.get('name', '')[:14]
        balt = db.get('baltimore', '?')
        gtype = db.get('genome_type', '')[:22]
        notes = db.get('notes', '')
        print(fmt.format(acc, virus, abbr, balt, gtype, notes))
        seg_accs = db.get('segment_accessions', [])
        if seg_accs:
            joined = ', '.join(seg_accs)
            print(f'                      NCBI segments: {joined}')

print()
print('  Install a database:  install_prebuilt_database.sh --install <ACCESSION>')
print('  Show installed:      install_prebuilt_database.sh --installed')
print()
"

    rm -f "$manifest"
}

# Show installed databases
show_installed() {
    echo -e "${BLUE}=== Installed Databases ===${NC}\n"

    if [[ ! -d "$SNPEFF_DATA" ]]; then
        echo "No databases installed (SnpEff data directory doesn't exist)"
        return
    fi

    local count=0
    for genome_dir in "$SNPEFF_DATA"/*; do
        if [[ -d "$genome_dir" ]] && [[ -f "$genome_dir/snpEffectPredictor.bin" ]]; then
            local genome=$(basename "$genome_dir")
            local size=$(du -sh "$genome_dir" 2>/dev/null | cut -f1)
            printf "  %-20s %10s\n" "$genome" "$size"
            ((count++))
        fi
    done

    if [[ $count -eq 0 ]]; then
        echo "  No databases installed"
    else
        echo ""
        echo "Total: $count database(s)"
    fi
}

# Install database
install_database() {
    local accession="$1"
    local force="${2:-false}"

    echo -e "${GREEN}Installing pre-built database: $accession${NC}"

    # Get manifest
    local manifest
    manifest=$(get_manifest) || exit 1

    # Find database in manifest
    local db_info
    if command -v jq &> /dev/null; then
        db_info=$(jq -r ".databases[] | select(.accession==\"$accession\")" "$manifest")
        if [[ -z "$db_info" ]]; then
            echo -e "${RED}Error: Database $accession not found in manifest${NC}" >&2
            rm -f "$manifest"
            return 1
        fi

        local url=$(echo "$db_info" | jq -r '.url')
        local md5=$(echo "$db_info" | jq -r '.md5')
        local name=$(echo "$db_info" | jq -r '.name')
    else
        # Python fallback
        local result=$(python3 -c "
import json, sys
with open('$manifest') as f:
    data = json.load(f)
    for db in data['databases']:
        if db['accession'] == '$accession':
            print(f'{db[\"url\"]}|{db[\"md5\"]}|{db[\"name\"]}')
            sys.exit(0)
    sys.exit(1)
" 2>/dev/null)

        if [[ $? -ne 0 ]]; then
            echo -e "${RED}Error: Database $accession not found${NC}" >&2
            rm -f "$manifest"
            return 1
        fi

        IFS='|' read -r url md5 name <<< "$result"
    fi

    rm -f "$manifest"

    # Check if already installed
    if [[ -d "$SNPEFF_DATA/$accession" ]] && [[ "$force" != "true" ]]; then
        echo -e "${YELLOW}Database $accession already installed${NC}"
        echo "Use --force to reinstall"
        return 0
    fi

    # Download database
    local temp_file="/tmp/${accession}.tar.gz"
    echo "Downloading $name ($accession)..."

    if command -v curl &> /dev/null; then
        curl -fL --progress-bar "$url" -o "$temp_file" || {
            echo -e "${RED}Download failed${NC}" >&2
            return 1
        }
    else
        wget --show-progress -q "$url" -O "$temp_file" || {
            echo -e "${RED}Download failed${NC}" >&2
            return 1
        }
    fi

    # Verify MD5 if available
    if [[ -n "$md5" ]] && [[ "$md5" != "null" ]]; then
        echo "Verifying checksum..."
        local actual_md5
        if command -v md5sum &> /dev/null; then
            actual_md5=$(md5sum "$temp_file" | cut -d' ' -f1)
        elif command -v md5 &> /dev/null; then
            actual_md5=$(md5 -q "$temp_file")
        fi

        if [[ -n "$actual_md5" ]] && [[ "$actual_md5" != "$md5" ]]; then
            echo -e "${RED}Checksum mismatch!${NC}" >&2
            echo "Expected: $md5" >&2
            echo "Got: $actual_md5" >&2
            rm -f "$temp_file"
            return 1
        fi
    fi

    # Extract database
    echo "Installing to $SNPEFF_DATA/$accession..."
    mkdir -p "$SNPEFF_DATA"
    tar -xzf "$temp_file" -C "$SNPEFF_DATA" || {
        echo -e "${RED}Extraction failed${NC}" >&2
        rm -f "$temp_file"
        return 1
    }

    rm -f "$temp_file"

    # Add genome to SnpEff config if not already there
    # Check custom config first (mounted volume), then built-in config
    local snpeff_config=""
    local custom_config="${SNPEFF_DATA}/snpEff.config"
    local builtin_config="${SNPEFF_HOME}/snpEff.config"

    if [[ -f "$custom_config" ]]; then
        snpeff_config="$custom_config"
    elif [[ -f "${SNPEFF_CONFIG_CUSTOM:-}" ]]; then
        snpeff_config="$SNPEFF_CONFIG_CUSTOM"
    elif [[ -f "$builtin_config" ]]; then
        snpeff_config="$builtin_config"
    fi

    if [[ -n "$snpeff_config" ]]; then
        if ! grep -q "^${accession}.genome" "$snpeff_config"; then
            echo "Adding $accession to SnpEff config: $snpeff_config"
            echo "" >> "$snpeff_config"
            echo "# ${accession} - Installed by install_prebuilt_database.sh" >> "$snpeff_config"
            echo "${accession}.genome : ${accession}" >> "$snpeff_config"
        else
            echo "Config entry for $accession already exists in $snpeff_config"
        fi
    else
        echo -e "${YELLOW}Warning: No SnpEff config file found to add genome entry${NC}"
        echo "You may need to manually add: ${accession}.genome : ${accession}"
    fi

    echo -e "${GREEN}✓ Successfully installed $accession${NC}"
    echo ""
    echo "Test with: snpeff dump $accession | head"
}

# Remove database
remove_database() {
    local accession="$1"

    if [[ ! -d "$SNPEFF_DATA/$accession" ]]; then
        echo -e "${YELLOW}Database $accession not installed${NC}"
        return 0
    fi

    echo -e "${YELLOW}Removing database: $accession${NC}"
    rm -rf "$SNPEFF_DATA/$accession"
    echo -e "${GREEN}✓ Removed${NC}"
}

# Install databases by tag
install_by_tag() {
    local tag="$1"
    local force="${2:-false}"

    echo -e "${GREEN}Installing all databases tagged: $tag${NC}\n"

    local manifest
    manifest=$(get_manifest) || exit 1

    # Get accessions with matching tag
    local accessions
    if command -v jq &> /dev/null; then
        accessions=$(jq -r ".databases[] | select(.tags[]? == \"$tag\") | .accession" "$manifest")
    else
        accessions=$(python3 -c "
import json
with open('$manifest') as f:
    data = json.load(f)
    for db in data['databases']:
        if '$tag' in db.get('tags', []):
            print(db['accession'])
")
    fi

    rm -f "$manifest"

    if [[ -z "$accessions" ]]; then
        echo -e "${YELLOW}No databases found with tag: $tag${NC}"
        return 1
    fi

    # Install each
    local count=0
    while IFS= read -r acc; do
        install_database "$acc" "$force"
        ((count++))
        echo ""
    done <<< "$accessions"

    echo -e "${GREEN}Installed $count database(s)${NC}"
}

# Main
main() {
    if [[ $# -eq 0 ]]; then
        usage
    fi

    local command=""
    local accession=""
    local force=false
    local tag=""

    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h)
                usage
                ;;
            --list)
                command="list"
                shift
                ;;
            --installed)
                command="installed"
                shift
                ;;
            --install)
                command="install"
                accession="$2"
                shift 2
                ;;
            --remove)
                command="remove"
                accession="$2"
                shift 2
                ;;
            --tag)
                tag="$2"
                shift 2
                ;;
            --force)
                force=true
                shift
                ;;
            *)
                echo -e "${RED}Unknown option: $1${NC}" >&2
                usage
                ;;
        esac
    done

    check_snpeff

    case "$command" in
        list)
            list_databases
            ;;
        installed)
            show_installed
            ;;
        install)
            if [[ -z "$accession" ]]; then
                echo -e "${RED}Error: Accession required${NC}" >&2
                usage
            fi
            install_database "$accession" "$force"
            ;;
        remove)
            if [[ -z "$accession" ]]; then
                echo -e "${RED}Error: Accession required${NC}" >&2
                usage
            fi
            remove_database "$accession"
            ;;
        *)
            if [[ -n "$tag" ]]; then
                install_by_tag "$tag" "$force"
            else
                usage
            fi
            ;;
    esac
}

main "$@"
