#!/bin/bash

# =============================================================================
# VICAST BLAST Contamination Database Setup
# Downloads and installs the pre-built contamination screening database
# =============================================================================

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default locations
DEFAULT_DB_DIR="/opt/vicast/blast_db"
DB_NAME="microbial_contaminants"
RELEASE_URL="https://github.com/mihinduk/VICAST/releases/download/blast-db-v1.0/microbial_contaminants.tar.gz"

usage() {
    cat << EOF
Usage: $(basename "$0") [OPTIONS] [TARGET_DIR]

Download and install the VICAST BLAST contamination screening database.

The database contains 18,804 sequences:
  - NCBI RefSeq viral genomes (18,760 sequences)
  - Common lab contaminants:
    - Escherichia coli K12 MG1655
    - Pseudomonas aeruginosa PAO1
    - Staphylococcus aureus NCTC8325
    - Mycoplasma hyorhinis
    - Candida albicans SC5314
    - Cryptococcus neoformans JEC21
    - Saccharomyces cerevisiae S288C

Arguments:
    TARGET_DIR      Directory to install database (default: $DEFAULT_DB_DIR)

Options:
    --info          Show information about installed database
    --force         Overwrite existing database
    --from-file F   Install from local tar.gz file instead of downloading
    --help          Show this help message

Examples:
    # Install database (default location)
    $(basename "$0")

    # Install to custom directory
    $(basename "$0") /path/to/blast_db

    # Check installed database
    $(basename "$0") --info

    # Install from local file (no download)
    $(basename "$0") --from-file /tmp/microbial_contaminants.tar.gz

EOF
    exit 0
}

# Show database info
show_info() {
    local db_dir="${1:-$DEFAULT_DB_DIR}"
    local db_path="${db_dir}/${DB_NAME}"

    echo -e "${BLUE}=== VICAST BLAST Contamination Database ===${NC}"
    echo ""

    if [[ -f "${db_path}.nhr" ]] || [[ -f "${db_path}.00.nhr" ]]; then
        echo -e "${GREEN}Status: INSTALLED${NC}"
        echo "Location: ${db_path}"
        echo ""

        # Show database stats
        if command -v blastdbcmd &> /dev/null; then
            blastdbcmd -db "${db_path}" -info 2>/dev/null || echo "  (blastdbcmd info unavailable)"
        fi

        # Show info file if present
        if [[ -f "${db_dir}/${DB_NAME}_info.txt" ]]; then
            echo ""
            cat "${db_dir}/${DB_NAME}_info.txt"
        fi

        # Show file sizes
        echo ""
        echo "Database files:"
        ls -lh "${db_path}".* 2>/dev/null | awk '{printf "  %-40s %s\n", $NF, $5}'

        # Show alias status
        echo ""
        if [[ -f "${db_dir}/vicast_combined.nal" ]]; then
            echo -e "${GREEN}Alias: vicast_combined (active)${NC}"
        else
            echo -e "${YELLOW}Alias: not created (run setup_blast_db.sh to create)${NC}"
        fi

        # Show user extensions
        if [[ -d "${db_dir}/user_extensions" ]]; then
            local ext_count=0
            for ext in "${db_dir}"/user_extensions/*.nhr; do
                [[ -f "$ext" ]] && ext_count=$((ext_count + 1))
            done
            if [[ $ext_count -gt 0 ]]; then
                echo "User extensions: ${ext_count} database(s)"
                for ext in "${db_dir}"/user_extensions/*.nhr; do
                    [[ -f "$ext" ]] || continue
                    local ext_name
                    ext_name=$(basename "$ext" .nhr)
                    local seq_count
                    seq_count=$(blastdbcmd -db "${db_dir}/user_extensions/${ext_name}" -info 2>/dev/null | grep -o '[0-9,]* sequences' | head -1 || echo "? sequences")
                    echo "  - ${ext_name}: ${seq_count}"
                done
            fi
        fi
    else
        echo -e "${YELLOW}Status: NOT INSTALLED${NC}"
        echo ""
        echo "Install with: $(basename "$0")"
    fi
}

# Create or update the vicast_combined alias DB
# This alias wraps the base DB (and any user extensions) so all consumers
# use a single consistent path: vicast_combined
create_alias() {
    local db_dir="$1"
    local alias_name="vicast_combined"

    if ! command -v blastdb_aliastool &> /dev/null; then
        echo -e "${YELLOW}blastdb_aliastool not found, skipping alias creation${NC}"
        return 0
    fi

    # Build list of databases: base DB + any user extensions
    local db_list="${DB_NAME}"
    if [[ -d "${db_dir}/user_extensions" ]]; then
        for ext_db in "${db_dir}"/user_extensions/*.nhr; do
            if [[ -f "$ext_db" ]]; then
                local ext_name
                ext_name=$(basename "$ext_db" .nhr)
                db_list="${db_list} user_extensions/${ext_name}"
            fi
        done
    fi

    echo "Creating BLAST alias database: ${alias_name} -> [${db_list}]"
    (cd "$db_dir" && blastdb_aliastool -dblist "${db_list}" \
        -dbtype nucl \
        -out "${alias_name}" \
        -title "VICAST Combined BLAST DB") 2>/dev/null || {
        echo -e "${YELLOW}Warning: Could not create alias DB${NC}"
        return 0
    }
    echo -e "${GREEN}Alias database created: ${db_dir}/${alias_name}${NC}"
}

# Download and install database
install_db() {
    local db_dir="$1"
    local force="$2"
    local from_file="${3:-}"

    # Check if already installed
    if [[ -f "${db_dir}/${DB_NAME}.nhr" ]] || [[ -f "${db_dir}/${DB_NAME}.00.nhr" ]]; then
        if [[ "$force" != "true" ]]; then
            echo -e "${YELLOW}Database already installed at ${db_dir}${NC}"
            echo "Use --force to reinstall"
            return 0
        fi
        echo "Removing existing database..."
        rm -f "${db_dir}/${DB_NAME}".*
    fi

    mkdir -p "$db_dir"

    local temp_file=""

    if [[ -n "$from_file" ]]; then
        # Install from explicit local file
        if [[ ! -f "$from_file" ]]; then
            echo -e "${RED}Error: File not found: $from_file${NC}" >&2
            return 1
        fi
        echo "Installing from local file: $from_file"
        temp_file="$from_file"
    elif [[ -f "${VICAST_HOME:-/opt/vicast}/microbial_contaminants.tar.gz" ]]; then
        # Found local copy in VICAST directory (bundled in Docker build context)
        echo "Found local database archive in VICAST directory"
        temp_file="${VICAST_HOME:-/opt/vicast}/microbial_contaminants.tar.gz"
        # Don't set from_file — allow cleanup after extraction to save image space
    else
        # Download from GitHub Release
        temp_file="/tmp/microbial_contaminants_$$.tar.gz"
        echo "Downloading contamination database..."
        echo "URL: $RELEASE_URL"

        if command -v curl &> /dev/null; then
            curl -fL --progress-bar "$RELEASE_URL" -o "$temp_file" || {
                echo -e "${RED}Download failed${NC}" >&2
                echo "Try: $(basename "$0") --from-file /path/to/microbial_contaminants.tar.gz" >&2
                rm -f "$temp_file"
                return 1
            }
        elif command -v wget &> /dev/null; then
            wget --show-progress -q "$RELEASE_URL" -O "$temp_file" || {
                echo -e "${RED}Download failed${NC}" >&2
                rm -f "$temp_file"
                return 1
            }
        else
            echo -e "${RED}Error: curl or wget required${NC}" >&2
            return 1
        fi
    fi

    # Extract
    echo "Extracting to ${db_dir}..."
    tar -xzf "$temp_file" -C "$db_dir" || {
        echo -e "${RED}Extraction failed${NC}" >&2
        [[ -z "$from_file" ]] && rm -f "$temp_file"
        return 1
    }

    # Clean up temp file (only if we downloaded it)
    [[ -z "$from_file" ]] && rm -f "$temp_file"

    # Verify installation
    if [[ -f "${db_dir}/${DB_NAME}.nhr" ]] || [[ -f "${db_dir}/${DB_NAME}.00.nhr" ]]; then
        echo -e "${GREEN}Database installed successfully${NC}"
        echo "Location: ${db_dir}/${DB_NAME}"

        # Create alias DB so all consumers use a consistent name (vicast_combined)
        # This enables user extensions via extend_blast_db.sh
        create_alias "$db_dir"

        # Show stats
        if command -v blastdbcmd &> /dev/null; then
            echo ""
            blastdbcmd -db "${db_dir}/${DB_NAME}" -info 2>/dev/null || true
        fi
    else
        echo -e "${RED}Error: Database files not found after extraction${NC}" >&2
        return 1
    fi
}

# Main
main() {
    local db_dir="$DEFAULT_DB_DIR"
    local force=false
    local from_file=""
    local action="install"

    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h)
                usage
                ;;
            --info)
                action="info"
                shift
                ;;
            --force)
                force=true
                shift
                ;;
            --from-file)
                from_file="$2"
                shift 2
                ;;
            -*)
                echo -e "${RED}Unknown option: $1${NC}" >&2
                usage
                ;;
            *)
                db_dir="$1"
                shift
                ;;
        esac
    done

    case "$action" in
        info)
            show_info "$db_dir"
            ;;
        install)
            install_db "$db_dir" "$force" "$from_file"
            ;;
    esac
}

main "$@"
