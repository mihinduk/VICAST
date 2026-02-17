#!/bin/bash
# =============================================================================
# Extend the VICAST BLAST contamination database with custom sequences
#
# Adds user-provided FASTA files on top of the base database using BLAST
# alias databases (blastdb_aliastool). No re-download of the base DB needed.
#
# Usage:
#   extend_blast_db.sh [OPTIONS] FASTA [FASTA ...]
#   extend_blast_db.sh --list
#   extend_blast_db.sh --reset
#
# Docker example:
#   docker run --rm --user $(id -u):$(id -g) \
#       -v $(pwd):/data -v blast_db:/opt/vicast/blast_db \
#       vicast:latest extend_blast_db.sh /data/custom_contaminants.fasta
# =============================================================================

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

DB_DIR="${BLAST_DB_DIR:-/opt/vicast/blast_db}"
BASE_DB_NAME="microbial_contaminants"
ALIAS_NAME="vicast_combined"
EXT_DIR="user_extensions"
EXT_NAME="user_custom"

usage() {
    cat << 'EOF'
Usage: extend_blast_db.sh [OPTIONS] FASTA [FASTA ...]

Extend the VICAST BLAST contamination database with custom sequences.
Uses BLAST alias databases for zero-overhead extension of the base DB.

Arguments:
    FASTA               One or more FASTA files to add to the database

Options:
    --name NAME         Name for extension database (default: user_custom)
    --db-dir DIR        BLAST database directory (default: $BLAST_DB_DIR)
    --list              Show current database composition
    --reset             Remove all user extensions, restore base DB only
    --help              Show this help message

Examples:
    # Add custom lab strain sequences
    extend_blast_db.sh my_lab_strains.fasta

    # Add with a named panel
    extend_blast_db.sh --name mycoplasma_panel extra_mycoplasma.fasta

    # Check database composition
    extend_blast_db.sh --list

    # Remove all extensions
    extend_blast_db.sh --reset

EOF
    exit 0
}

# Show current database composition
list_db() {
    echo -e "${BLUE}=== VICAST BLAST Database Composition ===${NC}"
    echo ""

    # Check base DB
    if [[ -f "${DB_DIR}/${BASE_DB_NAME}.nhr" ]] || [[ -f "${DB_DIR}/${BASE_DB_NAME}.00.nhr" ]]; then
        local base_info
        base_info=$(blastdbcmd -db "${DB_DIR}/${BASE_DB_NAME}" -info 2>/dev/null | head -3 || echo "  (info unavailable)")
        echo -e "${GREEN}Base database: ${BASE_DB_NAME}${NC}"
        echo "$base_info"
    else
        echo -e "${RED}Base database: NOT INSTALLED${NC}"
        echo "Run setup_blast_db.sh to install the base database first."
        return 1
    fi

    # Check extensions
    echo ""
    if [[ -d "${DB_DIR}/${EXT_DIR}" ]]; then
        local ext_count=0
        for ext_nhr in "${DB_DIR}/${EXT_DIR}"/*.nhr; do
            [[ -f "$ext_nhr" ]] || continue
            ext_count=$((ext_count + 1))
            local name
            name=$(basename "$ext_nhr" .nhr)
            echo -e "${GREEN}Extension: ${name}${NC}"
            blastdbcmd -db "${DB_DIR}/${EXT_DIR}/${name}" -info 2>/dev/null | head -3 || echo "  (info unavailable)"
            # Show source FASTA if preserved
            if [[ -f "${DB_DIR}/${EXT_DIR}/${name}.fasta" ]]; then
                local seq_count
                seq_count=$(grep -c "^>" "${DB_DIR}/${EXT_DIR}/${name}.fasta" 2>/dev/null || echo "?")
                echo "  Source: ${name}.fasta (${seq_count} sequences)"
            fi
            echo ""
        done
        if [[ $ext_count -eq 0 ]]; then
            echo "User extensions: none"
        fi
    else
        echo "User extensions: none"
    fi

    # Alias status
    echo ""
    if [[ -f "${DB_DIR}/${ALIAS_NAME}.nal" ]]; then
        echo -e "${GREEN}Active database: ${ALIAS_NAME} (alias)${NC}"
    else
        echo -e "${YELLOW}Active database: ${BASE_DB_NAME} (no alias)${NC}"
    fi
}

# Rebuild the vicast_combined alias from base + all extensions
rebuild_alias() {
    local db_list="${BASE_DB_NAME}"

    if [[ -d "${DB_DIR}/${EXT_DIR}" ]]; then
        for ext_nhr in "${DB_DIR}/${EXT_DIR}"/*.nhr; do
            [[ -f "$ext_nhr" ]] || continue
            local name
            name=$(basename "$ext_nhr" .nhr)
            db_list="${db_list} ${EXT_DIR}/${name}"
        done
    fi

    echo "Rebuilding alias: ${ALIAS_NAME} -> [${db_list}]"
    (cd "$DB_DIR" && blastdb_aliastool -dblist "${db_list}" \
        -dbtype nucl \
        -out "${ALIAS_NAME}" \
        -title "VICAST Combined BLAST DB") || {
        echo -e "${RED}Error: Failed to create alias database${NC}" >&2
        return 1
    }
    echo -e "${GREEN}Alias updated successfully${NC}"
}

# Remove all user extensions
reset_db() {
    echo "Removing all user extensions..."
    if [[ -d "${DB_DIR}/${EXT_DIR}" ]]; then
        rm -rf "${DB_DIR}/${EXT_DIR}"
        echo "Removed ${DB_DIR}/${EXT_DIR}/"
    else
        echo "No extensions to remove."
    fi

    # Rebuild alias with base only
    rebuild_alias
    echo -e "${GREEN}Database reset to base only${NC}"
}

# Add FASTA files as an extension
add_extension() {
    local name="$1"
    shift
    local fasta_files=("$@")

    # Validate base DB exists
    if [[ ! -f "${DB_DIR}/${BASE_DB_NAME}.nhr" ]] && [[ ! -f "${DB_DIR}/${BASE_DB_NAME}.00.nhr" ]]; then
        echo -e "${RED}Error: Base database not installed at ${DB_DIR}/${BASE_DB_NAME}${NC}" >&2
        echo "Run setup_blast_db.sh first." >&2
        return 1
    fi

    # Validate FASTA files
    for f in "${fasta_files[@]}"; do
        if [[ ! -f "$f" ]]; then
            echo -e "${RED}Error: File not found: $f${NC}" >&2
            return 1
        fi
        # Quick FASTA validation
        if ! head -1 "$f" | grep -q "^>"; then
            echo -e "${RED}Error: $f does not appear to be FASTA format${NC}" >&2
            return 1
        fi
    done

    # Create extensions directory
    mkdir -p "${DB_DIR}/${EXT_DIR}"

    # Combine input FASTAs
    local combined="${DB_DIR}/${EXT_DIR}/${name}.fasta"
    if [[ -f "$combined" ]]; then
        echo "Appending to existing extension: ${name}"
        for f in "${fasta_files[@]}"; do
            cat "$f" >> "$combined"
        done
    else
        echo "Creating new extension: ${name}"
        cat "${fasta_files[@]}" > "$combined"
    fi

    local seq_count
    seq_count=$(grep -c "^>" "$combined")
    echo "Extension ${name}: ${seq_count} sequences"

    # Build BLAST database from the extension FASTA
    echo "Building BLAST database for extension..."
    makeblastdb \
        -in "$combined" \
        -dbtype nucl \
        -out "${DB_DIR}/${EXT_DIR}/${name}" \
        -title "User extension: ${name}" \
        -parse_seqids 2>&1 || {
        # Retry without -parse_seqids (some FASTAs have non-standard headers)
        echo "Retrying without -parse_seqids..."
        makeblastdb \
            -in "$combined" \
            -dbtype nucl \
            -out "${DB_DIR}/${EXT_DIR}/${name}" \
            -title "User extension: ${name}" 2>&1
    }

    # Rebuild the combined alias
    rebuild_alias

    echo ""
    echo -e "${GREEN}Extension '${name}' added successfully (${seq_count} sequences)${NC}"
    echo "The database now includes base + all extensions."
    echo "Run 'extend_blast_db.sh --list' to see full composition."
}

# Main
main() {
    local action="add"
    local fasta_files=()

    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h)
                usage
                ;;
            --list)
                action="list"
                shift
                ;;
            --reset)
                action="reset"
                shift
                ;;
            --name)
                EXT_NAME="$2"
                shift 2
                ;;
            --db-dir)
                DB_DIR="$2"
                shift 2
                ;;
            -*)
                echo -e "${RED}Unknown option: $1${NC}" >&2
                usage
                ;;
            *)
                fasta_files+=("$1")
                shift
                ;;
        esac
    done

    case "$action" in
        list)
            list_db
            ;;
        reset)
            reset_db
            ;;
        add)
            if [[ ${#fasta_files[@]} -eq 0 ]]; then
                echo -e "${RED}Error: No FASTA files provided${NC}" >&2
                echo "Usage: extend_blast_db.sh [--name NAME] FASTA [FASTA ...]" >&2
                exit 1
            fi
            add_extension "$EXT_NAME" "${fasta_files[@]}"
            ;;
    esac
}

main "$@"
