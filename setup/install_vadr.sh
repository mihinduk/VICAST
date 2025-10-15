#!/bin/bash
# VADR Installation Script for Viral Annotation
# Version: 1.0.0 - Initial VADR installation
#
# VADR (Viral Annotation DefineR) is used to validate and improve viral genome annotations
# by transferring annotations from well-annotated reference genomes

set -e

SCRIPT_VERSION="1.0.0"

echo "========================================="
echo "VADR Installation Script"
echo "Version: $SCRIPT_VERSION"
echo "========================================="
echo ""
echo "VADR (Viral Annotation DefineR) validates and annotates viral genomes"
echo "using reference models and covariance models."
echo ""

# Check if installation path provided
if [ -z "$1" ]; then
    echo "Usage: $0 <software_directory_path> [--check-only]"
    echo ""
    echo "Examples:"
    echo "  $0 /ref/sahlab/software"
    echo "  $0 /opt/vadr"
    echo ""
    echo "This script will install:"
    echo "  - VADR: <software_dir>/vadr"
    echo "  - Required dependencies: Perl modules, Infernal, HMMER"
    echo "  - Viral model library"
    echo ""
    echo "Options:"
    echo "  --check-only  Only check dependencies, don't install"
    echo ""
    echo "Requirements:"
    echo "  - Perl 5.12 or later"
    echo "  - curl or wget"
    echo "  - At least 2GB free space"
    exit 1
fi

SOFTWARE_DIR="$1"
CHECK_ONLY="${2:-}"

echo "Software directory: $SOFTWARE_DIR"
[ "$CHECK_ONLY" = "--check-only" ] && echo "Mode: Check only (no installation)"

# Track installation status
INSTALL_SUMMARY=""

# Function to add to summary
add_to_summary() {
    INSTALL_SUMMARY="${INSTALL_SUMMARY}$1\n"
}

# Function to check system dependencies
check_system_deps() {
    echo ""
    echo "=== Checking System Dependencies ==="
    
    # Check curl/wget
    if command -v curl >/dev/null 2>&1; then
        echo "✓ curl is available"
        DOWNLOAD_CMD="curl -L -o"
        add_to_summary "✓ curl: available"
    elif command -v wget >/dev/null 2>&1; then
        echo "✓ wget is available"  
        DOWNLOAD_CMD="wget -O"
        add_to_summary "✓ wget: available"
    else
        echo "✗ Neither curl nor wget found"
        add_to_summary "✗ curl/wget: MISSING"
        if [ "$CHECK_ONLY" != "--check-only" ]; then
            echo "Please install curl or wget:"
            echo "  sudo yum install curl  # RHEL/CentOS"
            echo "  sudo apt-get install curl  # Ubuntu/Debian"
            exit 1
        fi
    fi
    
    # Check Perl
    if command -v perl >/dev/null 2>&1; then
        PERL_VERSION=$(perl -version | grep -oP 'v\d+\.\d+\.\d+' | head -1)
        echo "✓ Perl available: $PERL_VERSION"
        add_to_summary "✓ Perl: $PERL_VERSION"
    else
        echo "✗ Perl not found"
        add_to_summary "✗ Perl: MISSING"
        if [ "$CHECK_ONLY" != "--check-only" ]; then
            echo "Please install Perl:"
            echo "  sudo yum install perl  # RHEL/CentOS"
            echo "  sudo apt-get install perl  # Ubuntu/Debian"
            exit 1
        fi
    fi
    
    # Check unzip/tar
    if command -v tar >/dev/null 2>&1; then
        echo "✓ tar is available"
        add_to_summary "✓ tar: available"
    else
        echo "✗ tar not found"
        add_to_summary "✗ tar: MISSING"
    fi
    
    # Check make
    if command -v make >/dev/null 2>&1; then
        echo "✓ make is available"
        add_to_summary "✓ make: available"
    else
        echo "⚠ make not found (may be needed for dependencies)"
        add_to_summary "⚠ make: MISSING"
    fi
}

# Function to install VADR
install_vadr() {
    echo ""
    echo "=== Installing VADR ==="
    
    VADR_VERSION="1.6.3"
    VADR_DIR="$SOFTWARE_DIR/vadr"
    
    if [ -d "$VADR_DIR" ] && [ -f "$VADR_DIR/vadr.sh" ]; then
        echo "✓ VADR already installed at $VADR_DIR"
        add_to_summary "✓ VADR: already installed"
        return 0
    fi
    
    echo "✗ VADR not found at $VADR_DIR"
    add_to_summary "✗ VADR: NOT INSTALLED"
    
    if [ "$CHECK_ONLY" = "--check-only" ]; then
        return 1
    fi
    
    echo "Installing VADR v$VADR_VERSION..."
    mkdir -p "$SOFTWARE_DIR"
    
    # Download VADR
    VADR_URL="https://github.com/ncbi/vadr/archive/vadr-$VADR_VERSION.tar.gz"
    echo "Downloading VADR from $VADR_URL"
    
    if [ "$DOWNLOAD_CMD" = "curl -L -o" ]; then
        $DOWNLOAD_CMD "$SOFTWARE_DIR/vadr-$VADR_VERSION.tar.gz" "$VADR_URL"
    else
        $DOWNLOAD_CMD "$SOFTWARE_DIR/vadr-$VADR_VERSION.tar.gz" "$VADR_URL"
    fi
    
    # Extract VADR
    echo "Extracting VADR..."
    cd "$SOFTWARE_DIR"
    tar -xzf "vadr-$VADR_VERSION.tar.gz"
    mv "vadr-vadr-$VADR_VERSION" "$VADR_DIR"
    rm "vadr-$VADR_VERSION.tar.gz"
    
    echo "✓ VADR extracted to $VADR_DIR"
    
    # Run VADR installation script
    echo "Running VADR installation script..."
    cd "$VADR_DIR"
    
    # Set installation environment
    export VADRINSTALLDIR="$VADR_DIR"
    
    # Run the install script
    echo "This may take several minutes to compile dependencies..."
    if ./vadr-install.sh linux > vadr-install.log 2>&1; then
        echo "✓ VADR installation completed successfully"
        add_to_summary "✓ VADR: INSTALLED (v$VADR_VERSION)"
    else
        echo "✗ VADR installation failed. Check log: $VADR_DIR/vadr-install.log"
        echo "Common issues:"
        echo "  - Missing development tools (gcc, make)"
        echo "  - Insufficient disk space"
        echo "  - Network connectivity issues"
        add_to_summary "✗ VADR: INSTALLATION FAILED"
        return 1
    fi
    
    # Download VADR models
    echo ""
    echo "Downloading VADR models..."
    if [ "$DOWNLOAD_CMD" = "curl -L -o" ]; then
        $DOWNLOAD_CMD "$VADR_DIR/vadr-models.tar.gz" "https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/vadr-models-calici-1.4.2-1.tar.gz"
    else
        $DOWNLOAD_CMD "$VADR_DIR/vadr-models.tar.gz" "https://ftp.ncbi.nlm.nih.gov/pub/nawrocki/vadr-models/vadr-models-calici-1.4.2-1.tar.gz"
    fi
    
    # Extract models
    echo "Extracting VADR models..."
    tar -xzf vadr-models.tar.gz
    rm vadr-models.tar.gz
    
    echo "✓ VADR models installed"
}

# Function to create VADR configuration
create_vadr_config() {
    echo ""
    echo "=== Creating VADR Configuration ==="
    
    VADR_DIR="$SOFTWARE_DIR/vadr"
    CONFIG_DIR="$SOFTWARE_DIR/snpeff_configs"
    
    if [ "$CHECK_ONLY" = "--check-only" ]; then
        if [ -f "$CONFIG_DIR/vadr_env.sh" ]; then
            echo "✓ VADR configuration exists"
            add_to_summary "✓ VADR Config: exists"
        else
            echo "✗ VADR configuration not found"
            add_to_summary "✗ VADR Config: NOT FOUND"
        fi
        return 0
    fi
    
    mkdir -p "$CONFIG_DIR"
    
    # Create VADR environment configuration
    cat > "$CONFIG_DIR/vadr_env.sh" << CONFIG_EOF
#!/bin/bash
# VADR Environment Configuration
# Version: $SCRIPT_VERSION
# Generated: $(date)

export VADR_HOME="$VADR_DIR"
export PATH="\$VADR_HOME:\$PATH"

# Source VADR setup if available
if [ -f "\$VADR_HOME/vadr-env.sh" ]; then
    source "\$VADR_HOME/vadr-env.sh"
fi

# VADR utility functions
vadr_version() {
    echo "VADR Home: \$VADR_HOME"
    if command -v v-annotate.pl >/dev/null 2>&1; then
        v-annotate.pl -h 2>&1 | head -3
    else
        echo "VADR commands not found in PATH"
    fi
}

vadr_annotate() {
    if [ -z "\$1" ]; then
        echo "Usage: vadr_annotate <input.fasta> [output_dir] [model]"
        echo ""
        echo "Examples:"
        echo "  vadr_annotate OP713603.1.fasta"
        echo "  vadr_annotate OP713603.1.fasta my_output calici"
        return 1
    fi
    
    local input_fasta="\$1"
    local output_dir="\${2:-vadr_output}"
    local model="\${3:-calici}"
    
    echo "Running VADR annotation..."
    echo "  Input: \$input_fasta"
    echo "  Output: \$output_dir"
    echo "  Model: \$model"
    
    v-annotate.pl --split --cpu 4 --glsearch -s -r --nomisc \\
                  --mkey \$model --lowsim5seq 6 --lowsim3seq 6 \\
                  --alt_fail lowscore,insertnn,deletinn \\
                  --mdir "\$VADR_HOME/vadr-models-\$model" \\
                  "\$input_fasta" "\$output_dir"
}

# Display configuration when sourced
if [ -n "\$PS1" ]; then
    echo "VADR environment loaded!"
    echo "  VADR Home: \$VADR_HOME"
    echo ""
    echo "Functions: vadr_version, vadr_annotate"
    echo "Example: vadr_annotate my_genome.fasta"
fi
CONFIG_EOF
    
    chmod +x "$CONFIG_DIR/vadr_env.sh"
    
    echo "✓ VADR configuration created: $CONFIG_DIR/vadr_env.sh"
    add_to_summary "✓ VADR Config: CREATED"
}

# Function to display summary
display_summary() {
    echo ""
    echo "========================================="
    echo "Installation Summary"
    echo "========================================="
    echo -e "$INSTALL_SUMMARY"
    
    if [ "$CHECK_ONLY" != "--check-only" ]; then
        echo ""
        echo "To use VADR:"
        echo "  source $SOFTWARE_DIR/snpeff_configs/vadr_env.sh"
        echo ""
        echo "Test installation:"
        echo "  vadr_version"
        echo ""
        echo "Example annotation workflow:"
        echo "  # Download a poorly-annotated genome"
        echo "  download_ncbi_genome OP713603.1"
        echo "  "
        echo "  # Annotate using VADR with NC_009942.1 as reference"
        echo "  vadr_annotate OP713603.1.fasta OP713603_vadr calici"
        echo "  "
        echo "  # Check results in OP713603_vadr/ directory"
        echo "  ls OP713603_vadr/"
        echo ""
        echo "Integration with snpEff workflow:"
        echo "  # After VADR annotation, use the .gff file for snpEff"
        echo "  python3 step1_parse_viral_genome.py OP713603_vadr/OP713603.1.vadr.gff"
    fi
}

# Main execution
main() {
    echo "Script version: $SCRIPT_VERSION"
    
    # Check permissions
    if [ "$CHECK_ONLY" != "--check-only" ]; then
        if [ ! -w "$SOFTWARE_DIR" ] && [ ! -d "$SOFTWARE_DIR" ]; then
            echo "Creating software directory: $SOFTWARE_DIR"
            mkdir -p "$SOFTWARE_DIR" 2>/dev/null || {
                echo "Cannot create $SOFTWARE_DIR. You may need sudo:"
                echo "  sudo mkdir -p $SOFTWARE_DIR"
                echo "  sudo chown $USER:$(id -gn) $SOFTWARE_DIR"
                exit 1
            }
        fi
    fi
    
    check_system_deps
    install_vadr
    create_vadr_config
    display_summary
}

# Run main function
main