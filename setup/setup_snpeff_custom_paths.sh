#!/bin/bash
# Setup script for snpEff with custom paths and dependency installation
# Version: 2.0.0 - Enhanced with full pipeline dependencies
# 
# Git versioning strategy:
#   - Tag releases: git tag -a v2.0.0 -m "Add pipeline dependencies"
#   - Push tags: git push origin v2.0.0
#   - View version: git describe --tags

set -e

SCRIPT_VERSION="2.1.0"

echo "========================================="
echo "snpEff Setup Script with Custom Paths"
echo "Version: $SCRIPT_VERSION"
echo "========================================="
echo ""

# Check if software path provided
if [ -z "$1" ]; then
    echo "Usage: $0 <software_directory_path> [--check-only]"
    echo "   or: SCRATCH_DIR=<scratch_path> $0 <software_directory_path> [--check-only]"
    echo ""
    echo "Examples:"
    echo "  $0 /ref/sahlab/software"
    echo "  SCRATCH_DIR=/scratch/sahlab/kathie $0 /ref/sahlab/software"
    echo ""
    echo "This script will install/check:"
    echo "  - Java 21: <software_dir>/jdk-21.0.5+11"
    echo "  - snpEff: <software_dir>/snpEff" 
    echo "  - Python packages: pandas (scratch space install)"
    echo "  - System tools: curl/wget"
    echo "  - Pipeline configs: <software_dir>/snpeff_configs"
    echo ""
    echo "Options:"
    echo "  --check-only     Only check dependencies, don't install"
    echo "  SCRATCH_DIR=<path>  Specify scratch directory for Python packages"
    echo "                     (will prompt if not specified)"
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

# Function to check and install system packages
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
    
    # Check Python
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
        echo "✓ Python 3 available: $PYTHON_VERSION"
        add_to_summary "✓ Python: $PYTHON_VERSION"
    else
        echo "✗ Python 3 not found"
        add_to_summary "✗ Python 3: MISSING"
        if [ "$CHECK_ONLY" != "--check-only" ]; then
            echo "Please install Python 3:"
            echo "  sudo yum install python3  # RHEL/CentOS"
            echo "  sudo apt-get install python3  # Ubuntu/Debian"
            exit 1
        fi
    fi
    
    # Check unzip
    if command -v unzip >/dev/null 2>&1; then
        echo "✓ unzip is available"
        add_to_summary "✓ unzip: available"
    else
        echo "⚠ unzip not found (needed for snpEff installation)"
        add_to_summary "⚠ unzip: MISSING"
    fi
}

# Function to check/install Java 21
check_install_java21() {
    echo ""
    echo "=== Java 21 ==="
    
    JAVA_DIR="$SOFTWARE_DIR/jdk-21.0.5+11"
    
    if [ -f "$JAVA_DIR/bin/java" ]; then
        JAVA_VERSION=$("$JAVA_DIR/bin/java" -version 2>&1 | head -1 | awk -F'"' '{print $2}')
        echo "✓ Java already installed: $JAVA_VERSION at $JAVA_DIR"
        add_to_summary "✓ Java 21: $JAVA_VERSION"
        return 0
    fi
    
    echo "✗ Java 21 not found at $JAVA_DIR"
    add_to_summary "✗ Java 21: NOT INSTALLED"
    
    if [ "$CHECK_ONLY" = "--check-only" ]; then
        return 1
    fi
    
    echo "Installing Java 21..."
    mkdir -p "$SOFTWARE_DIR"
    
    if [ "$DOWNLOAD_CMD" = "curl -L -o" ]; then
        $DOWNLOAD_CMD "$SOFTWARE_DIR/openjdk21.tar.gz" "https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.5%2B11/OpenJDK21U-jdk_x64_linux_hotspot_21.0.5_11.tar.gz"
    else
        $DOWNLOAD_CMD "$SOFTWARE_DIR/openjdk21.tar.gz" "https://github.com/adoptium/temurin21-binaries/releases/download/jdk-21.0.5%2B11/OpenJDK21U-jdk_x64_linux_hotspot_21.0.5_11.tar.gz"
    fi
    
    echo "Extracting Java 21..."
    tar -xzf "$SOFTWARE_DIR/openjdk21.tar.gz" -C "$SOFTWARE_DIR"
    rm "$SOFTWARE_DIR/openjdk21.tar.gz"
    
    echo "✓ Java 21 installed to $JAVA_DIR"
    add_to_summary "✓ Java 21: INSTALLED"
}

# Function to check/install snpEff
check_install_snpeff() {
    echo ""
    echo "=== snpEff ==="
    
    SNPEFF_DIR="$SOFTWARE_DIR/snpEff"
    
    if [ -f "$SNPEFF_DIR/snpEff.jar" ]; then
        # Try to get version
        if [ -f "$JAVA_DIR/bin/java" ]; then
            SNPEFF_VERSION=$("$JAVA_DIR/bin/java" -jar "$SNPEFF_DIR/snpEff.jar" -version 2>&1 | head -1 || echo "unknown")
            echo "✓ snpEff already installed at $SNPEFF_DIR"
            echo "  Version: $SNPEFF_VERSION"
            add_to_summary "✓ snpEff: $SNPEFF_VERSION"
        else
            echo "✓ snpEff found at $SNPEFF_DIR (version check requires Java)"
            add_to_summary "✓ snpEff: found"
        fi
        return 0
    fi
    
    echo "✗ snpEff not found at $SNPEFF_DIR"
    add_to_summary "✗ snpEff: NOT INSTALLED"
    
    if [ "$CHECK_ONLY" = "--check-only" ]; then
        return 1
    fi
    
    echo "Installing snpEff..."
    
    if [ "$DOWNLOAD_CMD" = "curl -L -o" ]; then
        $DOWNLOAD_CMD "$SOFTWARE_DIR/snpEff_latest_core.zip" "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"
    else
        $DOWNLOAD_CMD "$SOFTWARE_DIR/snpEff_latest_core.zip" "https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip"
    fi
    
    echo "Extracting snpEff..."
    unzip -q "$SOFTWARE_DIR/snpEff_latest_core.zip" -d "$SOFTWARE_DIR"
    rm "$SOFTWARE_DIR/snpEff_latest_core.zip"
    
    echo "✓ snpEff installed to $SNPEFF_DIR"
    add_to_summary "✓ snpEff: INSTALLED"
}

# Function to check/install Python packages
check_install_python_deps() {
    echo ""
    echo "=== Python Dependencies ==="
    
    # Setup PIP cache and install location in scratch space to avoid home directory quota issues
    if [ -n "$SCRATCH_DIR" ]; then
        # Use user-provided scratch directory
        export PIP_CACHE_DIR="$SCRATCH_DIR/.cache/pip"
        export PIP_USER_DIR="$SCRATCH_DIR/.local"
    elif [ -n "$USER" ]; then
        # Fallback: ask user for scratch directory
        echo ""
        echo "⚠ Scratch directory not specified. Please provide your scratch space path."
        echo "Common examples:"
        echo "  - /scratch/\$USER"
        echo "  - /scratch/sahlab/\$USER"
        echo "  - /tmp/\$USER"
        read -p "Enter your scratch directory path: " SCRATCH_INPUT
        
        if [ -n "$SCRATCH_INPUT" ]; then
            export PIP_CACHE_DIR="$SCRATCH_INPUT/.cache/pip"
            export PIP_USER_DIR="$SCRATCH_INPUT/.local"
        else
            echo "⚠ No scratch directory provided, using home directory (may hit quota limits)"
            export PIP_USER_DIR="$HOME/.local"
            export PIP_CACHE_DIR="$HOME/.cache/pip"
        fi
    else
        echo "⚠ Warning: \$USER not set, using default PIP locations"
        export PIP_USER_DIR="$HOME/.local"
        export PIP_CACHE_DIR="$HOME/.cache/pip"
    fi
    
    export PYTHONUSERBASE="$PIP_USER_DIR"
    export PATH="$PIP_USER_DIR/bin:$PATH"
    
    echo "Setting PIP cache directory: $PIP_CACHE_DIR"
    echo "Setting PIP install directory: $PIP_USER_DIR"
    
    if [ "$CHECK_ONLY" != "--check-only" ]; then
        mkdir -p "$PIP_CACHE_DIR" "$PIP_USER_DIR" 2>/dev/null || {
            echo "⚠ Warning: Could not create PIP directories in scratch space"
            echo "  This may cause slower installs but won't prevent installation"
        }
    fi
    
    # Check pandas (check both standard location and scratch space)
    PANDAS_FOUND=false
    if python3 -c "import pandas" >/dev/null 2>&1; then
        PANDAS_VERSION=$(python3 -c "import pandas; print(pandas.__version__)")
        PANDAS_LOCATION=$(python3 -c "import pandas; print(pandas.__file__)")
        echo "✓ pandas available: $PANDAS_VERSION"
        echo "  Location: $PANDAS_LOCATION"
        add_to_summary "✓ pandas: $PANDAS_VERSION"
        PANDAS_FOUND=true
    fi
    
    # Check biopython
    if python3 -c "import Bio" >/dev/null 2>&1; then
        BIO_VERSION=$(python3 -c "import Bio; print(Bio.__version__)")
        echo "✓ biopython available: $BIO_VERSION"
        add_to_summary "✓ biopython: $BIO_VERSION"
    else
        echo "✗ biopython not available"
        add_to_summary "✗ biopython: NOT INSTALLED"
    fi
    
    # Install missing packages
    PACKAGES_TO_INSTALL=""
    if [ "$PANDAS_FOUND" = false ]; then
        echo "✗ pandas not available"
        add_to_summary "✗ pandas: NOT INSTALLED"
        PACKAGES_TO_INSTALL="pandas"
    fi
    
    # Check if biopython needs installation
    BIO_FOUND=false
    if python3 -c "import Bio" >/dev/null 2>&1; then
        BIO_FOUND=true
    else
        echo "✗ biopython needs installation"
        if [ -n "$PACKAGES_TO_INSTALL" ]; then
            PACKAGES_TO_INSTALL="$PACKAGES_TO_INSTALL biopython"
        else
            PACKAGES_TO_INSTALL="biopython"
        fi
    fi
    
    # Install packages if any are missing
    if [ -n "$PACKAGES_TO_INSTALL" ] && [ "$CHECK_ONLY" != "--check-only" ]; then
        echo "Installing $PACKAGES_TO_INSTALL to scratch space..."
        echo "Target directory: $PIP_USER_DIR"
        
        # Try different installation methods
        if command -v pip3 >/dev/null 2>&1; then
            echo "Using pip3 with scratch space installation..."
            env PYTHONUSERBASE="$PIP_USER_DIR" PIP_CACHE_DIR="$PIP_CACHE_DIR" pip3 install --user $PACKAGES_TO_INSTALL
        elif command -v pip >/dev/null 2>&1; then
            echo "Using pip with scratch space installation..."
            env PYTHONUSERBASE="$PIP_USER_DIR" PIP_CACHE_DIR="$PIP_CACHE_DIR" pip install --user $PACKAGES_TO_INSTALL  
        elif command -v conda >/dev/null 2>&1; then
            echo "Using conda..."
            conda install $PACKAGES_TO_INSTALL -y
        else
            echo "⚠ Warning: No package manager found (pip, conda)"
            echo "Please install packages manually:"
            echo "  pip3 install --user $PACKAGES_TO_INSTALL"
            return 1
        fi
        
        # Verify pandas installation if it was installed
        if [[ "$PACKAGES_TO_INSTALL" == *"pandas"* ]]; then
            if python3 -c "import pandas" >/dev/null 2>&1; then
                PANDAS_VERSION=$(python3 -c "import pandas; print(pandas.__version__)")
                PANDAS_LOCATION=$(python3 -c "import pandas; print(pandas.__file__)")
                echo "✓ pandas successfully installed: $PANDAS_VERSION"
                echo "  Installed to: $PANDAS_LOCATION"
                # Update summary to replace the NOT INSTALLED entry
                INSTALL_SUMMARY=$(echo -e "$INSTALL_SUMMARY" | sed 's/✗ pandas: NOT INSTALLED/✓ pandas: INSTALLED ('"$PANDAS_VERSION"')/')
            else
                echo "✗ pandas installation failed"
            fi
        fi
        
        # Verify Biopython installation if it was installed
        if [[ "$PACKAGES_TO_INSTALL" == *"biopython"* ]]; then
            if python3 -c "import Bio" >/dev/null 2>&1; then
                BIO_VERSION=$(python3 -c "import Bio; print(Bio.__version__)")
                echo "✓ biopython successfully installed: $BIO_VERSION"
                # Update summary to replace the NOT INSTALLED entry
                INSTALL_SUMMARY=$(echo -e "$INSTALL_SUMMARY" | sed 's/✗ biopython: NOT INSTALLED/✓ biopython: INSTALLED ('"$BIO_VERSION"')/')
            else
                echo "✗ biopython installation failed"
                add_to_summary "✗ biopython: FAILED"
            fi
        fi
    fi
    
    # Check other useful packages
    echo ""
    echo "Optional Python packages:"
    
    for pkg in requests numpy; do
        if python3 -c "import $pkg" >/dev/null 2>&1; then
            VERSION=$(python3 -c "import $pkg; print($pkg.__version__)" 2>/dev/null || echo "unknown")
            echo "  ✓ $pkg: $VERSION"
        else
            echo "  ○ $pkg: not installed (optional)"
        fi
    done
}

# Function to create/update pipeline configuration
create_update_config() {
    echo ""
    echo "=== Pipeline Configuration ==="
    
    CONFIG_DIR="$SOFTWARE_DIR/snpeff_configs"
    JAVA_DIR="$SOFTWARE_DIR/jdk-21.0.5+11"
    SNPEFF_DIR="$SOFTWARE_DIR/snpEff"
    
    if [ "$CHECK_ONLY" = "--check-only" ]; then
        if [ -f "$CONFIG_DIR/snpeff_env.sh" ]; then
            echo "✓ Configuration exists: $CONFIG_DIR/snpeff_env.sh"
            add_to_summary "✓ Config: exists"
        else
            echo "✗ Configuration not found"
            add_to_summary "✗ Config: NOT FOUND"
        fi
        return 0
    fi
    
    mkdir -p "$CONFIG_DIR"
    
    # Create versioned configuration with metadata
    cat > "$CONFIG_DIR/snpeff_env.sh" << CONFIG_EOF
#!/bin/bash
# snpEff Pipeline Configuration
# Version: $SCRIPT_VERSION
# Generated: $(date)
# Git commit: $(git rev-parse --short HEAD 2>/dev/null || echo "unknown")

export SNPEFF_PIPELINE_VERSION="$SCRIPT_VERSION"
export JAVA_HOME="$JAVA_DIR"
export SNPEFF_HOME="$SNPEFF_DIR"

# Python package management in scratch space
# Note: Set SCRATCH_DIR environment variable or script will prompt for scratch path
export PIP_CACHE_DIR="\${SCRATCH_DIR:-/scratch/\$USER}/.cache/pip"
export PIP_USER_DIR="\${SCRATCH_DIR:-/scratch/\$USER}/.local"
export PYTHONUSERBASE="\$PIP_USER_DIR"

# Update PATH to include both Java and Python packages
export PATH="\$JAVA_HOME/bin:\$PIP_USER_DIR/bin:\$PATH"

# Core functions
snpeff() {
    "\$JAVA_HOME/bin/java" -Xmx4g -jar "\$SNPEFF_HOME/snpEff.jar" "\$@"
}

snpsift() {
    "\$JAVA_HOME/bin/java" -Xmx4g -jar "\$SNPEFF_HOME/SnpSift.jar" "\$@"
}

# NCBI download function
download_ncbi_genome() {
    local accession="\$1"
    if [ -z "\$accession" ]; then
        echo "Usage: download_ncbi_genome <NCBI_ACCESSION>"
        return 1
    fi
    
    echo "Downloading \$accession from NCBI..."
    
    # Download FASTA
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\$accession&rettype=fasta&retmode=text" > "\$accession.fasta"
    
    # Download GenBank
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\$accession&rettype=gb&retmode=text" > "\$accession.gb"
    
    # Download GFF3
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=\$accession&rettype=gff3&retmode=text" > "\$accession.gff3"
    
    echo "Downloaded: \$accession.fasta, \$accession.gb, \$accession.gff3"
}

# Display version
snpeff_pipeline_version() {
    echo "snpEff Pipeline Version: \$SNPEFF_PIPELINE_VERSION"
    echo "Java: \$JAVA_HOME"
    echo "snpEff: \$SNPEFF_HOME"
    echo "Python packages: \$PYTHONUSERBASE"
}

# Display configuration when sourced
if [ -n "\$PS1" ]; then
    echo "snpEff Pipeline loaded! (v\$SNPEFF_PIPELINE_VERSION)"
    echo "  Java: \$JAVA_HOME"
    echo "  snpEff: \$SNPEFF_HOME"
    echo "  Python: \$PYTHONUSERBASE"
    echo ""
    echo "Functions: snpeff, snpsift, download_ncbi_genome, snpeff_pipeline_version"
fi
CONFIG_EOF
    
    chmod +x "$CONFIG_DIR/snpeff_env.sh"
    
    # Create/update current symlink
    ln -sf "$CONFIG_DIR/snpeff_env.sh" "$CONFIG_DIR/snpeff_current.sh"
    
    echo "✓ Configuration created/updated"
    add_to_summary "✓ Config: CREATED"
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
        echo "To use the pipeline:"
        echo "  source $SOFTWARE_DIR/snpeff_configs/snpeff_current.sh"
        echo ""
        echo "Note: Python packages are installed in scratch space (/scratch/\$USER/.local/)"
        echo "This avoids home directory quota issues on HPC systems."
        echo ""
        echo "Version this installation in Git:"
        echo "  git add setup_snpeff_custom_paths.sh"
        echo "  git commit -m 'Update snpEff setup v$SCRIPT_VERSION'"
        echo "  git tag -a v$SCRIPT_VERSION -m 'Scratch space support for Python packages'"
        echo "  git push origin main --tags"
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
    check_install_java21
    check_install_snpeff  
    check_install_python_deps
    create_update_config
    display_summary
}

# Function to test VICAST annotation pipeline
test_vicast_annotate() {
    echo ""
    echo "========================================="
    echo "Testing VICAST Annotation Pipeline"
    echo "========================================="
    
    # Check if required scripts exist
    local MISSING_SCRIPTS=""
    for script in step1_parse_viral_genome.py step2_add_to_snpeff.py vicast_validation.py; do
        if [ ! -f "$script" ]; then
            MISSING_SCRIPTS="$MISSING_SCRIPTS $script"
        fi
    done
    
    if [ -n "$MISSING_SCRIPTS" ]; then
        echo "[ERROR] Missing required scripts:$MISSING_SCRIPTS"
        echo "Please ensure all VICAST scripts are in the current directory"
        return 1
    fi
    
    # Test with a small, well-annotated virus
    local TEST_GENOME="${1:-NC_001477}"  # Dengue virus type 1 by default
    
    echo "Using test genome: $TEST_GENOME"
    echo ""
    
    # Create test directory
    local TEST_DIR="vicast_test_$(date +%Y%m%d_%H%M%S)"
    mkdir -p "$TEST_DIR"
    cd "$TEST_DIR" || return 1
    
    echo "Step 1: Downloading test genome ${TEST_GENOME}..."
    download_ncbi_genome "${TEST_GENOME}"
    
    if [ ! -f "${TEST_GENOME}.fasta" ] || [ ! -f "${TEST_GENOME}.gb" ]; then
        echo "[ERROR] Failed to download genome files"
        cd ..
        return 1
    fi
    
    echo "Step 2: Parsing genome to TSV..."
    python3 ../step1_parse_viral_genome.py "${TEST_GENOME}"
    
    if [ ! -f "${TEST_GENOME}_no_polyprotein.tsv" ]; then
        echo "[ERROR] TSV file not created"
        cd ..
        return 1
    fi
    
    echo "Step 3: Running validation..."
    if [ -f "../vicast_validation.py" ]; then
        python3 -c "import sys; sys.path.insert(0, '..'); from vicast_validation import validate_gff_for_snpeff; print('Validation module loaded successfully')"
        
        # Convert TSV to test GFF for validation
        python3 -c "
import sys
import pandas as pd
sys.path.insert(0, '..')

# Read TSV
df = pd.read_csv('${TEST_GENOME}_no_polyprotein.tsv', sep='\t')

# Write minimal GFF for testing
with open('test.gff3', 'w') as f:
    f.write('##gff-version 3\n')
    for _, row in df.iterrows():
        attrs = []
        if pd.notna(row.get('ID', '')):
            attrs.append(f\"ID={row['ID']}\")
        if pd.notna(row.get('gene_name', '')):
            attrs.append(f\"gene={row['gene_name']}\")
        line = f\"{row['seqid']}\t{row['source']}\t{row['type']}\t{row['start']}\t{row['end']}\t.\t{row['strand']}\t.\t{';'.join(attrs)}\n\"
        f.write(line)
"
        
        if [ -f "test.gff3" ]; then
            python3 ../vicast_validation.py test.gff3 "${TEST_GENOME}.fasta"
        fi
    fi
    
    echo ""
    echo "Test files created in $TEST_DIR:"
    ls -lh "${TEST_GENOME}"*
    
    echo ""
    echo "[SUCCESS] Pipeline test completed!"
    echo ""
    echo "Next steps:"
    echo "1. Review and edit: ${TEST_GENOME}_no_polyprotein.tsv"
    echo "2. Add to snpEff: python3 ../step2_add_to_snpeff.py ${TEST_GENOME} ${TEST_GENOME}_no_polyprotein.tsv --report"
    echo "3. Clean up test: rm -rf $TEST_DIR"
    
    cd ..
}

# Function to validate VICAST setup
validate_vicast_setup() {
    echo ""
    echo "========================================="
    echo "Validating VICAST Setup"
    echo "========================================="
    
    local ERRORS=0
    local WARNINGS=0
    
    # Check Python version
    echo "Checking Python..."
    if command -v python3 >/dev/null 2>&1; then
        PYTHON_VERSION=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}')")
        echo "  ✓ Python: $PYTHON_VERSION"
        
        # Check if version is adequate (3.6+)
        PYTHON_MAJOR=$(python3 -c "import sys; print(sys.version_info.major)")
        PYTHON_MINOR=$(python3 -c "import sys; print(sys.version_info.minor)")
        if [ "$PYTHON_MAJOR" -eq 3 ] && [ "$PYTHON_MINOR" -lt 6 ]; then
            echo "    ⚠ Warning: Python 3.6+ recommended (you have $PYTHON_VERSION)"
            ((WARNINGS++))
        fi
    else
        echo "  ✗ Python 3 not found"
        ((ERRORS++))
    fi
    
    # Check Python modules
    echo ""
    echo "Checking Python modules..."
    
    for module in pandas Bio; do
        if python3 -c "import $module" 2>/dev/null; then
            VERSION=$(python3 -c "import $module; print($module.__version__)" 2>/dev/null || echo "installed")
            echo "  ✓ $module: $VERSION"
        else
            echo "  ✗ $module: NOT INSTALLED"
            ((ERRORS++))
        fi
    done
    
    # Check optional modules
    for module in requests numpy; do
        if python3 -c "import $module" 2>/dev/null; then
            VERSION=$(python3 -c "import $module; print($module.__version__)" 2>/dev/null || echo "installed")
            echo "  ✓ $module: $VERSION (optional)"
        else
            echo "  ○ $module: not installed (optional)"
        fi
    done
    
    # Check snpEff environment
    echo ""
    echo "Checking snpEff environment..."
    
    if [ -n "$SNPEFF_HOME" ]; then
        echo "  ✓ SNPEFF_HOME: $SNPEFF_HOME"
        if [ -f "$SNPEFF_HOME/snpEff.jar" ]; then
            echo "  ✓ snpEff.jar found"
        else
            echo "  ✗ snpEff.jar not found at $SNPEFF_HOME"
            ((ERRORS++))
        fi
    else
        echo "  ✗ SNPEFF_HOME not set"
        echo "    Run: source $SOFTWARE_DIR/snpeff_configs/snpeff_current.sh"
        ((ERRORS++))
    fi
    
    if [ -n "$JAVA_HOME" ]; then
        echo "  ✓ JAVA_HOME: $JAVA_HOME"
        if [ -f "$JAVA_HOME/bin/java" ]; then
            JAVA_VERSION=$("$JAVA_HOME/bin/java" -version 2>&1 | head -1)
            echo "  ✓ Java: $JAVA_VERSION"
        else
            echo "  ✗ Java not found at $JAVA_HOME"
            ((ERRORS++))
        fi
    else
        echo "  ✗ JAVA_HOME not set"
        ((ERRORS++))
    fi
    
    # Check VICAST scripts
    echo ""
    echo "Checking VICAST scripts..."
    
    for script in step1_parse_viral_genome.py step2_add_to_snpeff.py vicast_validation.py; do
        if [ -f "$script" ]; then
            echo "  ✓ $script"
        else
            echo "  ⚠ $script not found in current directory"
            ((WARNINGS++))
        fi
    done
    
    # Check for download capability
    echo ""
    echo "Checking download tools..."
    if command -v curl >/dev/null 2>&1; then
        echo "  ✓ curl available"
    elif command -v wget >/dev/null 2>&1; then
        echo "  ✓ wget available"
    else
        echo "  ✗ No download tool (curl/wget) found"
        ((ERRORS++))
    fi
    
    # Check disk space in scratch directory
    echo ""
    echo "Checking disk space..."
    if [ -n "$SCRATCH_DIR" ]; then
        if [ -d "$SCRATCH_DIR" ]; then
            SPACE_AVAIL=$(df -h "$SCRATCH_DIR" | tail -1 | awk '{print $4}')
            echo "  ✓ Scratch directory: $SCRATCH_DIR"
            echo "    Available space: $SPACE_AVAIL"
        else
            echo "  ⚠ Scratch directory does not exist: $SCRATCH_DIR"
            ((WARNINGS++))
        fi
    else
        echo "  ⚠ SCRATCH_DIR not set (using home directory for Python packages)"
        ((WARNINGS++))
    fi
    
    # Summary
    echo ""
    echo "========================================="
    echo "Validation Summary"
    echo "========================================="
    
    if [ $ERRORS -eq 0 ] && [ $WARNINGS -eq 0 ]; then
        echo "[SUCCESS] All checks passed!"
        echo "Your VICAST setup is ready to use."
    elif [ $ERRORS -eq 0 ]; then
        echo "[SUCCESS] Setup is functional with $WARNINGS warning(s)"
        echo "Your VICAST setup should work but consider addressing warnings."
    else
        echo "[FAILED] Found $ERRORS error(s) and $WARNINGS warning(s)"
        echo "Please fix errors before using VICAST pipeline."
        return 1
    fi
    
    return 0
}

# Function to initialize a new VICAST project
init_vicast_project() {
    local PROJECT_NAME="${1:-vicast_project}"
    
    echo ""
    echo "========================================="
    echo "Initializing VICAST Project: $PROJECT_NAME"
    echo "========================================="
    
    # Create project directory
    if [ -d "$PROJECT_NAME" ]; then
        echo "[ERROR] Directory already exists: $PROJECT_NAME"
        return 1
    fi
    
    mkdir -p "$PROJECT_NAME"/{genomes,annotations,reports,logs}
    
    # Copy VICAST scripts if they exist
    for script in step1_parse_viral_genome.py step2_add_to_snpeff.py vicast_validation.py; do
        if [ -f "$script" ]; then
            cp "$script" "$PROJECT_NAME/"
            echo "  ✓ Copied $script"
        fi
    done
    
    # Create project README
    cat > "$PROJECT_NAME/README.md" << 'README_EOF'
# VICAST Annotation Project

## Directory Structure
- `genomes/` - Downloaded genome files (FASTA, GenBank, GFF)
- `annotations/` - Edited TSV files and final GFF files
- `reports/` - Curation reports and validation logs
- `logs/` - Processing logs

## Workflow

### 1. Setup Environment
```bash
source /path/to/snpeff_configs/snpeff_current.sh
```

### 2. Download Genome
```bash
cd genomes
download_ncbi_genome NC_XXXXXX
cd ..
```

### 3. Parse Genome
```bash
python3 step1_parse_viral_genome.py genomes/NC_XXXXXX
mv *.tsv annotations/
```

### 4. Edit TSV File
Edit the TSV file in `annotations/` to curate annotations

### 5. Add to snpEff
```bash
python3 step2_add_to_snpeff.py NC_XXXXXX annotations/NC_XXXXXX_no_polyprotein.tsv \
    --fasta genomes/NC_XXXXXX.fasta \
    --gb genomes/NC_XXXXXX.gb \
    --report
mv *_report.txt reports/
```

## Notes
- Always validate your edits before adding to snpEff
- Keep original files in `genomes/` directory
- Document changes in commit messages
README_EOF
    
    # Create example batch script
    cat > "$PROJECT_NAME/batch_process.sh" << 'BATCH_EOF'
#!/bin/bash
# Batch processing script for multiple genomes

# List of genome accessions to process
GENOMES=(
    # Add your genome accessions here, one per line
    # NC_001477  # Dengue virus 1
    # NC_001474  # Dengue virus 2
)

# Process each genome
for GENOME in "${GENOMES[@]}"; do
    echo "Processing $GENOME..."
    
    # Download
    cd genomes
    download_ncbi_genome "$GENOME"
    cd ..
    
    # Parse
    python3 step1_parse_viral_genome.py "genomes/$GENOME"
    mv "${GENOME}"*.tsv annotations/
    
    echo "Ready for manual curation: annotations/${GENOME}_no_polyprotein.tsv"
done
BATCH_EOF
    chmod +x "$PROJECT_NAME/batch_process.sh"
    
    echo ""
    echo "[SUCCESS] Project initialized: $PROJECT_NAME"
    echo ""
    echo "Project structure:"
    tree -L 2 "$PROJECT_NAME" 2>/dev/null || ls -la "$PROJECT_NAME"
    echo ""
    echo "Next steps:"
    echo "1. cd $PROJECT_NAME"
    echo "2. source /path/to/snpeff_configs/snpeff_current.sh"
    echo "3. Start processing genomes!"
}

# Add these functions to the help/usage display
if [ "$1" = "--help" ] || [ "$1" = "-h" ]; then
    echo ""
    echo "Additional VICAST Functions (after setup):"
    echo "  test_vicast_annotate [genome_id]  - Test the annotation pipeline"
    echo "  validate_vicast_setup              - Validate all dependencies"
    echo "  init_vicast_project [name]        - Initialize a new project"
    echo ""
    echo "Examples:"
    echo "  source $SOFTWARE_DIR/snpeff_configs/snpeff_current.sh"
    echo "  test_vicast_annotate              # Test with default genome"
    echo "  test_vicast_annotate NC_001477    # Test with specific genome"
    echo "  validate_vicast_setup              # Check all dependencies"
    echo "  init_vicast_project my_viruses    # Start new project"
fi

# Run main function
main
