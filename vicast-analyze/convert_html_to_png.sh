#!/bin/bash
# =============================================================================
# Convert HTML Report to PNG Image
# =============================================================================
# Converts VICAST diagnostic HTML reports to PNG images for publication
#
# Usage:
#   ./convert_html_to_png.sh <input.html> [output.png]
#   ./convert_html_to_png.sh diagnostic_report.html
#   ./convert_html_to_png.sh report.html figure.png
#
# Requirements:
#   - wkhtmltoimage (installed via conda or system package)
# =============================================================================

set -e

# Check arguments
if [ $# -lt 1 ]; then
    echo "Usage: $0 <input.html> [output.png]"
    echo ""
    echo "Examples:"
    echo "  $0 diagnostic_SRR5992153_presentation_ready_report.html"
    echo "  $0 report.html custom_output.png"
    echo ""
    echo "Converts HTML diagnostic reports to PNG images."
    echo "Output filename defaults to input with .png extension if not specified."
    exit 1
fi

INPUT_HTML="$1"

# Check if input file exists
if [ ! -f "$INPUT_HTML" ]; then
    echo "Error: Input file not found: $INPUT_HTML"
    exit 1
fi

# Determine output filename
if [ $# -ge 2 ]; then
    OUTPUT_PNG="$2"
else
    # Replace .html extension with .png
    OUTPUT_PNG="${INPUT_HTML%.html}.png"
fi

# Check if wkhtmltoimage is available
if ! command -v wkhtmltoimage &> /dev/null; then
    echo "Error: wkhtmltoimage not found"
    echo ""
    echo "Please install wkhtmltoimage:"
    echo "  conda: conda install -c conda-forge wkhtmltopdf"
    echo "  macOS: brew install wkhtmltopdf"
    echo "  Ubuntu: sudo apt-get install wkhtmltopdf"
    exit 1
fi

echo "Converting HTML to PNG..."
echo "  Input:  $INPUT_HTML"
echo "  Output: $OUTPUT_PNG"

# Convert HTML to PNG
# --quality: Image quality (0-100, default 94)
# --width: Page width in pixels (1200px = good for figures)
# --enable-local-file-access: Allow loading local resources
# The tool automatically handles scrolling/long pages
wkhtmltoimage \
    --quality 95 \
    --width 1200 \
    --enable-local-file-access \
    "$INPUT_HTML" \
    "$OUTPUT_PNG"

# Check if conversion succeeded
if [ -f "$OUTPUT_PNG" ]; then
    FILE_SIZE=$(ls -lh "$OUTPUT_PNG" | awk '{print $5}')
    echo "✓ Conversion complete: $OUTPUT_PNG ($FILE_SIZE)"
else
    echo "✗ Conversion failed"
    exit 1
fi
