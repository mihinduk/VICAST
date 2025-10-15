#!/bin/bash
# snpEff Environment Setup
# Source this file to use snpEff with Java 21

export JAVA_HOME="/ref/sahlab/software/jdk-21.0.5+11"
export SNPEFF_HOME="/home/mihindu/software/snpEff"
export PATH="$JAVA_HOME/bin:$PATH"

# Function to run snpEff
snpeff() {
    "$JAVA_HOME/bin/java" -Xmx4g -jar "$SNPEFF_HOME/snpEff.jar" "$@"
}

# Function to run snpSift
snpsift() {
    "$JAVA_HOME/bin/java" -Xmx4g -jar "$SNPEFF_HOME/SnpSift.jar" "$@"
}

echo "snpEff environment loaded!"
echo "  Java: $JAVA_HOME"
echo "  snpEff: $SNPEFF_HOME"
echo ""
echo "Usage:"
echo "  snpeff -version    # Run snpEff"
echo "  snpsift -h         # Run SnpSift"
