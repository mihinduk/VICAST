# =============================================================================
# VICAST Docker Image
# Viral Cultured-virus Annotation and SnpEff Toolkit
# =============================================================================
#
# Build:
#   docker build -t vicast:latest .
#
# Run interactive shell:
#   docker run -it -v $(pwd):/data vicast:latest bash
#
# Run annotation pipeline:
#   docker run -v $(pwd):/data vicast:latest \
#     python /opt/vicast/vicast-annotate/step1_parse_viral_genome.py NC_001477
#
# =============================================================================

FROM mambaorg/micromamba:1.5-jammy

LABEL maintainer="Kathie Mihindukulasuriya, Scott Handley"
LABEL description="VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit"
LABEL version="2.2.0"

# Set environment variables
ENV VICAST_HOME=/opt/vicast
ENV SNPEFF_HOME=/opt/snpEff
ENV SNPEFF_JAR=${SNPEFF_HOME}/snpEff.jar
ENV SNPEFF_DATA=${SNPEFF_HOME}/data
ENV NCBI_EMAIL=vicast_docker@example.com
ENV VICAST_THREADS=4

USER root

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget \
    unzip \
    git \
    && rm -rf /var/lib/apt/lists/*

# Create vicast user
ARG MAMBA_USER=mambauser
ARG MAMBA_USER_ID=1000
ARG MAMBA_USER_GID=1000

# Copy environment file
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install conda environment
USER $MAMBA_USER
RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Switch back to root for SnpEff installation
USER root

# Download and install SnpEff
RUN mkdir -p /opt && \
    cd /opt && \
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 3 \
        https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip || \
    curl -L -o snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip -q snpEff_latest_core.zip && \
    rm snpEff_latest_core.zip && \
    chmod -R 755 /opt/snpEff

# Copy VICAST source code
COPY --chown=$MAMBA_USER:$MAMBA_USER . ${VICAST_HOME}

# Install VICAST Python package
USER $MAMBA_USER
RUN cd ${VICAST_HOME} && \
    /opt/conda/bin/pip install -e .

# Set working directory
WORKDIR /data

# Add VICAST to PATH
ENV PATH="${VICAST_HOME}/vicast-annotate:${VICAST_HOME}/vicast-analyze:${PATH}"

# Default command
CMD ["bash"]
