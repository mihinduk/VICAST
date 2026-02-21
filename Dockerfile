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
# SnpEff installed via conda
ENV SNPEFF_HOME=/opt/conda/share/snpeff-5.4.0a-0
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

# SnpEff is already installed via conda in environment.yml
# Create SnpEff data directory and set permissions
USER root
RUN mkdir -p ${SNPEFF_DATA} && \
    chown -R $MAMBA_USER:$MAMBA_USER ${SNPEFF_HOME} && \
    chmod -R 755 ${SNPEFF_HOME}

# Copy VICAST source code
COPY --chown=$MAMBA_USER:$MAMBA_USER . ${VICAST_HOME}

# Install VICAST Python package
USER $MAMBA_USER
RUN cd ${VICAST_HOME} && \
    /opt/conda/bin/pip install -e .

# Make all Python scripts executable
USER root
RUN find ${VICAST_HOME}/vicast-annotate -name "*.py" -exec chmod +x {} \; && \
    find ${VICAST_HOME}/vicast-analyze -name "*.py" -exec chmod +x {} \;

# Create wrapper scripts for easy access
RUN printf '#!/bin/bash\npython /opt/vicast/vicast-annotate/step1_parse_viral_genome.py "$@"\n' > /usr/local/bin/step1_parse_viral_genome.py && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-annotate/step1_blastx_annotate.py "$@"\n' > /usr/local/bin/step1_blastx_annotate.py && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-annotate/step2_add_to_snpeff.py "$@"\n' > /usr/local/bin/step2_add_to_snpeff.py && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-annotate/step0_check_snpeff.py "$@"\n' > /usr/local/bin/step0_check_snpeff.py && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-annotate/vicast_annotate_segmented.py "$@"\n' > /usr/local/bin/vicast_annotate_segmented.py && \
    chmod +x /usr/local/bin/step*.py && \
    chmod +x /usr/local/bin/vicast_annotate_segmented.py

# Create wrappers for analyze scripts
RUN printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/run_vicast_analyze_full.sh "$@"\n' > /usr/local/bin/run_vicast_analyze_full.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/run_vicast_analyze_qc_only.sh "$@"\n' > /usr/local/bin/run_vicast_analyze_qc_only.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/run_vicast_analyze_annotate_only.sh "$@"\n' > /usr/local/bin/run_vicast_analyze_annotate_only.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/download_sra_data.sh "$@"\n' > /usr/local/bin/download_sra_data.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/install_prebuilt_database.sh "$@"\n' > /usr/local/bin/install_prebuilt_database.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/setup_blast_db.sh "$@"\n' > /usr/local/bin/setup_blast_db.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/build_contamination_db.sh "$@"\n' > /usr/local/bin/build_contamination_db.sh && \
    chmod +x /usr/local/bin/run_vicast_analyze*.sh && \
    chmod +x /usr/local/bin/download_sra_data.sh && \
    chmod +x /usr/local/bin/install_prebuilt_database.sh && \
    chmod +x /usr/local/bin/setup_blast_db.sh && \
    chmod +x /usr/local/bin/build_contamination_db.sh && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-analyze/vcf_to_tsv.py "$@"\n' > /usr/local/bin/vcf_to_tsv.py && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/run_vicast_analyze_postprocess.sh "$@"\n' > /usr/local/bin/run_vicast_analyze_postprocess.sh && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-analyze/parse_snpeff_tsv.py "$@"\n' > /usr/local/bin/parse_snpeff_tsv.py && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-analyze/generate_realistic_haplotype_consensus.py "$@"\n' > /usr/local/bin/generate_realistic_haplotype_consensus.py && \
    printf '#!/bin/bash\npython /opt/vicast/vicast-analyze/parse_blast_results.py "$@"\n' > /usr/local/bin/parse_blast_results.py && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/test_blast_diagnostic.sh "$@"\n' > /usr/local/bin/test_blast_diagnostic.sh && \
    printf '#!/bin/bash\nbash /opt/vicast/vicast-analyze/extend_blast_db.sh "$@"\n' > /usr/local/bin/extend_blast_db.sh && \
    chmod +x /usr/local/bin/vcf_to_tsv.py && \
    chmod +x /usr/local/bin/run_vicast_analyze_postprocess.sh && \
    chmod +x /usr/local/bin/parse_snpeff_tsv.py && \
    chmod +x /usr/local/bin/generate_realistic_haplotype_consensus.py && \
    chmod +x /usr/local/bin/parse_blast_results.py && \
    chmod +x /usr/local/bin/test_blast_diagnostic.sh && \
    chmod +x /usr/local/bin/extend_blast_db.sh

# Create snpeff wrapper script
RUN printf '#!/bin/bash\njava -jar %s "$@"\n' "${SNPEFF_JAR}" > /usr/local/bin/snpeff && \
    chmod +x /usr/local/bin/snpeff

# Create helpful MOTD
RUN echo '#!/bin/bash\ncat << "EOF"\n\n╔══════════════════════════════════════════════════════════════╗\n║                    VICAST Docker v2.3.0                      ║\n║        Viral Cultured-virus Annotation & SnpEff Toolkit      ║\n╚══════════════════════════════════════════════════════════════╝\n\nPre-built Genomes (Ready to Use):\n  ✓ DENV-2      (NC_001474.2)\n  ✓ SARS-CoV-2  (NC_045512.2)\n  ✓ Ebola       (NC_002549.1)\n\nQuick Start - Pre-built Genome:\n  download_sra_data.sh SRR5992153\n  run_vicast_analyze_full.sh SRR5992153_1.fastq.gz SRR5992153_2.fastq.gz NC_001474.2 8\n\nCustom Genome Workflow:\n  step1_parse_viral_genome.py AF252854.1  # Download & parse\n  step2_add_to_snpeff.py AF252854.1 AF252854.1.tsv  # Add to SnpEff\n  run_vicast_analyze_full.sh R1.fq.gz R2.fq.gz AF252854.1 8\n\nContamination Screening:\n  ✓ BLAST database: 18,804 sequences (viruses + lab contaminants)\n  ✓ UniVec vector filtering: ~3,058 cloning vectors (automatic)\n  Run: setup_blast_db.sh --info\n\nMounted Locations:\n  /data              → Your working directory\n  /opt/vicast        → VICAST installation\n  $SNPEFF_DATA_CUSTOM → Custom databases (if mounted)\n\nDocumentation: /opt/vicast/README.md\n\nEOF' > /etc/profile.d/vicast_motd.sh && \
    chmod +x /etc/profile.d/vicast_motd.sh

# Create conda symlink to micromamba for compatibility
RUN ln -s /usr/bin/micromamba /usr/local/bin/conda

# Initialize micromamba for all users
RUN micromamba shell init --shell=bash --root-prefix=/opt/conda && \
    echo 'eval "$(micromamba shell hook --shell bash)"' >> /etc/bash.bashrc && \
    echo 'micromamba activate base' >> /etc/bash.bashrc

# =============================================================================
# Database Architecture: Hybrid Approach
# - Built-in databases: Pre-installed common genomes (fast startup)
# - Custom databases: User-mounted directory for curated genomes (persistent)
# =============================================================================

# Create directory for custom databases (user will mount here)
USER root
RUN mkdir -p /opt/vicast/snpeff_data_custom && \
    chown -R $MAMBA_USER:$MAMBA_USER /opt/vicast/snpeff_data_custom

# Pre-build common viral genomes (DISABLED - install at runtime instead)
# Users can run: install_prebuilt_database.sh --install NC_001474.2
# Commented out because manifest download fails during Docker build
# USER $MAMBA_USER
# RUN cd ${VICAST_HOME}/vicast-analyze && \
#     bash install_prebuilt_database.sh --install NC_001474.2 && \
#     echo "Pre-built database: DENV-2 (NC_001474.2)" && \
#     bash install_prebuilt_database.sh --install NC_045512.2 && \
#     echo "Pre-built database: SARS-CoV-2 (NC_045512.2)"

# Set environment variables for database locations
# SNPEFF_DATA_BUILTIN: Read-only pre-built databases in image (reference only)
# SNPEFF_DATA_CUSTOM: User-mounted writable location for custom genomes
# SNPEFF_DATA: Points to custom directory (SnpEff doesn't support colon-separated paths)
ENV SNPEFF_DATA_BUILTIN=${SNPEFF_DATA}
ENV SNPEFF_DATA_CUSTOM=/opt/vicast/snpeff_data_custom
ENV SNPEFF_DATA=${SNPEFF_DATA_CUSTOM}

# Fix SnpEff config to use absolute path for data directory
# Make config writable so users can add genome entries when running with --user flag
USER root
RUN sed -i 's|^data\.dir = \./data/$|data.dir = /opt/vicast/snpeff_data_custom|' ${SNPEFF_HOME}/snpEff.config && \
    chmod 666 ${SNPEFF_HOME}/snpEff.config && \
    grep "^data.dir" ${SNPEFF_HOME}/snpEff.config

# Create writable config location (for when running with --user flag)
ENV SNPEFF_CONFIG_BUILTIN=${SNPEFF_HOME}/snpEff.config
ENV SNPEFF_CONFIG_CUSTOM=/opt/vicast/snpeff_data_custom/snpEff.config

# =============================================================================
# BLAST Contamination Screening Database
# Pre-built database with 18,804 sequences: RefSeq viral genomes + common
# lab contaminants (E. coli, Pseudomonas, Staph, Mycoplasma, Candida, etc.)
# =============================================================================

ENV BLAST_DB_DIR=/opt/vicast/blast_db
ENV BLAST_DB=${BLAST_DB_DIR}/vicast_combined

USER root
RUN mkdir -p ${BLAST_DB_DIR}/user_extensions && \
    chown -R $MAMBA_USER:$MAMBA_USER ${BLAST_DB_DIR} && \
    chmod -R 775 ${BLAST_DB_DIR}

# Download pre-built contamination database from GitHub Release
# If download fails (e.g., restricted Docker network), place the tarball in
# the build context before building:
#   curl -fL -o microbial_contaminants.tar.gz \
#     https://github.com/mihinduk/VICAST/releases/download/blast-db-v1.0/microbial_contaminants.tar.gz
#   docker build -t vicast:latest .
# Or install at runtime: docker run --rm -v blastdb:/opt/vicast/blast_db vicast setup_blast_db.sh
USER $MAMBA_USER
RUN bash ${VICAST_HOME}/vicast-analyze/setup_blast_db.sh ${BLAST_DB_DIR} || \
    echo "WARNING: BLAST DB download failed during build - run setup_blast_db.sh at runtime"

# =============================================================================
# UniVec Cloning Vector Database
# NCBI UniVec (~3,058 vector sequences) for pre-mapping read filtering.
# Removes cloning vector artifacts (e.g., NotI sites from infectious clones).
# Used automatically by Step 4 (map_and_call_variants) before BWA alignment.
# =============================================================================
RUN curl -fL -o ${BLAST_DB_DIR}/cloning_vectors.fasta \
        "https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec" && \
    bwa index ${BLAST_DB_DIR}/cloning_vectors.fasta && \
    echo "UniVec database installed: $(grep -c '^>' ${BLAST_DB_DIR}/cloning_vectors.fasta) vectors" || \
    echo "WARNING: UniVec download failed - vector filtering will be skipped at runtime"

# Set working directory
WORKDIR /data

# Add VICAST to PATH
ENV PATH="${VICAST_HOME}/vicast-annotate:${VICAST_HOME}/vicast-analyze:${PATH}"

# Set friendly prompt for when running with --user flag
ENV PS1="(vicast) \\w$ "

# Custom entrypoint that suppresses micromamba warnings with --user
USER root
RUN printf '#!/bin/bash\neval "$(micromamba shell hook --shell bash)" 2>/dev/null\nmicromamba activate base 2>/dev/null\nexec "$@"\n' > /usr/local/bin/vicast_entrypoint.sh && \
    chmod +x /usr/local/bin/vicast_entrypoint.sh
ENTRYPOINT ["/usr/local/bin/vicast_entrypoint.sh"]

# Switch back to regular user for runtime
USER $MAMBA_USER

# Show welcome message on container start
ENV BASH_ENV=/etc/profile.d/vicast_motd.sh

# Default command
CMD ["bash"]
