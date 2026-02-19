---
title: 'VICAST: Viral Cultured-virus Annotation and SnpEff Toolkit for Passage Study Analysis'
tags:
  - Python
  - virology
  - viral genomics
  - genome annotation
  - variant calling
  - passage studies
  - SnpEff
authors:
  - name: Kathie A. Mihindukulasuriya
    orcid: 0000-0001-9372-3758
    affiliation: 1
  - name: Luis Alberto Chica Cardenas
    orcid: 0009-0007-7408-9385
    affiliation: 1
  - name: Scott A. Handley
    orcid: 0000-0002-2143-6570
    corresponding: true
    affiliation: 1
affiliations:
  - name: Department of Pathology and Immunology, Washington University School of Medicine, St. Louis, MO 63110, USA
    index: 1
    ror: 01yc7t268
date: 10 February 2026
bibliography: paper.bib
---

# Summary

VICAST (Viral Cultured-virus Annotation and SnpEff Toolkit) is a comprehensive software suite for viral genome annotation and variant analysis, specifically designed for tissue culture passage studies. The toolkit addresses critical gaps in viral genomics workflows by combining semi-automated genome annotation with quality curation checkpoints and low-frequency variant calling optimized for viral populations. VICAST integrates with SnpEff to provide functional annotation of variants, supporting complex viral genome architectures including polyproteins and multi-segment viruses. The software is distributed as a Docker container for reproducible deployment and as a conda-based local installation, ensuring reproducibility across computing environments including high-performance computing clusters.

# Statement of Need

Cultured virus passage studies are fundamental to understanding viral evolution, attenuation, and host adaptation [@Duffy2008; @Stern2017]. However, analyzing genomic changes across passages requires two distinct capabilities that existing tools do not adequately address together: (1) accurate functional annotation of poorly-characterized or newly-sequenced viral genomes [@Lauber2017], and (2) detection and annotation of low-frequency variants that emerge during tissue culture adaptation [@Acevedo2014].

Current viral genome annotation tools face significant challenges with cultured viruses. Reference genomes from GenBank often contain incomplete or inconsistent annotations, particularly for polyproteins that are post-translationally cleaved into multiple functional proteins [@Lauber2017]. Automated annotation pipelines like VADR [@Schaffer2020] and VIGOR4 [@Wang2010; @Wang2012] work well for well-characterized virus families but struggle with novel or poorly-annotated genomes. Meanwhile, variant calling pipelines designed for clinical diagnostics (e.g., nf-core/viralrecon [@Patel2020]) focus on consensus sequences rather than the low-frequency variants (0.5-5% frequency) that are biologically meaningful in passage studies [@Acevedo2014].

VICAST bridges this gap by providing: (1) a semi-automated annotation curation workflow with manual quality checkpoints that ensures accurate functional annotation before variant analysis, (2) integration with SnpEff for consistent variant effect prediction, (3) low-frequency variant calling using lofreq [@Wilm2012] with validated filtering strategies, and (4) specialized support for complex viral architectures including segmented genomes and polyproteins. The target audience includes virologists conducting passage studies, researchers working with poorly-characterized viruses, and groups establishing new viral culture systems where accurate functional annotation is essential for interpreting evolutionary trajectories.

# State of the Field

Several tools address aspects of viral genome annotation and analysis, but none provide the integrated curation-to-variant workflow required for passage studies (Table 1).

**Viral annotation tools:** VADR [@Schaffer2020] provides automated validation and annotation for submissions to NCBI's GenBank, with pre-built models for well-characterized virus families. However, VADR's fixed models cannot accommodate novel viruses or incorporate manual curation, and it lacks variant calling capabilities. VIGOR4 [@Wang2010; @Wang2012] similarly provides automated annotation but is limited to influenza and other specific virus families. Neither tool integrates with variant effect prediction frameworks.

**Viral analysis pipelines:** nf-core/viralrecon [@Patel2020] and ViroProfiler [@Ru2023] provide end-to-end workflows for viral genome analysis including assembly and variant calling. However, these pipelines focus on clinical consensus sequences rather than low-frequency variants, lack annotation curation capabilities, and do not integrate variant functional annotation. MVP (Modular Viromics Pipeline) [@Coclet2024] offers flexible metagenomics analysis but is designed for environmental samples rather than cultured virus systems.

**VICAST's unique contributions:** Unlike existing tools, VICAST combines (1) flexible annotation pathways accommodating both well-annotated and novel genomes, (2) manual curation checkpoints ensuring annotation quality before downstream analysis, (3) SnpEff database integration for consistent variant effect prediction across genomes, (4) low-frequency variant calling (0.5-5% frequency range) optimized for detecting passage-associated changes, (5) specialized handling of polyproteins with mature peptide tracking, and (6) native multi-segment genome support for viruses like influenza. This integrated approach ensures that variants detected during passage are annotated with accurate functional context.

**Table 1: Comparison of viral genomics tools**

| Feature | VICAST | VADR | nf-core/viralrecon | ViroProfiler |
|---------|--------|------|-------------------|--------------|
| Annotation curation | ✓ | ✗ | ✗ | ✗ |
| Novel virus support | ✓ | Limited | ✗ | ✗ |
| Low-frequency variants | ✓ | ✗ | ✗ | Limited |
| Variant functional annotation | ✓ | ✗ | Limited | ✗ |
| Polyprotein handling | ✓ | Limited | ✗ | ✗ |
| Multi-segment genomes | ✓ | Limited | ✗ | ✗ |
| Passage study optimization | ✓ | ✗ | ✗ | ✗ |

# Software Design

VICAST consists of two integrated components reflecting the distinct phases of passage study analysis:

## VICAST-annotate: Genome Annotation Pipeline

VICAST-annotate implements four annotation pathways to accommodate the diverse annotation quality of viral genomes:

**Pathway 1** checks if genomes already exist in SnpEff databases, avoiding redundant work, but allowing overwriting if the user wants to update the annotation. **Pathway 2** processes well-annotated GenBank files by parsing CDS features and mat_peptides for polyproteins and integrating them into SnpEff with viral-optimized settings (5' and 3' UTR extensions, shifted reading frames). **Pathway 3** uses BLASTx homology searches against custom or public protein databases to annotate poorly-characterized genomes, producing curated TSV files for manual review. **Pathway 4** handles segmented viruses (e.g., influenza, rotavirus) by concatenating segments with unique identifiers and managing multi-FASTA references.

The critical design decision is **mandatory manual curation checkpoints** between automated steps. After parsing or BLAST annotation, researchers review TSV files to verify gene names, resolve overlapping features, and handle polyproteins before SnpEff database creation. This semi-automated approach balances efficiency with accuracy, acknowledging that viral genome annotation requires domain expertise that cannot be fully automated.

## VICAST-analyze: Variant Calling Pipeline

VICAST-analyze implements a nine-step workflow: (1) reference genome preparation, (2) read statistics calculation, (3) quality control with fastp [@Chen2018], (4) mapping with BWA-MEM [@Li2013] and variant calling with lofreq [@Wilm2012], (5) coverage depth profiling, (6) quality diagnostic report generation, (7) two-stage variant filtering (high-confidence ≥1% frequency, low-frequency 0.5-5%), (8) SnpEff functional annotation, and (9) TSV output generation with HGVS nomenclature.

Key architectural decisions include: **Two-stage filtering** distinguishes dominant variants (≥1% frequency, ≥200× depth, quality ≥1000) from low-frequency variants (0.5-5% frequency, ≥200× depth, quality ≥1000), enabling analysis of quasispecies dynamics. **Contamination screening** via de novo assembly (MEGAHIT [@Li2015]) and BLAST searches identifies non-target sequences before investing computational resources. **BAM read-level validation** confirms variant co-occurrence patterns to distinguish true haplotypes from alignment artifacts. The pipeline is structured as QC-first (steps 1-6) followed by annotation (steps 7-9), allowing manual review before computationally expensive variant annotation.

## Integration and Reproducibility

Both components share a unified configuration system supporting Docker, conda, and HPC environments. Docker containers ensure exact reproducibility for publications, while conda-based installations support development workflows. The modular design allows users to run annotation independently when adding new genomes, or full end-to-end analysis for viral culture experiments.

**Development evolution:** VICAST originated as two independent pipelines—VICAST-annotate (first released September 2025) and VICAST-analyze (first released July 2025)^[Original repositories: https://github.com/mihinduk/VICAST-annotate and https://github.com/mihinduk/VICAST-analyze]—developed separately to address annotation curation and variant calling respectively. In October 2025, these were unified into a single cohesive toolkit (https://github.com/mihinduk/VICAST) with integrated configuration, shared utilities, and coordinated workflows. This architectural evolution reflects the recognition that annotation quality directly impacts variant interpretation, necessitating tight integration between curation and analysis phases. The unified repository maintains 7 months of public development history (July 2025 - February 2026) with 292 commits across the three repositories, demonstrating sustained, iterative refinement based on real-world usage in passage study experiments.

# Research Impact Statement

VICAST's utility is demonstrated through comprehensive validation with three virus families representing distinct annotation challenges:

**SARS-CoV-2 (DRR878516):** Validates polyprotein handling by correctly annotating ORF1ab's 16 mature peptides (nsp1-nsp16), detecting D614G in the spike protein and properly mapping mutations to individual nonstructural proteins within the polyprotein. This addresses a common issue where polyprotein annotations fail to resolve variants to specific functional domains.

**Dengue virus 2 (SRR24480393):** Demonstrates standard genome annotation and low-frequency variant detection in a well-characterized flavivirus system. The pipeline successfully identifies variants in structural (E, prM, C) and nonstructural (NS1-NS5) proteins with correct functional annotations.

**Influenza A H1N1 (SRR36836026):** Validates multi-segment genome handling by processing all eight segments (PB2, PB1, PA, HA, NP, NA, M, NS) independently while maintaining consistent variant annotation. This capability is essential for segmented viruses where reassortment and segment-specific evolution must be tracked.

All validation datasets use publicly available SRA data and are distributed with the software as reproducible examples, enabling users to verify installation and understand expected outputs before analyzing their own viral culture experiments.

The comprehensive benchmark comparison with VADR (available at `benchmarks/vadr_comparison/`) demonstrates VICAST's advantages for passage study workflows: flexible annotation curation, low-frequency variant detection, and integrated functional annotation that VADR's validation-focused approach does not provide. While VADR excels at standardizing submissions for GenBank, VICAST optimizes for the iterative annotation-analysis cycles required in experimental virology.

Future research applications include multi-passage evolutionary trajectory analysis, identification of tissue culture adaptation signatures, and comparative analysis of attenuation strategies across virus families. The software has been actively used for internal passage studies (data not yet published), demonstrating real-world applicability beyond validation datasets.

# Generative AI Usage Disclosure

Development of VICAST was assisted by Claude (Anthropic, Sonnet 4.5) in the following capacities:

**Documentation:** Comprehensive user guides (8 guides, >3,000 lines) were drafted with AI assistance to ensure consistent formatting, complete coverage of installation scenarios (Docker/Conda/HPC), and clear workflow explanations. All content was reviewed and validated by the authors against actual software behavior.

**Testing framework:** Test suite structure and pytest configurations were developed with AI assistance. Test cases themselves were designed by the authors based on biological validation requirements.

**Code optimization:** Specific improvements to Docker containerization, configuration management, and error handling were implemented with AI suggestions. All code changes were reviewed, tested, and validated by the authors.

**Paper drafting:** Initial structure and technical descriptions for this manuscript were drafted with AI assistance, then substantially revised by the authors to ensure accurate representation of research context and biological significance.

No AI-generated content appears in the final software without human review, validation, and often substantial modification. All scientific decisions, architectural designs, and biological interpretations are the work of the named authors.

# Figures

![Overview of the VICAST pipeline. VICAST comprises two modules. VICAST-Annotate (left) prepares viral reference genomes for variant annotation through four pathways (pre-built database, GenBank parsing, BLASTx homology search, or segmented genome assembly), with a manual curation checkpoint before building a custom SnpEff database. VICAST-Analyze (right) processes paired-end sequencing reads through parallel tracks: read-level QC, alignment, and variant calling (left track) and de novo assembly with BLAST-based contamination screening against a curated database of 18,804 viral and microbial sequences (right track). Two manual QC checkpoints allow users to review contamination results and variant quality before proceeding to SnpEff annotation. Outputs include annotated variants in multiple formats (VCF, TSV, HTML) and an optional advanced analysis module for haplotype consensus generation and read-level co-occurrence validation. Tool names are shown in italics beneath each step.](figures/figure1.png)

# Acknowledgments

[PLACEHOLDER - Funding information to be provided]

We thank the Washington University Research Computing Facility (HTCF) for computational resources. We acknowledge the NCBI Sequence Read Archive and contributing researchers for making public viral sequencing datasets available for validation.

# References
