# VICAST User Guides

Complete documentation for all VICAST features and workflows.

---

## 🔄 VICAST Workflow Overview

```mermaid
flowchart TB
    subgraph "VICAST-ANNOTATE: Genome Preparation"
        A1[Reference Genome<br/>GenBank/FASTA] --> A2{Annotation<br/>Pathway}
        A2 -->|Pathway 1| A3[Check SnpEff<br/>Already exists?]
        A2 -->|Pathway 2| A4[Parse GenBank<br/>Standard annotation]
        A2 -->|Pathway 3| A5[BLASTx<br/>Homology search]
        A2 -->|Pathway 4| A6[Segmented Virus<br/>Combine segments]

        A3 --> A7[Ready to use!]
        A4 --> A8[Manual Curation<br/>QC checkpoint]
        A5 --> A8
        A6 --> A8
        A8 --> A9[Build SnpEff<br/>Database]
        A9 --> A10[(Custom Viral<br/>Database)]
        A7 --> A10
    end

    subgraph "VICAST-ANALYZE: Variant Calling Pipeline"
        B1[Raw Reads<br/>FASTQ] --> B2[Quality Control<br/>fastp]
        B2 --> B3[Read Alignment<br/>bwa mem]
        B2 --> B4[De novo Assembly<br/>MEGAHIT]
        B4 --> B5[Contamination Screening<br/>BLAST contigs]
        B3 --> B6[Variant Calling<br/>lofreq]
        B6 --> B7[Manual Filtering<br/>QC checkpoint]
        B7 --> B8[SnpEff Annotation<br/>& Effect Prediction]
        A10 -.Database input.-> B8
        B8 --> B9[Annotated Variants<br/>VCF + TSV + Reports]
    end

    subgraph "Advanced Analysis: Population Structure"
        C1[Variants + BAM] --> C2[Frequency-Stratified<br/>Analysis]
        C2 --> C3[Haplotype Consensus<br/>Generation]
        C3 --> C4{Validate<br/>Linkage?}
        C4 -->|Yes| C5[BAM Co-Occurrence<br/>Read-level evidence]
        C4 -->|No| C6[Dominant Consensus<br/>Alternative consensi]
        C5 --> C7[Validated Haplotypes<br/>+ Evidence metrics]
        C6 --> C7
    end

    B9 --> C1
    B5 --> C8[Contamination<br/>Report]

    style A8 fill:#fff4a3
    style B7 fill:#fff4a3
    style A10 fill:#e6b3ff
    style C8 fill:#ffcccc
    style C7 fill:#b3ffb3
    style B9 fill:#b3ffb3
```

**Legend:**
- 🔵 **Blue boxes** = Input data
- 🟠 **Orange boxes** = Processing steps
- 🟡 **Yellow diamonds** = Manual QC checkpoints
- 🟣 **Purple** = Database
- 🟢 **Green** = Final outputs
- 🔴 **Red** = Quality control reports

---

## 📚 Quick Navigation

### Getting Started
- **[Getting Started Guide](GETTING_STARTED.md)** - Installation, setup, and first analysis

### Core Workflows
- **[VICAST-Annotate Guide](VICAST_ANNOTATE_GUIDE.md)** - Annotate viral genomes for SnpEff
- **[VICAST-Analyze Guide](VICAST_ANALYZE_GUIDE.md)** - Complete variant calling pipeline
- **[Contamination Screening Guide](CONTAMINATION_SCREENING_GUIDE.md)** - De novo assembly + BLAST pipeline

### Advanced Features
- **[Haplotype Consensus Guide](HAPLOTYPE_CONSENSUS_GUIDE.md)** - Generate frequency-stratified consensus genomes
- **[BAM Co-Occurrence Guide](BAM_COOCCURRENCE_GUIDE.md)** - Read-level validation of variant linkage

### Reference
- **[Parameter Reference](PARAMETER_REFERENCE.md)** - Complete parameter documentation
- **[Troubleshooting Guide](TROUBLESHOOTING.md)** - Common issues and solutions

---

## By Use Case

### 🧬 New to VICAST?
1. [Getting Started Guide](GETTING_STARTED.md)
2. [VICAST-Analyze Guide](VICAST_ANALYZE_GUIDE.md)
3. [Contamination Screening Guide](CONTAMINATION_SCREENING_GUIDE.md)

### 📊 Passage Studies
1. [VICAST-Analyze Guide](VICAST_ANALYZE_GUIDE.md)
2. [Haplotype Consensus Guide](HAPLOTYPE_CONSENSUS_GUIDE.md)
3. [BAM Co-Occurrence Guide](BAM_COOCCURRENCE_GUIDE.md) - for validation

### 🆕 Adding New Virus
1. [VICAST-Annotate Guide](VICAST_ANNOTATE_GUIDE.md)
2. [VICAST-Analyze Guide](VICAST_ANALYZE_GUIDE.md)

### 🔬 Publication-Quality Analysis
1. [Contamination Screening Guide](CONTAMINATION_SCREENING_GUIDE.md) - critical for publications
2. [Haplotype Consensus Guide](HAPLOTYPE_CONSENSUS_GUIDE.md) - population structure
3. [BAM Co-Occurrence Guide](BAM_COOCCURRENCE_GUIDE.md) - validate claims

---

## Guide Descriptions

### Core Workflows

#### Getting Started Guide
Complete installation and setup instructions for VICAST on different systems (local, HPC, Docker). Includes quickstart example and verification steps.

#### VICAST-Annotate Guide (VICAST-Annotate)
How to prepare viral genomes for variant annotation:
- Download genomes from NCBI
- Parse GFF/GenBank files
- Build SnpEff databases
- Validate annotations
- Troubleshoot common issues

#### VICAST-Analyze Guide (VICAST-Analyze)
Complete variant calling pipeline from FASTQ to annotated VCF:
- Quality control (fastp)
- Read alignment (bwa)
- Variant calling (lofreq)
- Variant annotation (SnpEff)
- Filtering and reporting

#### Contamination Screening Guide
De novo assembly and BLAST-based contamination detection:
- When to use contamination screening
- De novo assembly with megahit
- BLAST against contamination databases
- Interpreting results
- Reporting for publications

### Advanced Features

#### Haplotype Consensus Guide
Generate frequency-stratified consensus genomes:
- Understanding variant frequencies
- Parameter selection by virus type
- Biological interpretation
- Publication language
- Limitations and caveats

#### BAM Co-Occurrence Guide
Validate variant co-occurrence with read-level evidence:
- When direct evidence is needed
- Running co-occurrence analysis
- Interpreting results
- Integration with haplotype consensus
- Strengthening publication claims

---

## Quick Reference Tables

### Command Cheat Sheet

| Task | Command | Guide |
|------|---------|-------|
| Install VICAST | `pip install -e .` | [Getting Started](GETTING_STARTED.md) |
| Annotate genome | `python step1_parse_viral_genome.py NC_001477` | [Annotation](VICAST_ANNOTATE_GUIDE.md) |
| Run variant calling | `bash run_pipeline.sh R1.fq.gz R2.fq.gz NC_001477` | [Variant Calling](VICAST_ANALYZE_GUIDE.md) |
| Screen contamination | `bash viral_diagnostic.sh reads.fq.gz` | [Contamination](CONTAMINATION_SCREENING_GUIDE.md) |
| Generate consensus | `python generate_realistic_haplotype_consensus.py` | [Haplotype](HAPLOTYPE_CONSENSUS_GUIDE.md) |
| Check co-occurrence | `python check_read_cooccurrence.py` | [BAM](BAM_COOCCURRENCE_GUIDE.md) |

### Typical Workflows

**Complete Analysis Pipeline:**
```bash
# 1. Annotate genome (once per virus)
python vicast-annotate/step1_parse_viral_genome.py NC_001477

# 2. Run variant calling
cd vicast-analyze
bash run_pipeline_htcf_consolidated.sh R1.fastq.gz R2.fastq.gz NC_001477 4

# 3. Screen for contamination
bash viral_diagnostic.sh sample_R1.fastq.gz sample_R2.fastq.gz

# 4. Generate consensus genomes
python generate_realistic_haplotype_consensus.py \
    --vcf sample_variants.vcf \
    --reference NC_001477.fasta \
    --accession NC_001477

# 5. Validate with BAM (optional)
python check_read_cooccurrence.py \
    --bam sample_aligned.bam \
    --vcf sample_variants.vcf \
    --output cooccurrence.tsv
```

**Quick Passage Study:**
```bash
# For each passage, run variant calling
for passage in P0 P5 P10; do
    bash run_pipeline.sh ${passage}_R1.fq.gz ${passage}_R2.fq.gz NC_001477
done

# Generate consensus for each
for passage in P0 P5 P10; do
    python generate_realistic_haplotype_consensus.py \
        --vcf ${passage}_variants.vcf \
        --reference NC_001477.fasta \
        --accession NC_001477 \
        --output-prefix ${passage}
done

# Compare frequency changes across passages
```

---

## Support & Troubleshooting

- **Common Issues:** See [Troubleshooting Guide](TROUBLESHOOTING.md)
- **Bug Reports:** https://github.com/mihinduk/VICAST/issues
- **Questions:** Check existing issues or create new one

---

## Contributing to Documentation

Found an error or have a suggestion? Please:
1. Open an issue on GitHub
2. Submit a pull request with corrections
3. Share your use case for new guide ideas

---

**Last Updated:** 2026-02-05
**VICAST Version:** 2.2.0
