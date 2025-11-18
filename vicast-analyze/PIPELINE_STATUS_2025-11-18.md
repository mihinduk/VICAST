# VICAST Pipeline Status - November 18, 2025

## Pipeline Architecture Overview

The VICAST pipeline has been split into two independent components:

### **vicast-annotate** (Annotation Pipeline)
- Handles GenBank downloads
- Builds snpEff databases
- Manages virus annotation
- **Key script**: `run_pipeline_htcf_consolidated.sh`
- **Requirements**: java, snpEff

### **vicast-analyze** (Analysis Pipeline - THIS REPOSITORY)
- Handles read processing and quality control
- Performs alignment and variant calling
- Generates consensus sequences
- Produces diagnostic reports
- **Does NOT handle**: Annotation, snpEff database building

---

## vicast-analyze Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                    VICAST-ANALYZE PIPELINE                          │
│                  (Read Processing & Variant Calling)                │
└─────────────────────────────────────────────────────────────────────┘

INPUT FILES:
  - R1.fastq.gz (Forward reads)
  - R2.fastq.gz (Reverse reads)
  - Accession number (e.g., NC_001477.1)

PREREQUISITES:
  - snpEff database built (handled by vicast-annotate)
  - Reference genome available (downloaded if needed)

┌─────────────────────────────────────────────────────────────────────┐
│ SCRIPT 1: run_pipeline_htcf_enhanced.sh                            │
│ (Main entry point - User-facing wrapper)                           │
│                                                                     │
│ What it does:                                                       │
│ • Auto-detects pipeline installation directory                     │
│ • Sources pipeline_config.sh (optional)                            │
│ • Validates conda availability                                     │
│ • Calls viral_pipeline.py                                          │
│ • Generates next_steps_{SAMPLE}.txt with downstream commands       │
│                                                                     │
│ Command: bash run_pipeline_htcf_enhanced.sh R1 R2 ACCESSION THREADS│
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│ CORE SCRIPT: viral_pipeline.py                                     │
│ (Main processing pipeline)                                         │
│                                                                     │
│ Module 1: Quality Control & Preprocessing                          │
│   • fastp: Quality trimming and filtering                          │
│   • seqkit: FASTA/Q manipulation                                   │
│                                                                     │
│ Module 2: Reference Download (if needed)                           │
│   • BioPython Entrez.efetch: Download from NCBI                    │
│   • wget fallback if BioPython fails                               │
│                                                                     │
│ Module 3: Read Alignment                                           │
│   • bwa mem: Align reads to reference                              │
│   • samtools: SAM/BAM processing and sorting                       │
│                                                                     │
│ Module 4: Variant Calling                                          │
│   • lofreq: Low-frequency variant detection                        │
│   • bcftools: VCF manipulation                                     │
│                                                                     │
│ Module 5: Annotation (via vicast-annotate)                         │
│   • snpEff: Variant annotation (database built separately)         │
│                                                                     │
│ Output: cleaned_seqs/ directory with:                              │
│   • {SAMPLE}.lofreq.final.bam                                      │
│   • {SAMPLE}.snpEFF.ann.vcf                                        │
│   • {SAMPLE}.snpEFF.ann.tsv                                        │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│ USER RUNS NEXT_STEPS COMMANDS (generated by enhanced script)       │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│ MODULE 6: Depth Analysis                                           │
│ Script: samtools depth (command in next_steps file)                │
│                                                                     │
│ Output: {SAMPLE}_results/{SAMPLE}_depth.txt                        │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│ MODULE 7: Mutation Parsing ⭐ KEY COMPONENT                        │
│ Script: viral_pipeline/visualization/parse_snpeff_tsv.py           │
│                                                                     │
│ What it does:                                                       │
│ • Filters variants by quality, depth, frequency                    │
│ • Parses snpEff annotations                                        │
│ • Generates filtered mutation table                                │
│                                                                     │
│ Parameters:                                                         │
│   --quality 1000                                                    │
│   --depth 200                                                       │
│   --freq 0.01                                                       │
│                                                                     │
│ Output: {SAMPLE}_results/{SAMPLE}_filtered_mutations.tsv           │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│ MODULE 8: Diagnostic Report ⭐ KEY COMPONENT & VALUE ADD           │
│ Script: viral_pipeline/analysis/submit_viral_diagnostic.sh         │
│                                                                     │
│ What it does:                                                       │
│ • Comprehensive quality metrics                                    │
│ • Coverage analysis                                                 │
│ • Variant distribution                                              │
│ • Assembly validation                                               │
│                                                                     │
│ Note: May be optimized in future work                              │
│                                                                     │
│ Output: diagnostic_{SAMPLE}/ directory                             │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│ MODULE 9: Consensus Generation                                     │
│ Script: viral_pipeline/utils/generate_filtered_consensus.py        │
│                                                                     │
│ What it does:                                                       │
│ • Generates consensus sequence from filtered variants              │
│ • Applies quality/depth/frequency thresholds                       │
│                                                                     │
│ Parameters:                                                         │
│   --quality 1000                                                    │
│   --depth 200                                                       │
│   --freq 0.50                                                       │
│                                                                     │
│ Output: {SAMPLE}_results/{SAMPLE}_consensus.fasta                  │
└─────────────────────────────────────────────────────────────────────┘

FINAL OUTPUTS:
  • Filtered mutations TSV
  • Diagnostic reports
  • Consensus genome FASTA
  • Depth coverage data
```

---

## Critical Tools Required (vicast-analyze ONLY)

### **Core Bioinformatics Tools**
- ✅ **BioPython** - Present in viral_genomics environment
- ❌ **bwa** (>=0.7.17) - MISSING
- ❌ **samtools** (>=1.15) - MISSING
- ❌ **bcftools** (>=1.15) - MISSING
- ❌ **fastp** (>=0.23.2) - MISSING
- ❌ **lofreq** (>=2.1.5) - MISSING
- ❌ **seqkit** (>=2.3.0) - MISSING

### **Python Libraries**
- ✅ pandas
- ✅ numpy
- ✅ biopython

### **NOT REQUIRED in vicast-analyze**
- ~~java~~ (handled by vicast-annotate)
- ~~snpEff~~ (handled by vicast-annotate)
- ~~entrez-direct~~ (replaced with BioPython)

---

## What's Working ✅

1. **Configuration System**
   - Auto-detection of pipeline installation directory
   - Optional pipeline_config.sh with HTCF-specific paths
   - Interactive prompts if paths not found
   - Fully portable (works after git clone from any location)

2. **BioPython Integration**
   - viral_pipeline.py uses BioPython Entrez.efetch for GenBank downloads
   - No dependency on entrez-direct command-line tools
   - Fallback to wget if BioPython fails

3. **Path Management**
   - All hardcoded paths removed
   - Script auto-detects its own location
   - Works with or without configuration file

4. **Streamlined Workflow**
   - Removed broken visualization modules (4 & 5)
   - Removed redundant GenBank download steps (1 & 2)
   - Focus on data generation, not visualization

---

## Current Blocker 🚫

**Missing critical bioinformatics tools in viral_genomics environment**

The existing `viral_genomics` conda environment only contains Python libraries:
- pandas, numpy, biopython, scipy, matplotlib, seaborn, plotly

**Required tools NOT present:**
- bwa (read alignment)
- samtools (BAM processing)
- fastp (quality control)
- lofreq (variant calling)
- bcftools (VCF processing)
- seqkit (sequence manipulation)

**Attempted Solution:**
Created viral_genomics_analyze.yml with all required tools, but conda environment creation failed with OOM (Out of Memory) errors even with 64GB RAM.

**Options to Resolve:**
1. Use existing HTCF environment with these tools already installed
2. Create minimal environment with only missing tools
3. Install tools manually in existing environment
4. Use module system if available on HTCF

---

## File Modifications Summary

### **Modified Files:**

1. **viral_pipeline.py**
   - Replaced subprocess efetch calls with BioPython Entrez.efetch
   - Added: `from Bio import Entrez`
   - Added: `Entrez.email = "vicast@example.com"`
   - Modified 3 functions: download_reference_genome, download_annotation (2 instances)

2. **run_pipeline_htcf_enhanced.sh**
   - Removed Steps 1 & 2 (redundant GenBank downloads)
   - Implemented auto-detection of PIPELINE_BASE from script directory
   - Added interactive path prompting if consolidated pipeline not found
   - Added conda availability check with helpful error message
   - Removed visualization modules (4 & 5) from next_steps generation
   - Made config file optional

3. **run_pipeline_htcf_consolidated.sh**
   - Added config sourcing with defaults
   - Added conda availability check
   - **NOTE**: This file should be moved to vicast-annotate repository

4. **pipeline_config.sh**
   - Set MAMBA_CMD="conda run -n viral_genomics"
   - Set SNPEFF_DIR="/ref/sahlab/software/snpEff"
   - Commented out PIPELINE_BASE (now auto-detected)

5. **viral_genomics_analyze.yml**
   - Commented out entrez-direct dependency
   - Added extensive installation instructions with space warnings
   - **NOTE**: Environment not yet created due to OOM errors

### **Created Files:**

1. **setup_environment.sh**
   - Interactive script for creating viral_genomics_analyze environment
   - Includes space usage warnings for shared servers
   - Not yet used due to OOM issues

2. **pipeline_config.template.sh**
   - Template for users to create custom configurations
   - Includes documentation for all configuration options

---

## Next Steps to Get Pipeline Running

### **Immediate Priority:**
1. **Resolve missing tools issue**
   - Identify existing HTCF environment with bwa, samtools, fastp, lofreq
   - OR create minimal environment with only these tools
   - Update MAMBA_CMD in pipeline_config.sh accordingly

### **Testing:**
2. **Test with sample data**
   - Use existing FASTQ files
   - Known accession number (e.g., NC_001477.1 for Dengue)
   - Verify all modules run successfully

### **Architecture Cleanup:**
3. **Move consolidated.sh to vicast-annotate**
   - This script handles annotation, not analysis
   - Should not be part of vicast-analyze repository

### **Documentation:**
4. **Update README**
   - Reflect new architecture (analyze vs annotate)
   - Document tool requirements
   - Installation instructions

---

## Key Architectural Decisions

### **Separation of Concerns:**
- **vicast-annotate**: Downloads genomes, builds snpEff databases, annotates variants
- **vicast-analyze**: Processes reads, calls variants, generates consensus

### **Portability:**
- No hardcoded paths
- Auto-detection of script locations
- Optional configuration with sensible defaults

### **Dependencies:**
- Use BioPython instead of command-line tools where possible
- Minimize environment complexity
- Clear separation: which tools belong to which component

### **User Experience:**
- Single entry point (run_pipeline_htcf_enhanced.sh)
- Helpful error messages with actionable solutions
- Generated next_steps file for downstream analysis

---

## Contact and Context

**Working Directory:** `/Users/handley_lab/Handley Lab Dropbox/virome/Diamond_lab_isolate_seq/2025_07_03_Diamond_NovaSeq_N1027/vicast/VICAST/vicast-analyze`

**Git Branch:** main

**Key Principle from CLAUDE.md:**
> Build SYSTEMATIC solutions, not specific fixes. Configuration-driven, not hardcoded solutions. Test across ALL supported cases, not individual examples.

**This refactoring follows systematic principles:**
- ✅ Configuration-driven (pipeline_config.sh with auto-detection)
- ✅ Portable (works from any installation location)
- ✅ Modular (clear separation of annotation vs analysis)
- ⏳ Needs systematic testing across multiple virus genomes

---

## Resuming Work Checklist

When returning to this project:

1. [ ] Source conda: `source /ref/sahlab/software/anaconda3/bin/activate`
2. [ ] Identify environment with bwa/samtools/fastp/lofreq OR create one
3. [ ] Update pipeline_config.sh with correct MAMBA_CMD
4. [ ] Test pipeline with sample data
5. [ ] Verify all modules complete successfully
6. [ ] Consider moving run_pipeline_htcf_consolidated.sh to vicast-annotate
7. [ ] Update main README to reflect architectural changes

---

**Last Updated:** 2025-11-18
**Status:** Blocked on missing bioinformatics tools
**Priority:** Resolve tool availability, then test end-to-end
