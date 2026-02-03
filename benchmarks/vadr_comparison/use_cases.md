# Use Case Guide: VICAST vs VADR

## Decision Tree

```
                    What is your primary goal?
                              │
          ┌───────────────────┼───────────────────┐
          │                   │                   │
          ▼                   ▼                   ▼
    GenBank              Variant             Novel Virus
    Submission           Analysis            Annotation
          │                   │                   │
          ▼                   ▼                   ▼
       VADR               VICAST              VICAST
                              │             (Pathway 3)
                              │
                    Is virus in VADR?
                    ┌─────────┴─────────┐
                   YES                  NO
                    │                   │
                    ▼                   ▼
              Use Both:            VICAST only
              VADR → VICAST
```

## Detailed Use Cases

### Use Case 1: Viral Passage Study

**Scenario**: You're studying how a virus adapts during serial passage in cell culture.

**Recommended Tool**: **VICAST**

**Why**:
- Native SnpEff integration for variant effect prediction
- Tracks synonymous vs. non-synonymous changes
- Manual curation ensures accurate gene boundaries
- Integrated variant calling pipeline

**Workflow**:
```bash
# 1. Annotate reference genome
python vicast_annotate.py NC_001477

# 2. Manual curation of TSV
# (edit NC_001477_no_polyprotein.tsv)

# 3. Build SnpEff database
python step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv

# 4. Run variant analysis
./run_vicast_analyze_full.sh reads_R1.fq reads_R2.fq NC_001477
```

---

### Use Case 2: Polyprotein Virus Analysis (Flavivirus, Picornavirus, etc.)

**Scenario**: Studying adaptation in Dengue, Zika, HCV, or other polyprotein-encoding viruses.

**Recommended Tool**: **VICAST**

**Why**:
- Specialized polyprotein → mature peptide annotation
- Variants annotated at individual protein level, not "polyprotein"
- Critical for interpreting NS3, NS5, and other functional domains
- SnpEff reports meaningful gene-level effects

**Workflow**:
```bash
# 1. Annotate with polyprotein skipping (default)
python vicast_annotate.py NC_001477  # Dengue

# 2. Review TSV - verify cleavage sites
# File shows individual mature peptides:
#   C, prM, E, NS1, NS2A, NS2B, NS3, NS4A, NS4B, NS5

# 3. Build SnpEff database
python step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv

# 4. Variant analysis - now shows meaningful annotations
./run_vicast_analyze_full.sh reads_R1.fq reads_R2.fq NC_001477
# Output: "p.Glu135Val in NS5" instead of "p.Glu2457Val in polyprotein"
```

**Result**: Variants mapped to functional domains (e.g., NS3 helicase, NS5 RdRp) rather than generic polyprotein positions.

---

### Use Case 3: GenBank Submission

**Scenario**: You've sequenced a novel virus isolate and want to submit to GenBank.

**Recommended Tool**: **VADR**

**Why**:
- Designed specifically for GenBank submission
- Generates required feature tables
- Validates against NCBI standards
- Detects assembly and annotation errors

**Workflow**:
```bash
# Run VADR validation and annotation
v-annotate.pl --mdir vadr-models novel_virus.fasta output_dir

# Review alerts and fix issues
# Submit to GenBank
```

---

### Use Case 4: Novel Virus Without Reference

**Scenario**: You've discovered a novel virus with no closely related sequences in databases.

**Recommended Tool**: **VICAST** (Pathway 3)

**Why**:
- BLASTx can find distant homologs
- Manual curation for novel features
- No pre-built model required
- Flexible annotation editing

**Workflow**:
```bash
# 1. BLASTx-based annotation
python step1_blastx_annotate.py novel_virus.fasta --blast-db viral_proteins

# 2. Extensive manual curation
# (review BLAST hits, define gene boundaries)

# 3. Build SnpEff database
python step2_add_to_snpeff.py novel_virus novel_virus_blastx.tsv
```

---

### Use Case 5: High-Throughput Surveillance

**Scenario**: Processing hundreds of viral sequences for surveillance.

**Recommended Tool**: **VADR** (with VICAST for selected samples)

**Why VADR**:
- Fully automated, no manual intervention
- Consistent, reproducible results
- Quality metrics for flagging issues

**Why add VICAST**:
- For samples requiring detailed variant analysis
- For samples with unusual features

**Workflow**:
```bash
# 1. Batch VADR processing
for fasta in samples/*.fasta; do
    v-annotate.pl --mdir models $fasta results/
done

# 2. Identify samples needing detailed analysis
# (based on VADR alerts or research interest)

# 3. VICAST for selected samples
python vicast_annotate.py sample_of_interest.fasta
```

---

### Use Case 6: Segmented Virus Analysis

**Scenario**: Analyzing influenza virus evolution across passages.

**Recommended Tool**: **VICAST**

**Why**:
- Native segmented virus support
- Unified database for all segments
- Variant calling across all segments
- Manual curation for reassortant analysis

**Workflow**:
```bash
# 1. Create combined database
python vicast_annotate_segmented.py influenza_h1n1 \
    --segments CY121680,CY121681,CY121682,CY121683,CY121684,CY121685,CY121686,CY121687 \
    --names PB2,PB1,PA,HA,NP,NA,M,NS

# 2. Run analysis for each passage
./run_vicast_analyze_full.sh p0_R1.fq p0_R2.fq influenza_h1n1
./run_vicast_analyze_full.sh p5_R1.fq p5_R2.fq influenza_h1n1
./run_vicast_analyze_full.sh p10_R1.fq p10_R2.fq influenza_h1n1

# 3. Compare variants across passages
```

---

### Use Case 7: Contamination Screening for Passage Studies

**Scenario**: You need to verify sample identity and check for contamination before analyzing viral passages.

**Recommended Tool**: **VICAST**

**Why**:
- Integrated contamination screening pipeline
- De novo assembly identifies unexpected sequences
- Multi-kingdom detection (viruses, mycoplasma, bacteria, fungi)
- Mapping statistics verify correct reference
- Presentation-ready HTML reports

**Workflow**:
```bash
# 1. Run contamination diagnostic
./vicast-analyze/viral_diagnostic.sh sample_R1.fq sample_R2.fq NC_001477 sample_001

# 2. Review diagnostic report
cat diagnostic_sample_001/sample_001_diagnostic_report.txt

# Key metrics to check:
#   - Mapping % to expected reference (should be >70%)
#   - Contaminating viruses (any unexpected viral contigs?)
#   - Mycoplasma contamination (common in cell culture)
#   - Duplication rate (high is OK for deep sequencing)

# 3. Open HTML report for presentation
open diagnostic_sample_001/sample_001_presentation_ready_report.html

# 4. If clean, proceed with variant analysis
./run_vicast_analyze_full.sh sample_R1.fq sample_R2.fq NC_001477
```

**Example Findings**:
```
MAPPING STATISTICS
==================
Deduplicated Mapping: 85.2% ← Good, correct virus

CONFIRMED VIRAL CONTAMINANTS: None ← Clean sample

CONTAMINATION SUMMARY:
  Virus: 1 contigs (expected reference)
  Mycoplasma: 0 contigs ← No mycoplasma
  Bacteria: 0 contigs
```

**When Contamination is Found**:
- Investigate source of contamination
- Consider re-culturing from clean stock
- Document in analysis notes
- May need to filter reads before variant calling

---

### Use Case 8: Quality Control Before Analysis

**Scenario**: Want to ensure genome assembly quality before in-depth analysis.

**Recommended Approach**: **VADR then VICAST**

**Workflow**:
```bash
# 1. VADR for QC
v-annotate.pl --mdir models assembly.fasta vadr_output/

# 2. Review VADR alerts
cat vadr_output/*.alt

# 3. If assembly passes QC, use VICAST for analysis
python vicast_annotate.py assembly.fasta

# 4. Curate and build database
python step2_add_to_snpeff.py assembly assembly_edited.tsv
```

---

### Use Case 9: Teaching and Training

**Scenario**: Teaching students about viral annotation.

**Recommended Tool**: **VICAST**

**Why**:
- Manual curation step teaches annotation concepts
- TSV format easy to understand and edit
- Step-by-step pathway makes process transparent
- Multiple pathways demonstrate different approaches

**Workflow**:
```bash
# Demonstrate different annotation sources
python step0_check_snpeff.py NC_045512  # Show SnpEff database
python step1_parse_viral_genome.py NC_001477  # Show NCBI parsing
python step1_blastx_annotate.py novel.fasta  # Show homology search

# Students edit TSV files to understand annotation
# Build databases and run variant analysis
```

---

## Quick Reference Table

| Scenario | Primary Tool | Add Secondary? |
|----------|--------------|----------------|
| Passage study | VICAST | Optional: VADR for QC |
| **Polyprotein viruses** (Flavi-, Picorna-, Corona-) | **VICAST** | No |
| **Contamination screening** | **VICAST** | No (VADR lacks this) |
| GenBank submission | VADR | No |
| Novel virus | VICAST (Pathway 3) | No |
| Surveillance | VADR | VICAST for select samples |
| Segmented viruses | VICAST | Optional: VADR for QC |
| Assembly quality control | VADR | VICAST for analysis |
| Teaching | VICAST | VADR for comparison |
| Quick annotation | VICAST (Pathway 1-2) | No |
| Standardized results | VADR | No |
| Custom curation | VICAST | No |
