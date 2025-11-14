# Pipeline Simplification Instructions

## Files Created

1. **pipeline_config.template.sh** - Template configuration file with all paths
2. **viral_genomics_analyze.yml** - Conda environment specification

## Setup Steps

### 1. Create the Conda Environment

```bash
# Navigate to where you saved the YAML file
cd /path/to/your/pipeline

# Create the environment
mamba env create -f viral_genomics_analyze.yml

# Activate it
mamba activate viral_genomics_analyze

# Verify installation
bwa
samtools --version
python --version
efetch -help
```

### 2. Configure Your Pipeline Paths

```bash
# Copy the template to create your config
cp pipeline_config.template.sh pipeline_config.sh

# Edit with your actual paths
nano pipeline_config.sh  # or use your preferred editor
```

**Update these key paths in pipeline_config.sh:**
- `MAMBA_CMD` → Your mamba installation path
- `SNPEFF_DIR` → Where snpEff is installed
- `PIPELINE_BASE` → Your pipeline directory

**For HTCF, your config might look like:**
```bash
MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics_analyze"
SNPEFF_DIR="/home/mihindu/software/snpEff"
PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
```

### 3. Validate Configuration

```bash
# Source the config and run validation
source pipeline_config.sh
validate_config
```

## What Gets Removed from the Pipeline

The simplification removes **Modules 4 and 5** from the next_steps file:
- ❌ **Module 4**: Mutation visualization (visualize_mutations_python_outdir.py)
- ❌ **Module 5**: Depth visualization (visualize_depth.py)

**What stays:**
- ✅ **Module 1**: Main viral assembly pipeline
- ✅ **Module 2**: Depth file generation (data only)
- ✅ **Module 3**: Mutation parsing (data only)
- ✅ **Module 7**: Diagnostic report (optional)
- ✅ **Module 8**: Consensus generation

## Key Benefits

1. **No more broken pipelines** - New genomes won't break visualization
2. **Focus on data** - All TSV/FASTA outputs still generated
3. **Flexible** - Easy to add visualization later in R or Python as needed
4. **Portable** - Config file makes it easy to move between systems

## Next Steps for Claude Code

Once you've reviewed these files and are happy with them, you can tell Claude Code to:

1. Use pipeline_config.template.sh to create the configuration system
2. Modify run_pipeline_htcf_enhanced.sh to:
   - Source the config file
   - Remove Modules 4 & 5 from next_steps template
3. Modify run_pipeline_htcf_consolidated.sh to:
   - Source the config file with fallbacks
4. Test that the pipeline still generates all data outputs

## Testing the Simplified Pipeline

After Claude Code makes the changes, test with:
```bash
# Run the enhanced pipeline
bash run_pipeline_htcf_enhanced.sh \
  sample_R1.fastq.gz \
  sample_R2.fastq.gz \
  NC_001477.1 \
  4

# Check that you get:
# ✅ Depth files (_depth.txt)
# ✅ Filtered mutation files (_filtered_mutations.tsv)
# ✅ Consensus sequences (_consensus.fasta)
# ✅ Haplotype files (_realistic_haplotypes.fasta)
# ❌ NO visualization attempts that could break
```
