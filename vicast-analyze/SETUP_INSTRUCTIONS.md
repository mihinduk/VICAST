# VICAST-Analyze Setup Instructions

## Quick Start

### Option 1: Automated Setup (Recommended)

```bash
# Clone the repository
git clone https://github.com/mihinduk/VICAST.git
cd VICAST/vicast-analyze

# Activate conda (HTCF users)
source /ref/sahlab/software/anaconda3/bin/activate

# Run the setup script
bash setup_environment.sh
```

The setup script will:
- ⚠️ Warn you about home directory space limits
- Help you choose an appropriate installation location
- Create the `viral_genomics_analyze` environment
- Tell you how to update `pipeline_config.sh`

### Option 2: Manual Setup

**⚠️ IMPORTANT: Do NOT install in your home directory on shared servers!**
Conda environments are large (>2GB) and will exceed quota limits.

#### 1. Create the Conda Environment

```bash
# Navigate to the pipeline directory
cd /path/to/VICAST/vicast-analyze

# Activate conda
source /ref/sahlab/software/anaconda3/bin/activate

# Set environment location (HTCF shared lab - RECOMMENDED)
export CONDA_ENVS_DIRS="/ref/sahlab/software/envs:$HOME/.conda/envs"

# Create the environment
conda env create -f viral_genomics_analyze.yml

# Activate it
conda activate viral_genomics_analyze

# Verify installation
bwa
samtools --version
python --version
which efetch
```

#### 2. Configure Your Pipeline Paths

```bash
# The default pipeline_config.sh uses the existing HTCF shared environment
# If you created viral_genomics_analyze, update it:
nano pipeline_config.sh
```

**Update MAMBA_CMD if you created the new environment:**
```bash
# Change from:
MAMBA_CMD="conda run -n viral_genomics"

# To:
MAMBA_CMD="conda run -n viral_genomics_analyze"
```

**Other paths to verify (usually correct by default on HTCF):**
- `SNPEFF_DIR` → `/ref/sahlab/software/snpEff`
- `PIPELINE_BASE` → Your pipeline directory path

#### 3. Run Your First Test

**Prerequisites:**
- Activate conda in your session: `source /ref/sahlab/software/anaconda3/bin/activate`

**Run the pipeline:**
```bash
cd /path/to/your/data

bash /path/to/VICAST/vicast-analyze/run_pipeline_htcf_enhanced.sh \
  sample_R1.fastq.gz \
  sample_R2.fastq.gz \
  NC_001477.1 \
  4
```

## Environment Options

### For HTCF Users

**Option A: Use existing shared environment** (Quick start)
- Environment: `viral_genomics` (already installed)
- No setup needed
- Config already set correctly

**Option B: Create dedicated environment** (Independent, systematic)
- Environment: `viral_genomics_analyze` (you create)
- Run `bash setup_environment.sh`
- Update `pipeline_config.sh` with `MAMBA_CMD="conda run -n viral_genomics_analyze"`
- Advantage: Complete independence from other pipelines

### For New Systems

Use `setup_environment.sh` to guide installation with space warnings.

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
