# Complete Instructions for Claude Code

## Files You Need to Upload to Claude Code

### Configuration Files (reference/templates):
1. `pipeline_config.template.sh` - Template configuration file
2. `viral_genomics_analyze.yml` - Conda environment spec
3. `SETUP_INSTRUCTIONS.md` - User documentation

### Scripts to Modify:
1. `run_pipeline_htcf_enhanced.sh` - Main wrapper script
2. `run_pipeline_htcf_consolidated.sh` - Consolidated pipeline script

### Reference (don't modify, just for context):
1. `generate_realistic_haplotype_consensus.py`
2. `parse_snpeff_tsv.py`
3. `generate_filtered_consensus.py`
4. `visualize_mutations_python_outdir.py` (to understand what we're removing)
5. `visualize_depth.py` (to understand what we're removing)

---

## What to Tell Claude Code

Copy and paste this to Claude Code:

```
I need to simplify my viral genomics pipeline by removing visualization steps and parameterizing hardcoded paths. Here are the changes needed:

**CONTEXT:**
I have a viral genomics pipeline that currently breaks when adding new virus genomes because the visualization code has hardcoded assumptions. I want to strip out the visualization and focus on data generation only. I've created config files and now need you to modify the shell scripts.

**FILES PROVIDED:**
- pipeline_config.template.sh (reference - DO NOT MODIFY)
- run_pipeline_htcf_enhanced.sh (MODIFY THIS)
- run_pipeline_htcf_consolidated.sh (MODIFY THIS)
- Other Python scripts (reference only, don't touch)

**CHANGES NEEDED:**

**1. Modify `run_pipeline_htcf_enhanced.sh`:**
   
   a) At the top of the script (around line 30-40, after the argument parsing), add:
   ```bash
   # Source configuration file if it exists
   CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"
   if [ -f "$CONFIG_FILE" ]; then
       echo "Loading configuration from: $CONFIG_FILE"
       source "$CONFIG_FILE"
   else
       echo "⚠️  Configuration file not found: $CONFIG_FILE"
       echo "   Using hardcoded paths (consider creating pipeline_config.sh from template)"
       # Fallback to hardcoded paths for backward compatibility
       MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics_analyze"
       PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
   fi
   ```
   
   b) Replace ALL instances of the hardcoded path:
   - `/home/mihindu/miniforge3/bin/mamba run -n viral_genomics` 
   with:
   - `$MAMBA_CMD`
   
   c) Replace ALL instances of:
   - `/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline`
   with:
   - `${PIPELINE_BASE}`
   
   d) Find the section that creates the "next_steps" file (starts around line 147 with `cat > "next_steps_${SAMPLE_NAME}.txt"`).
      In that section:
      - COMPLETELY REMOVE Module 4 (visualize_mutations_python_outdir.py section)
      - COMPLETELY REMOVE Module 5 (visualize_depth.py section)
      - Keep Modules 2, 3, 7, and 8
      - Add a comment where Modules 4 & 5 were:
        ```bash
        # Modules 4 & 5 (visualization) removed - focus on data generation
        # Visualization broke with new genomes due to hardcoded assumptions
        # Use generated TSV files for custom visualization as needed
        ```
   
   e) Update the usage/help message to mention the configuration file:
      Add a line: "  Configuration: Create pipeline_config.sh from pipeline_config.template.sh"

**2. Modify `run_pipeline_htcf_consolidated.sh`:**
   
   a) At the top of the script (around line 30, before "Set up environment"), add:
   ```bash
   # Source configuration file if it exists, with fallbacks
   SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
   CONFIG_FILE="${SCRIPT_DIR}/pipeline_config.sh"
   
   if [ -f "$CONFIG_FILE" ]; then
       echo "Loading configuration from: $CONFIG_FILE"
       source "$CONFIG_FILE"
   else
       echo "⚠️  No configuration file found. Using hardcoded paths."
       # Fallback defaults for backward compatibility
       MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics_analyze"
       SNPEFF_DIR="/home/mihindu/software/snpEff"
       SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
       JAVA_PATH="java"
   fi
   ```
   
   b) Remove or comment out these hardcoded lines (they're now in config):
   - `MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics"`
   - `SNPEFF_JAR="/home/mihindu/software/snpEff/snpEff.jar"`
   - `SNPEFF_DIR="/home/mihindu/software/snpEff"`
   - `JAVA_PATH="java"`
   
   c) Update the usage/help message to mention the configuration file.

**3. Create the actual `pipeline_config.sh` file:**
   - Copy pipeline_config.template.sh to pipeline_config.sh
   - Fill in the HTCF-specific paths:
     ```bash
     MAMBA_CMD="/home/mihindu/miniforge3/bin/mamba run -n viral_genomics_analyze"
     SNPEFF_DIR="/home/mihindu/software/snpEff"
     SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
     JAVA_PATH="java"
     PIPELINE_BASE="/scratch/sahlab/kathie/viral_genomics_pipeline_dev/viral-genomics-pipeline"
     CONSOLIDATED_PIPELINE="${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh"
     VIRAL_PIPELINE_SCRIPT="${PIPELINE_BASE}/viral_pipeline.py"
     ```

**KEY PRINCIPLES:**
- Don't touch any Python scripts
- Keep all data generation steps intact
- Ensure backward compatibility (scripts work even without config file)
- Remove ONLY the visualization command generation, not any data processing
- Make sure the changes are consistent across both shell scripts

**VERIFICATION:**
After making changes, please:
1. Show me the sections you modified
2. Confirm that Modules 2, 3, 7, and 8 are still in the next_steps file
3. Confirm that Modules 4 and 5 are completely removed
4. Verify all instances of hardcoded paths are replaced with variables
```

---

## How to Work with Claude Code

### Step 1: Start Claude Code session
```bash
cd /path/to/your/viral-genomics-pipeline
```

### Step 2: Upload files to Claude Code
When Claude Code asks what files you want to work with, upload:
- `run_pipeline_htcf_enhanced.sh`
- `run_pipeline_htcf_consolidated.sh`
- `pipeline_config.template.sh` (for reference)

### Step 3: Paste the instructions
Copy the instructions from the box above and paste them into Claude Code.

### Step 4: Review the changes
Claude Code will show you the modifications. Review them carefully and ask it to:
- Show you the exact lines it removed from the next_steps template
- Confirm Module 4 and Module 5 are gone
- Show you that Modules 2, 3, 7, and 8 remain

### Step 5: Test
After Claude Code makes the changes, test with:
```bash
# Make sure scripts are executable
chmod +x run_pipeline_htcf_enhanced.sh
chmod +x run_pipeline_htcf_consolidated.sh

# Check the configuration loads
bash run_pipeline_htcf_enhanced.sh --help
```

---

## Expected Results

After Claude Code finishes, you should have:

### Modified Files:
1. `run_pipeline_htcf_enhanced.sh` - Sources config, removed visualization modules
2. `run_pipeline_htcf_consolidated.sh` - Sources config with fallbacks
3. `pipeline_config.sh` - New file with your HTCF paths

### Unchanged Files:
- All Python scripts remain untouched
- All data generation logic intact

### What Works:
✅ Pipeline runs and generates all TSV/FASTA outputs
✅ Works on new virus genomes without breaking
✅ Portable (easy to move to different systems via config file)
✅ Backward compatible (runs with or without config file)

### What's Removed:
❌ Module 4: Automatic mutation visualization PNG generation
❌ Module 5: Automatic depth visualization PNG/HTML generation

You can still manually run those scripts later if needed, but they won't be part of the automated pipeline.

---

## Troubleshooting

If Claude Code gets confused:
1. **Ask it to show you specific sections**: "Show me the next_steps template section where you removed Modules 4 and 5"
2. **Be specific about line numbers**: "The next_steps section starts around line 147"
3. **Ask for verification**: "List all the MAMBA_CMD replacements you made"
4. **One script at a time**: Start with run_pipeline_htcf_enhanced.sh, then do consolidated

If you need to redo:
```bash
# Restore from backup or git
git checkout run_pipeline_htcf_enhanced.sh
git checkout run_pipeline_htcf_consolidated.sh
# Then try again with more specific instructions
```
