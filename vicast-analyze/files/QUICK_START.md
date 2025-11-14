# QUICK START: Pipeline Simplification with Claude Code

## 📦 Files You Have

### 1. Configuration & Documentation (Already Created ✅)
- **pipeline_config.template.sh** - Template with placeholder paths
- **viral_genomics_analyze.yml** - Conda environment definition
- **SETUP_INSTRUCTIONS.md** - Full setup guide for users
- **CLAUDE_CODE_GUIDE.md** - Complete instructions for Claude Code
- **VISUAL_GUIDE_DELETIONS.md** - Shows exactly what to delete

### 2. Scripts to Upload to Claude Code
- **run_pipeline_htcf_enhanced.sh** - Needs modification
- **run_pipeline_htcf_consolidated.sh** - Needs modification

### 3. Reference Scripts (Optional)
- **generate_realistic_haplotype_consensus.py** - Context only
- **parse_snpeff_tsv.py** - Context only
- **visualize_mutations_python_outdir.py** - Shows what we're removing
- **visualize_depth.py** - Shows what we're removing

---

## 🚀 Three-Step Process

### STEP 1: Read the Visual Guide (2 minutes)
Open **VISUAL_GUIDE_DELETIONS.md** to see exactly what gets deleted from the next_steps template.

### STEP 2: Use Claude Code (10-15 minutes)
1. Open Claude Code in your terminal
2. Navigate to your pipeline directory
3. Upload these two scripts:
   - `run_pipeline_htcf_enhanced.sh`
   - `run_pipeline_htcf_consolidated.sh`
4. Copy/paste the instructions from **CLAUDE_CODE_GUIDE.md**
5. Review Claude Code's changes
6. Accept and save

### STEP 3: Test (5 minutes)
```bash
# Verify the scripts load properly
bash run_pipeline_htcf_enhanced.sh --help

# Run a test
bash run_pipeline_htcf_enhanced.sh \
  test_R1.fastq.gz \
  test_R2.fastq.gz \
  NC_001477.1 \
  4
```

---

## 🎯 What Gets Changed

### In `run_pipeline_htcf_enhanced.sh`:
1. Add config file sourcing at the top
2. Replace hardcoded paths with `${MAMBA_CMD}` and `${PIPELINE_BASE}`
3. **DELETE** Module 4 (mutation visualization) from next_steps
4. **DELETE** Module 5 (depth visualization) from next_steps
5. **KEEP** Modules 2, 3, 7, 8 (all data generation)

### In `run_pipeline_htcf_consolidated.sh`:
1. Add config file sourcing at the top
2. Replace hardcoded `MAMBA_CMD`, `SNPEFF_DIR`, `SNPEFF_JAR`
3. Add fallback defaults for backward compatibility

### New File Created:
- `pipeline_config.sh` - Your actual configuration (from template)

---

## ✅ Expected Results

### After Claude Code finishes:
- ✅ Pipeline generates all TSV and FASTA files
- ✅ Works with new virus genomes without breaking
- ✅ No automatic visualization attempts
- ✅ Backward compatible (runs even without config file)
- ✅ Easy to move between systems

### What You'll Still Get:
- ✅ Depth files (`_depth.txt`)
- ✅ Filtered mutations (`_filtered_mutations.tsv`)
- ✅ Consensus sequences (`_consensus.fasta`)
- ✅ Haplotype files (`_realistic_haplotypes.fasta`)
- ✅ Protein sequences (`_realistic_proteins.fasta`)

### What Won't Auto-Generate:
- ❌ Mutation visualization PNGs (can run manually later if needed)
- ❌ Depth visualization PNGs/HTML (can run manually later if needed)

---

## 💡 Pro Tips

### If Claude Code Gets Confused:
1. Show it **VISUAL_GUIDE_DELETIONS.md** - point to specific line numbers
2. Ask it to work on one script at a time
3. Ask it to show you what it changed before accepting

### After Changes:
1. Commit to git: `git add -A && git commit -m "Simplified pipeline - removed visualization"`
2. Create pipeline_config.sh: `cp pipeline_config.template.sh pipeline_config.sh`
3. Edit paths in pipeline_config.sh for your system
4. Test with a small sample

### For Collaborators:
1. Give them: `pipeline_config.template.sh` + `SETUP_INSTRUCTIONS.md`
2. They copy template → config
3. They edit paths for their system
4. They create the conda environment from the YAML file

---

## 📋 Checklist

Before starting with Claude Code:
- [ ] Read VISUAL_GUIDE_DELETIONS.md
- [ ] Have run_pipeline_htcf_enhanced.sh ready
- [ ] Have run_pipeline_htcf_consolidated.sh ready
- [ ] Read CLAUDE_CODE_GUIDE.md instructions

After Claude Code finishes:
- [ ] Verify Module 4 and 5 are gone from next_steps
- [ ] Verify Modules 2, 3, 7, 8 remain
- [ ] Check all hardcoded paths are replaced
- [ ] Create pipeline_config.sh from template
- [ ] Edit pipeline_config.sh with your paths
- [ ] Test pipeline runs without errors
- [ ] Commit changes to git

---

## 🆘 Need Help?

If something doesn't work:
1. Check the config file exists: `ls -la pipeline_config.sh`
2. Source it manually: `source pipeline_config.sh`
3. Check variables are set: `echo $MAMBA_CMD`
4. Verify paths exist: `ls $SNPEFF_JAR`

The pipeline should work even without the config file (using hardcoded fallbacks).
