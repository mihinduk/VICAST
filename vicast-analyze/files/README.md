# Pipeline Simplification - Documentation Index

## 📖 Read These In Order

### 1️⃣ START HERE: [QUICK_START.md](QUICK_START.md)
**Read this first!** 3-step guide to using Claude Code, expected results, and checklists.

### 2️⃣ VISUAL REFERENCE: [VISUAL_GUIDE_DELETIONS.md](VISUAL_GUIDE_DELETIONS.md)
Shows exactly what gets deleted from the next_steps template. Lines 178-193 before/after.

### 3️⃣ CLAUDE CODE INSTRUCTIONS: [CLAUDE_CODE_GUIDE.md](CLAUDE_CODE_GUIDE.md)
Complete instructions to copy/paste to Claude Code. Includes troubleshooting.

---

## 📁 Configuration Files

### [pipeline_config.template.sh](pipeline_config.template.sh)
Template configuration file with placeholder paths. Copy this to `pipeline_config.sh` and customize.

### [viral_genomics_analyze.yml](viral_genomics_analyze.yml)
Conda environment specification. Install with: `mamba env create -f viral_genomics_analyze.yml`

---

## 👥 For Collaborators

### [SETUP_INSTRUCTIONS.md](SETUP_INSTRUCTIONS.md)
User-facing documentation. Give this + template files to anyone setting up the pipeline.

---

## 🎯 Summary

**Goal:** Remove visualization steps from pipeline so it doesn't break with new virus genomes.

**What you need:**
- The 2 shell scripts you want to modify
- 15 minutes with Claude Code
- These documentation files

**What you get:**
- Pipeline that generates all data files
- No broken visualization code
- Easy to move between systems
- Backward compatible

**Files to send to Claude Code:**
1. `run_pipeline_htcf_enhanced.sh` (to modify)
2. `run_pipeline_htcf_consolidated.sh` (to modify)
3. Instructions from CLAUDE_CODE_GUIDE.md (copy/paste)

**Files for reference only:**
- Python scripts (don't modify)
- Configuration templates (don't modify)
- This documentation

---

## ✨ Quick Decision Tree

**"I want to use Claude Code right now"**
→ Read QUICK_START.md, then CLAUDE_CODE_GUIDE.md

**"I want to see exactly what gets deleted"**
→ Read VISUAL_GUIDE_DELETIONS.md

**"I need to give this to a collaborator"**
→ Give them: SETUP_INSTRUCTIONS.md + pipeline_config.template.sh + viral_genomics_analyze.yml

**"I want to understand the whole picture"**
→ Read all files in the numbered order above

---

## 📊 What Changes

| Component | Before | After |
|-----------|--------|-------|
| **Paths** | Hardcoded | Variables from config file |
| **Visualization** | Auto-generated | Removed (manual only) |
| **Data generation** | All modules | All modules (unchanged) |
| **Portability** | Difficult | Easy (just edit config) |
| **Breaking on new genomes** | Yes | No |
| **Backward compatible** | N/A | Yes (fallback to defaults) |

---

## 🎓 Learning Points

This simplification demonstrates:
1. **Separation of concerns**: Visualization vs. data generation
2. **Configuration management**: External config files
3. **Backward compatibility**: Fallback defaults
4. **Failure isolation**: Removing failure points without losing functionality
5. **Portability**: Easy system migration

---

Last updated: 2025-11-14
Pipeline: viral-genomics-pipeline
Environment: viral_genomics_analyze
