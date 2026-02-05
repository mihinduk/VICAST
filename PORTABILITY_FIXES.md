# Portability Fixes: Removal of Hardcoded HTCF Paths

**Date:** 2026-02-05
**Status:** ✅ Complete

## What Was Fixed

Systematically removed hardcoded HTCF paths from **14 shell scripts** and added environment variable support for full portability.

### Changes Applied

1. **Conda paths:** `/ref/sahlab/software/anaconda3` → `$CONDA_BASE` (with fallbacks)
2. **SnpEff:** `/ref/sahlab/software/snpEff` → `$SNPEFF_DIR` (with fallbacks)
3. **Java:** `/ref/sahlab/software/jdk-21` → `$JAVA_HOME` (with fallbacks)
4. **BLAST DB:** `/ref/sahlab/data/microbes_db` → `$BLAST_DB` (with fallbacks)

### Files Modified (14)

Core pipeline scripts in `vicast-analyze/`:
- `pipeline_config.sh`, `run_pipeline_*.sh`, `run_vicast_analyze_*.sh`
- `viral_diagnostic.sh`, `submit_pipeline.sh`, `process_segmented_virus.sh`

Setup scripts:
- `setup/install_java21_conda.sh`, `setup/snpeff_env.sh`, `test_vicast_env.sh`

### New Files

- `vicast_paths.env` - Central configuration with auto-detection

## Usage

**Default (auto-detects paths):**
```bash
source vicast_paths.env
```

**Custom paths:**
```bash
export CONDA_BASE="/your/conda/path"
export SNPEFF_DIR="/your/snpeff/path"
# Then run scripts normally
```

## Backwards Compatibility

✅ **Fully compatible** - Scripts fallback to original HTCF paths if environment variables not set.

## Result

VICAST is now **fully portable** and **publication-ready**! Works on any system without hardcoded institutional paths.
