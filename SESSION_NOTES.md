# VICAST Testing Session Notes

## Date: 2025-10-20

## Testing Status Summary

### Completed Testing
- ✅ **Pathway 1**: Already in SnpEff - WORKING
- ✅ **Pathway 2**: Well-annotated genome (Dengue NC_001477) - **COMPLETE SUCCESS**
- ⚠️ **Pathway 3**: VADR-enhanced annotation - **NEEDS RE-EVALUATION**
- ❌ **Pathway 4**: BLASTx annotation - NOT YET TESTED

## Pathway 2 Success Details

Successfully tested complete workflow:
1. Automatic genome download from NCBI (NC_001477 - Dengue)
2. GenBank parsing with polyprotein skipping
3. TSV generation for manual curation
4. Proper GFF3 with mRNA → CDS hierarchy
5. CDS and protein FASTA generation with correct headers
6. SnpEff database build with validation
7. Database ready for variant annotation

**Location:** `/scratch/sahlab/kathie/vicast/test/pathway2_dengue`

**Key fixes implemented:**
- Automatic NCBI download integration
- Duplicate ID resolution using gene names
- mRNA feature generation for proper hierarchy
- CDS/protein FASTA with matching IDs
- Variable scope fixes

## Pathway 3 (VADR) Issues - CRITICAL DECISION NEEDED

### Problems Encountered:

1. **Memory Requirements:**
   - VADR requires >16GB RAM for cmalign step
   - Gets killed even with 16GB SLURM allocation
   - Too memory-intensive for login nodes
   - Would require dedicated compute jobs

2. **Environment Complexity:**
   - Requires separate conda environment (vicast-vadr)
   - BLAST version conflicts (2.13 vs 2.15+)
   - Complex installation with model files
   - Path configuration issues with VADR_DIR

3. **Conda Environment Creation:**
   - Gets killed on login node during `conda env create`
   - Requires SLURM job with substantial memory
   - Long solve times for dependency resolution

4. **Practical Limitations:**
   - Cannot run interactively on HTCF login nodes
   - Adds significant complexity to pipeline
   - May not be worth the overhead for typical use cases

### VADR Fallback Works Correctly
- Script correctly falls back to standard parsing when VADR fails
- TSV still generated for manual curation
- No data loss when VADR unavailable

## ✅ DECISION MADE: VADR Pathway Removed (2025-10-20)

### Final Decision: Remove VADR (Option B - IMPLEMENTED)

**Rationale:**
- Too memory-intensive for typical HPC usage (>16GB required)
- Adds significant complexity with separate conda environment
- Most viral genomes are adequately annotated in NCBI
- Standard parsing with manual curation TSV is sufficient
- Fallback mechanism already works correctly

**Benefits Realized:**
- ✅ Single conda environment (simpler for users)
- ✅ Faster installation and testing
- ✅ Clearer documentation
- ✅ More maintainable codebase
- ✅ No functionality lost (manual curation catches issues)

### NEW Pathway Structure (Implemented)

1. **Pathway 1**: Already in SnpEff → Use existing database
2. **Pathway 2**: Well-annotated (NCBI) → Standard pipeline ✅ TESTED
3. **Pathway 3**: BLASTx annotation → Homology-based annotation
4. **Pathway 4**: Segmented viruses → Multi-chromosome genomes

### Files Updated (All Complete)

**Code Changes:**
- ✅ `vicast-annotate/README.md` - Removed Pathway 3 VADR, renumbered
- ✅ `vicast-annotate/step0_check_snpeff.py` - Updated pathway detection logic
- ✅ `vicast-annotate/step1_parse_viral_genome.py` - Removed --use-vadr flag and all VADR code
- ✅ `README.md` - Removed VADR references, simplified environment setup

**Environment Changes:**
- ✅ Deleted `environment_vadr.yml`
- ✅ Deleted `setup/install_vadr.sh`
- ✅ Deleted `setup/create_vadr_env.slurm`

**Documentation Changes:**
- ✅ `docs/ENVIRONMENT_SETUP.md` - Completely rewritten for single environment
- ✅ `docs/TESTING_GUIDE.md` - Already clean, no changes needed
- ✅ Updated validation sections to remove `validate_vadr`

## Environment Setup on HTCF

### Current Working Environment
**Location:** `/ref/sahlab/software/envs/vadr_env`

**Contains:**
- ✅ VADR 1.6.4 (if we keep it)
- ✅ pandas, numpy, biopython
- ✅ Python 3.x
- ✅ All necessary Python packages

**Conda activation:**
```bash
source /ref/sahlab/software/anaconda3/bin/activate
conda activate /ref/sahlab/software/envs/vadr_env
```

### SnpEff Configuration
**Location:** `/ref/sahlab/software/snpEff/`
- Working correctly
- Database building successful
- Validation passing

## Testing Locations on HTCF

```
/scratch/sahlab/kathie/vicast/
├── VICAST/                           # Git repository
├── test/
│   ├── pathway1_check/              # Not created yet
│   ├── pathway2_dengue/             # ✅ COMPLETE
│   │   ├── NC_001477.gb
│   │   ├── NC_001477.fasta
│   │   ├── NC_001477_no_polyprotein.tsv
│   │   ├── NC_001477_final.gff3
│   │   └── SnpEff database built
│   ├── pathway3_zika/               # ⚠️ VADR failed - memory
│   │   ├── NC_012532.gb
│   │   ├── NC_012532.fasta
│   │   └── NC_012532_vadr_curated.tsv
│   └── pathway4_blastx/             # Not started
```

## Key Technical Achievements

1. **Automatic NCBI Downloads**
   - Integrated Entrez.efetch into step1
   - Downloads both GenBank and FASTA
   - No manual download needed

2. **Proper GFF3 Structure**
   - Creates mRNA → CDS hierarchy
   - Prevents SnpEff from auto-generating TRANSCRIPT_ IDs
   - Clean gene names without superfluous counters

3. **CDS/Protein FASTA Generation**
   - Matches mRNA IDs from GFF3
   - Translates CDS to proteins
   - Enables SnpEff validation

4. **Overlapping Gene Support**
   - Correctly handles viral polyprotein cleavage products
   - Preserves overlapping features (normal in viruses)
   - SnpEff accepts overlapping mRNA/CDS

## Next Steps After Restart

1. **DECISION POINT:** Keep or remove VADR Pathway?
   - Discuss with team
   - Consider user base and typical use cases
   - Evaluate cost/benefit of complexity

2. **If Removing VADR:**
   - Renumber pathways (4→3, segmented→4)
   - Remove VADR-related code and docs
   - Simplify to single conda environment
   - Update all documentation

3. **If Keeping VADR:**
   - Document memory requirements clearly
   - Provide SLURM submission script for VADR jobs
   - Set realistic expectations about when to use
   - Consider making it truly optional/advanced

4. **Complete Testing:**
   - Test Pathway 1 (check existing genome)
   - Test Pathway 3/4 (BLASTx) if keeping that pathway
   - Test segmented virus handling
   - Create comprehensive testing documentation

5. **Documentation Updates:**
   - Update README with clear pathway descriptions
   - Document memory requirements
   - Provide example workflows
   - Add troubleshooting guide

## Current Git Status

**Repository:** https://github.com/mihinduk/VICAST
**Branch:** main
**Last commit:** Fix VADR model directory path

**Recent commits:**
- Add proper mRNA features to GFF3
- Fix duplicate ID errors
- Add CDS/protein FASTA generation
- Add automatic NCBI downloads
- Add explicit pip dependency

## Known Working Components

✅ Step 0: Pathway detection
✅ Step 1: GenBank parsing, polyprotein skipping
✅ Step 1: NCBI automatic downloads
✅ Step 1: TSV generation for manual curation
✅ Step 2: GFF3 conversion with mRNA hierarchy
✅ Step 2: CDS/protein FASTA generation
✅ Step 2: SnpEff database building
✅ Step 2: Database validation

## Issues Resolved During Session

1. ✅ Duplicate IDs in GFF3 → Use gene names directly
2. ✅ SnpEff TRANSCRIPT_ prefix mismatch → Add mRNA features
3. ✅ Missing CDS/protein FASTA → Auto-generation integrated
4. ✅ Variable scope error → Fixed snpeff_home initialization
5. ✅ Conda pip warning → Added explicit pip dependency
6. ✅ GenBank file missing → Added auto-download

## Outstanding Issues

1. ⚠️ VADR memory requirements (decision needed)
2. ❌ Pathway 4 (BLASTx) not implemented
3. ❌ Segmented virus pipeline not tested
4. ❌ Need comprehensive testing documentation

---

## Implementation Summary (2025-10-20)

### Changes Successfully Applied

All VADR references have been removed from VICAST, resulting in a streamlined 4-pathway system:

**Code Simplifications:**
- Removed ~100 lines of VADR-specific code from `step1_parse_viral_genome.py`
- Simplified pathway detection logic in `step0_check_snpeff.py`
- Updated all help text and documentation strings

**User-Facing Improvements:**
- Single conda environment (`vicast`) instead of two
- Simpler installation instructions
- Clearer pathway decision tree
- Reduced disk space requirements (~5-7 GB vs ~12 GB)

**Documentation Clarity:**
- Rewrote ENVIRONMENT_SETUP.md for single environment
- Updated all README files with new pathway structure
- Removed confusing two-environment workflow

### Ready for Next Steps

With VADR removed, the codebase is now:
- ✅ Simpler and more maintainable
- ✅ Easier for users to understand and install
- ✅ Ready for Pathway 3 (BLASTx) implementation
- ✅ Ready for testing with real-world genomes

**Recommended next actions:**
1. Test complete Pathway 2 workflow on HTCF
2. Implement Pathway 3 (BLASTx annotation) if needed
3. Test Pathway 4 (segmented viruses) thoroughly
4. Create comprehensive testing documentation
5. Consider creating example workflows

### Lessons Learned

**Design Principle Validated:** "Build SYSTEMATIC solutions, not specific fixes"
- Rather than patching VADR's memory issues, we evaluated the root need
- Removing VADR simplified the entire system without losing functionality
- Manual curation TSV provides the validation users actually need

**Complexity Budget:**
- Every optional feature has a cost in documentation, testing, and support
- VADR's cost outweighed its benefit for typical use cases
- Sometimes less is more
