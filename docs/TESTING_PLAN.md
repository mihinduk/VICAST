# VICAST-annotate Testing Plan

## Overview

This document outlines the testing strategy for all 4 annotation pathways in VICAST-annotate.

## Test Genomes

### Pathway 1: Already in SnpEff
- **NC_045512** (SARS-CoV-2) - Widely available in SnpEff
- **Expected**: Detected as available, no action needed

### Pathway 2: Well-Annotated (NCBI)
- **NC_001477** (Dengue virus 1) - Well-annotated flavivirus with polyprotein
- **NC_001802** (HIV-1) - Well-annotated retrovirus
- **Expected**: Parse successfully, skip polyproteins, create editable TSV

### Pathway 3: VADR-Enhanced
- **NC_009942** (H1N1 influenza A) - Moderately annotated, can be improved with VADR
- **NC_003977** (HCV) - Flavivirus with annotation issues
- **Expected**: VADR validation improves annotations, generates alerts

### Pathway 4: BLASTx-Based
- Custom/novel viral genome with poor/no annotations
- Synthetic test case with minimal GenBank features
- **Expected**: BLASTx finds homologs, assigns genes based on similarity

## Test Cases

### Test 1: Pathway Detection (step0)

**Objective**: Verify automatic pathway detection works correctly

**Commands**:
```bash
cd vicast-annotate/

# Test pathway 1 detection
python3 step0_check_snpeff.py NC_045512

# Test pathway 2 detection  
python3 step0_check_snpeff.py NC_001477

# Test pathway 3 detection
python3 step0_check_snpeff.py NC_009942

# Test pathway 4 detection (use poorly annotated genome)
python3 step0_check_snpeff.py NC_XXXXXX
```

**Expected Output**:
- Pathway number returned as exit code
- Correct recommendation for each genome
- Annotation quality assessment displayed

**Success Criteria**:
- [ ] Correctly identifies genomes already in SnpEff
- [ ] Accurately assesses NCBI annotation quality
- [ ] Provides appropriate pathway recommendation
- [ ] Exit code matches recommended pathway (1-4)

---

### Test 2: Pathway 2 - Standard Annotation

**Objective**: Test standard annotation pipeline with well-annotated genome

**Test Genome**: NC_001477 (Dengue virus 1)

**Commands**:
```bash
# Step 1: Parse GenBank
python3 step1_parse_viral_genome.py NC_001477

# Verify outputs
ls -lh NC_001477*

# Check TSV contents
head -20 NC_001477_no_polyprotein.tsv

# Step 2: Add to SnpEff (after manual review)
python3 step2_add_to_snpeff.py NC_001477 NC_001477_no_polyprotein.tsv

# Verify database
snpeff dump NC_001477 | head -20
```

**Expected Output**:
- GenBank file downloaded
- Polyprotein features skipped
- TSV file created with individual genes
- GFF3 file generated
- SnpEff database built successfully

**Success Criteria**:
- [ ] Polyprotein correctly identified and skipped
- [ ] Individual mature peptides preserved
- [ ] Gene names extracted correctly
- [ ] TSV editable in spreadsheet software
- [ ] SnpEff database builds without errors
- [ ] Test variant annotation works

---

### Test 3: Pathway 3 - VADR-Enhanced Annotation

**Objective**: Test VADR integration and validation

**Test Genome**: NC_009942 (H1N1 influenza A)

**Pre-requisites**:
- VADR installed: `bash setup/install_vadr.sh /path/to/software`
- VADR models available
- `$VADR_DIR` environment variable set

**Commands**:
```bash
# Step 1: Parse with VADR
python3 step1_parse_viral_genome.py NC_009942 --use-vadr --vadr-model flua

# Check VADR results
ls -lh NC_009942_vadr/
cat NC_009942_vadr/NC_009942_vadr.vadr.alt

# Review curated TSV
head -20 NC_009942_vadr_curated.tsv

# Step 2: Add to SnpEff
python3 step2_add_to_snpeff.py NC_009942 NC_009942_vadr_curated.tsv
```

**Expected Output**:
- VADR annotation completes
- Alert file generated with validation issues
- Curated TSV incorporates VADR improvements
- SnpEff database built

**Success Criteria**:
- [ ] VADR runs successfully
- [ ] Alerts file generated and readable
- [ ] Annotation improvements incorporated
- [ ] VADR warnings addressed in TSV
- [ ] SnpEff database validates

---

### Test 4: Pathway 4 - BLASTx Annotation

**Objective**: Test homology-based annotation for poorly annotated genomes

**Test Genome**: Custom viral genome or NC_XXXXXX (poorly annotated)

**Pre-requisites**:
- BLAST+ installed: `conda install -c bioconda blast`
- BLAST database available (nr, refseq_protein, or custom viral DB)

**Commands**:
```bash
# Step 1: BLASTx annotation
python3 step1_blastx_annotate.py test_genome.fasta --blast-db nr --evalue 1e-5

# Check BLAST results
cat test_genome_blastx.txt | head -20

# Review TSV
head -20 test_genome_blastx.tsv

# Step 2: Add to SnpEff
python3 step2_add_to_snpeff.py test_genome test_genome_blastx.tsv
```

**Expected Output**:
- BLASTx completes successfully
- Tabular BLAST results generated
- Hits merged and deduplicated
- TSV created with homology-based annotations
- SnpEff database built

**Success Criteria**:
- [ ] BLASTx finds relevant hits
- [ ] Overlapping hits merged appropriately
- [ ] Gene names extracted from BLAST titles
- [ ] Product descriptions informative
- [ ] E-values and identity scores reported
- [ ] Manual curation checkpoint works

---

### Test 5: Master Controller

**Objective**: Test automatic pathway selection and execution

**Commands**:
```bash
# Automatic pathway selection
python3 vicast_annotate.py NC_001477

# Force specific pathway
python3 vicast_annotate.py NC_009942 --pathway 3

# Auto-proceed (skip manual curation)
python3 vicast_annotate.py NC_001477 --auto
```

**Expected Output**:
- Correct pathway detected and executed
- Manual curation checkpoint presented (unless --auto)
- Complete end-to-end execution
- Genome added to SnpEff successfully

**Success Criteria**:
- [ ] Pathway detection works
- [ ] Forced pathway respected
- [ ] Manual curation checkpoint functions
- [ ] --auto flag bypasses checkpoint
- [ ] Error handling appropriate
- [ ] Final SnpEff database verified

---

### Test 6: Validation Functions

**Objective**: Test comprehensive validation suite

**Commands**:
```bash
cd vicast-annotate/

# Validate complete setup
python3 -c "from vicast_validation import validate_vicast_setup; validate_vicast_setup()"

# Test individual components
python3 -c "from vicast_validation import validate_snpeff; validate_snpeff()"
python3 -c "from vicast_validation import validate_vadr; validate_vadr()"
python3 -c "from vicast_validation import validate_python_packages; validate_python_packages()"

# Test annotation pipeline
python3 -c "from vicast_validation import test_vicast_annotate; test_vicast_annotate('NC_001477')"
```

**Expected Output**:
- All validation checks pass
- Clear error messages for missing components
- Test annotation completes successfully

**Success Criteria**:
- [ ] Setup validation comprehensive
- [ ] Individual component checks work
- [ ] Test pipeline executes
- [ ] Helpful error messages
- [ ] Installation verification accurate

---

## Integration Testing

### Full Workflow Test

**Objective**: Test complete workflow from detection to variant annotation

**Steps**:
1. Run pathway detection
2. Execute appropriate pathway
3. Perform manual curation
4. Add to SnpEff
5. Annotate test variants
6. Verify output

**Commands**:
```bash
# 1. Detect pathway
python3 step0_check_snpeff.py NC_001477

# 2-4. Run full pipeline
python3 vicast_annotate.py NC_001477

# 5. Create test VCF
echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nNC_001477\t100\t.\tA\tG\t.\t.\t." > test.vcf

# 6. Annotate variants
snpeff NC_001477 test.vcf > annotated.vcf

# Verify annotations
grep "ANN=" annotated.vcf
```

**Success Criteria**:
- [ ] Complete workflow executes without errors
- [ ] Variant annotation includes gene names
- [ ] Effect prediction accurate
- [ ] Output VCF properly formatted

---

## Edge Cases and Error Handling

### Test 7: Error Conditions

**Test Cases**:

1. **Missing GenBank file**
   ```bash
   python3 step1_parse_viral_genome.py NC_INVALID
   # Expected: Error message, instructions to download
   ```

2. **Invalid VADR model**
   ```bash
   python3 step1_parse_viral_genome.py NC_001477 --use-vadr --vadr-model invalid
   # Expected: Error message, list available models
   ```

3. **BLAST database not found**
   ```bash
   python3 step1_blastx_annotate.py test.fasta --blast-db /invalid/path
   # Expected: Error message, installation instructions
   ```

4. **Empty TSV (all features deleted)**
   ```bash
   # Create empty TSV, run step2
   python3 step2_add_to_snpeff.py NC_001477 empty.tsv
   # Expected: Warning, prompt to continue
   ```

5. **SnpEff environment not loaded**
   ```bash
   unset SNPEFF_JAR SNPEFF_DATA
   python3 step2_add_to_snpeff.py NC_001477 test.tsv
   # Expected: Error message, instructions to source environment
   ```

**Success Criteria**:
- [ ] Graceful error handling
- [ ] Helpful error messages
- [ ] Instructions for resolution
- [ ] No crashes or stack traces shown to user

---

## Performance Testing

### Test 8: Processing Time

**Genomes**:
- Small genome (<5kb): NC_001477 (Dengue, 10.7kb)
- Medium genome (5-15kb): NC_045512 (SARS-CoV-2, 29.9kb)
- Large genome (>15kb): NC_001802 (HIV-1, 9.2kb)

**Metrics**:
- Step1 execution time
- VADR execution time (pathway 3)
- BLASTx execution time (pathway 4)
- Step2 execution time
- Total end-to-end time

**Targets**:
- Small genome: <2 minutes (pathway 2)
- Medium genome: <5 minutes (pathway 2)
- VADR annotation: <10 minutes
- BLASTx (with viral DB): <15 minutes

---

## Regression Testing

After any code changes, run:

```bash
# Quick regression test
python3 vicast_annotate.py NC_001477 --pathway 2 --auto

# Verify database
snpeff dump NC_001477 | head

# Full test suite
bash tests/run_all_tests.sh
```

---

## Test Documentation

### For Each Test:
- [ ] Test ID and name
- [ ] Objective clearly stated
- [ ] Pre-requisites listed
- [ ] Commands documented
- [ ] Expected output described
- [ ] Success criteria defined
- [ ] Actual results recorded
- [ ] Issues noted

### Test Log Template:
```
Test: [Test Name]
Date: [YYYY-MM-DD]
Tester: [Name]
Version: v2.1.0
Result: [PASS/FAIL]
Notes: [Any observations]
Issues: [Link to GitHub issues if any]
```

---

## Continuous Integration

Recommended CI tests:
1. Syntax checking (python linting)
2. Pathway 2 with test genome
3. Validation functions
4. Mock pathway detection
5. Documentation build

---

## Conclusion

This testing plan ensures all 4 pathways are thoroughly validated before release and after any updates. Each pathway should be tested independently and as part of the integrated workflow.
