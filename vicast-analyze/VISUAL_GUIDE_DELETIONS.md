# VISUAL GUIDE: What to Remove from next_steps Template

## Location in run_pipeline_htcf_enhanced.sh
Lines 162-211: The `cat > "next_steps_${SAMPLE_NAME}.txt" << EOF` section

---

## BEFORE (Current - Lines 178-193):

```bash
# Module 4: Visualize mutations
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \
  --output ${SAMPLE_NAME}_mutations.png \
  --accession ${ACCESSION} \
  --cutoff 0.01 \
  --outdir ${SAMPLE_NAME}_results

# Module 5: Visualize depth
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/visualize_depth.py \
  --depth ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt \
  --output ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.png \
  --output-html ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.html \
  --accession ${ACCESSION}
```

---

## AFTER (What it should become):

```bash
# Modules 4 & 5 (visualization) removed - focus on data generation
# Visualization broke with new genomes due to hardcoded assumptions
# Use generated TSV files for custom visualization as needed
```

---

## COMPLETE AFTER STRUCTURE:

The entire next_steps section should look like this after modifications:

```bash
    # Create next steps file
    cat > "next_steps_${SAMPLE_NAME}.txt" << EOF
# Next steps for ${SAMPLE_NAME} (${ACCESSION})

# Module 2: Generate depth file
${MAMBA_CMD} bash -c '
  mkdir -p ${SAMPLE_NAME}_results &&
  echo -e "chrom\\tposition\\tdepth" > ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt &&
  samtools depth cleaned_seqs/mapping/${SAMPLE_NAME}.lofreq.final.bam >> ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt'

# Module 3: Parse mutations
${MAMBA_CMD} \\
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/parse_snpeff_tsv.py \\
  "cleaned_seqs/variants/${SAMPLE_NAME}.snpEFF.ann.tsv" \\
  "${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv" \\
  --quality 1000 --depth 200 --freq 0.01

# Modules 4 & 5 (visualization) removed - focus on data generation
# Visualization broke with new genomes due to hardcoded assumptions
# Use generated TSV files for custom visualization as needed

# Module 7: Diagnostic Report (Optional)
sbatch ${PIPELINE_BASE}/viral_pipeline/analysis/submit_viral_diagnostic.sh \\
  "${R1}" \\
  "${R2}" \\
  ${ACCESSION} \\
  diagnostic_${SAMPLE_NAME} \\
  4

# Module 8: Generate consensus
${MAMBA_CMD} \\
  python3 ${PIPELINE_BASE}/viral_pipeline/utils/generate_filtered_consensus.py \\
  --vcf ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \\
  --reference cleaned_seqs/${ACCESSION}.fasta \\
  --accession ${ACCESSION} \\
  --quality 1000 --depth 200 --freq 0.50 \\
  --output-prefix ${SAMPLE_NAME}_results/${SAMPLE_NAME}_consensus
EOF
```

---

## KEY CHANGES TO VERIFY:

1. ✅ Lines 178-193 (Module 4 and Module 5) completely removed
2. ✅ Replaced with 3-line comment block
3. ✅ All `/home/mihindu/miniforge3/bin/mamba run -n viral_genomics` replaced with `${MAMBA_CMD}`
4. ✅ Module 2 kept (depth file generation)
5. ✅ Module 3 kept (mutation parsing)
6. ✅ Module 7 kept (diagnostic report)
7. ✅ Module 8 kept (consensus generation)
8. ✅ `${PIPELINE_BASE}` variable used throughout

---

## EXACT LINES TO DELETE:

Delete these exact lines (178-193):

```
# Module 4: Visualize mutations
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/visualize_mutations_python_outdir.py \
  --input ${SAMPLE_NAME}_results/${SAMPLE_NAME}_filtered_mutations.tsv \
  --output ${SAMPLE_NAME}_mutations.png \
  --accession ${ACCESSION} \
  --cutoff 0.01 \
  --outdir ${SAMPLE_NAME}_results

# Module 5: Visualize depth
/home/mihindu/miniforge3/bin/mamba run -n viral_genomics \
  python3 ${PIPELINE_BASE}/viral_pipeline/visualization/visualize_depth.py \
  --depth ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.txt \
  --output ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.png \
  --output-html ${SAMPLE_NAME}_results/${SAMPLE_NAME}_depth.html \
  --accession ${ACCESSION}
```

And replace with:

```
# Modules 4 & 5 (visualization) removed - focus on data generation
# Visualization broke with new genomes due to hardcoded assumptions
# Use generated TSV files for custom visualization as needed
```
