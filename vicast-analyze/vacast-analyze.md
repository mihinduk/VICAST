I need to simplify my viral genomics pipeline by removing visualization steps and parameterizing hardcoded paths. Here are the changes needed:

**1. Create a new configuration file: `pipeline_config.sh`**
This should contain all environment-specific paths as variables:
- MAMBA_CMD="/path/to/mamba"
- SNPEFF_DIR="/path/to/snpEff"
- SNPEFF_JAR="${SNPEFF_DIR}/snpEff.jar"
- PIPELINE_BASE="/path/to/viral-genomics-pipeline"
- CONSOLIDATED_PIPELINE="${PIPELINE_BASE}/run_pipeline_htcf_consolidated.sh"

**2. Modify `run_pipeline_htcf_enhanced.sh`:**
- Add `source pipeline_config.sh` at the top
- Replace all hardcoded paths with config variables
- In the "next_steps" template (starts around line 147), REMOVE Module 4 and Module 5 entirely
- Keep Modules 2, 3, 7, and 8
- Add a comment where Modules 4 & 5 were: `# Modules 4 & 5 (visualization) removed - focus on data generation`

**3. Modify `run_pipeline_htcf_consolidated.sh`:**
- Add `source pipeline_config.sh` at the top (if the file exists)
- Replace hardcoded MAMBA_CMD, SNPEFF_DIR, and SNPEFF_JAR with config variables
- Add fallback defaults if config file doesn't exist (for backward compatibility)

**4. Add a template config file: `pipeline_config.template.sh`**
Create this with placeholder values and instructions for users to copy and customize.

**5. Update the usage/help messages** in both wrapper scripts to mention the configuration file.

**Key principles:**
- Don't touch the Python scripts - they're fine as-is
- Keep all data generation intact
- Only remove the visualization module calls from the next_steps template
- Make sure the pipeline can still run if someone has the old hardcoded paths (backward compatibility)
- The consolidated pipeline should complete successfully and generate all TSV/FASTA outputs