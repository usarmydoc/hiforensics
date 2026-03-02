#!/usr/bin/env bash
# =============================================================================
# 02_download_geo_validation.sh
# HIForensics Phase 3 — Download independent postmortem GEO datasets
#
# These are used to validate the PMI prediction model trained on GTEx V8.
# All are open access (no account required).
#
# Datasets:
#   GSE53987  — postmortem human brain, well-characterised PMI (microarray — kept for comparison)
#   GSE216281 — postmortem human frontal cortex, Parkinson's + controls, NovaSeq RNA-seq (PRIMARY)
#
# Excluded:
#   GSE45642  — PMI not in GEO series matrix (in Stanley Brain Bank supplementary table only)
#   GSE116754 — GPL13534 = Illumina 450K methylation array, hESC source (wrong dataset)
#
# NOTE: Run this during Phase 3 (validation), not Phase 1.
# Requires: GEOquery R package
# =============================================================================

set -euo pipefail

PROJECT="/media/ross/New Volume1/HIForensics"
GEO_DIR="${PROJECT}/data/validation"
LOG="${GEO_DIR}/download.log"

mkdir -p "${GEO_DIR}"
exec > >(tee -a "${LOG}") 2>&1
echo "=== GEO validation download started: $(date) ==="

# --- Datasets ----------------------------------------------------------------
# Format: accession | tissue | PMI characterisation | notes
#
# GSE45642  | brain         | multiple PMI timepoints    | Blalock et al.
# GSE53987  | brain         | well-characterised PMI     | Lanz et al.
# GSE116754 | multi-tissue  | RNA degradation timing     | designed for this

ACCESSIONS=("GSE53987" "GSE216281")

# --- Check for GEOquery ------------------------------------------------------
Rscript -e "
if (!requireNamespace('GEOquery', quietly=TRUE)) {
    stop('GEOquery not installed. Run: BiocManager::install(\"GEOquery\")')
}
cat('GEOquery available\n')
"

# --- Download via GEOquery ---------------------------------------------------
for acc in "${ACCESSIONS[@]}"; do
    out_dir="${GEO_DIR}/${acc}"
    if [[ -d "${out_dir}" ]]; then
        echo "[SKIP] ${acc} already downloaded"
        continue
    fi

    echo "[DOWN] ${acc}"
    mkdir -p "${out_dir}"

    Rscript -e "
library(GEOquery)
gse <- getGEO('${acc}', destdir='${out_dir}', GSEMatrix=TRUE, AnnotGPL=FALSE)
saveRDS(gse, '${out_dir}/${acc}_gse.rds')
cat('Saved: ${out_dir}/${acc}_gse.rds\n')
" || { echo "[FAIL] ${acc}"; continue; }

    echo "[OK]   ${acc}"
done

echo ""
echo "=== Downloaded validation datasets ==="
ls -lh "${GEO_DIR}"

echo ""
echo "=== Done: $(date) ==="
echo "Next step: python/04_validate.py"
