#!/usr/bin/env bash
# =============================================================================
# 00_download_gtex.sh
# HIForensics Phase 1 — Download GTEx V8 expression data and annotations
#
# Using V8 (not V10) — fully open access, no dbGaP required.
#
# Downloads:
#   1. Combined gene TPM matrix (all tissues, ~14 GB compressed)
#   2. Sample attributes (contains SMTSISCH = ischemic time, minutes)
#   3. Subject phenotypes (age, sex, Hardy death classification)
#
# Note: gene symbol mapping uses the Description column in the GCT file
# directly — no separate GTF required.
#
# GCS bucket: gs://adult-gtex/
# Portal:     https://gtexportal.org/home/datasets
# =============================================================================

set -euo pipefail

# --- Paths -------------------------------------------------------------------
PROJECT="/media/ross/New Volume1/HIForensics"
GTEX_DIR="${PROJECT}/data/gtex"
LOG="${GTEX_DIR}/download.log"

mkdir -p "${GTEX_DIR}"
exec > >(tee -a "${LOG}") 2>&1
echo "=== GTEx V8 download started: $(date) ==="

# --- GTEx V8 URLs (verified 2026-03-01) --------------------------------------
BASE_EXPR="https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq"
BASE_ANNO="https://storage.googleapis.com/adult-gtex/annotations/v8/metadata-files"

declare -A FILES=(
    ["GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"]="${BASE_EXPR}"
    ["GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"]="${BASE_ANNO}"
    ["GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"]="${BASE_ANNO}"
)

# --- Download function -------------------------------------------------------
download_file() {
    local filename="$1"
    local base_url="$2"
    local dest="${GTEX_DIR}/${filename}"

    if [[ -f "${dest}" ]]; then
        echo "[SKIP] ${filename} already exists ($(du -sh "${dest}" | cut -f1))"
        return 0
    fi

    echo "[DOWN] ${filename}"
    wget \
        --no-verbose \
        --show-progress \
        --continue \
        --retry-connrefused \
        --tries=5 \
        --timeout=60 \
        -O "${dest}" \
        "${base_url}/${filename}" \
        || { echo "[FAIL] ${filename}"; rm -f "${dest}"; return 1; }

    echo "[OK]   ${filename} ($(du -sh "${dest}" | cut -f1))"
}

# --- Run downloads -----------------------------------------------------------
for filename in "${!FILES[@]}"; do
    download_file "${filename}" "${FILES[$filename]}"
done

# --- Verify ------------------------------------------------------------------
echo ""
echo "=== Downloaded files ==="
ls -lh "${GTEX_DIR}"

echo ""
echo "=== Key column check: SMTSISCH in sample attributes ==="
SA_FILE="${GTEX_DIR}/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
if [[ -f "${SA_FILE}" ]]; then
    head -1 "${SA_FILE}" \
        | tr '\t' '\n' \
        | grep -n "SMTSISCH\|SMTSD\|SAMPID\|SMAFRZE\|DTHHRDY" \
        || echo "WARNING: Expected columns not found — check file format"
fi

echo ""
echo "=== Download complete: $(date) ==="
echo "Next step: Rscript R/02_score_hif.R"
