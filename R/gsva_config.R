# =============================================================================
# gsva_config.R
# HIForensics — shared GSVA scoring configuration
#
# SOURCE THIS FILE in every scoring script, after PROJECT is defined.
# Every parameter that affects how expression data is converted to GSVA
# scores lives here. If you change anything here it propagates to both
# GTEx scoring (02_score_hif.R) and GEO scoring (02b_score_geo_validation.R).
#
# Usage:
#   PROJECT <- "/media/ross/New Volume1/HIForensics"
#   source(file.path(PROJECT, "R", "gsva_config.R"))
# =============================================================================

# --- Gene sets ---------------------------------------------------------------
# Single source of truth: the RDS written by scripts/01_download_msigdb.R
GENESETS_RDS <- file.path(PROJECT, "data", "msigdb", "hif_genesets.rds")

# --- GSVA parameters ---------------------------------------------------------
# These must be identical between GTEx and GEO scoring runs.
GSVA_PARAMS <- list(
    min_size            = 5L,    # minimum genes from set found in expression matrix
    max_size            = 500L,  # maximum genes (caps very large sets)
    method              = "gsva" # gsvaParam default; change to "ssgsea" etc. here only
)

# --- Expression pre-filtering ------------------------------------------------
# Applied before GSVA to reduce matrix size and noise.
# Must be identical between GTEx and GEO runs.
EXPR_FILTER <- list(
    min_tpm             = 0.1,   # TPM threshold (original scale, pre-log)
    min_expressed_frac  = 0.10   # fraction of samples that must exceed min_tpm
)

# --- Reproducibility log -----------------------------------------------------
GSVA_VERSION  <- as.character(packageVersion("GSVA"))
CONFIG_LOADED <- format(Sys.time())

cat(sprintf(
    "  [gsva_config] GSVA v%s | minSize=%d maxSize=%d | min_tpm=%.2f min_frac=%.0f%% | gene sets: %s\n",
    GSVA_VERSION,
    GSVA_PARAMS$min_size,
    GSVA_PARAMS$max_size,
    EXPR_FILTER$min_tpm,
    EXPR_FILTER$min_expressed_frac * 100,
    GENESETS_RDS
))

if (!file.exists(GENESETS_RDS)) {
    stop("Gene sets RDS not found: ", GENESETS_RDS,
         "\nRun scripts/01_download_msigdb.R first.")
}
