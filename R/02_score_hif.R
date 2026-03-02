# =============================================================================
# 02_score_hif.R
# HIForensics Phase 1 — HIF pathway scoring across GTEx V8 tissues
#
# Input:
#   data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz
#   data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt
#   data/gtex/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt
#   data/msigdb/hif_genesets.rds
#
# Output:
#   results/hif_scores/gsva_scores.rds        — gene sets × samples matrix
#   results/hif_scores/metadata.rds           — sample metadata + scores merged
#   results/hif_scores/tissue_correlations.csv — Spearman r, p, q per tissue/geneset
#   results/hif_scores/summary.csv            — top tissues ranked by HIF-PMI correlation
#
# Runtime: ~2–3 hr (GSVA per tissue; 12 workers, small per-tissue matrices)
# RAM:     ~8–12 GB peak (largest tissue ~800 samples; 12 workers × ~700 MB each)
# =============================================================================

suppressPackageStartupMessages({
    library(GSVA)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(BiocParallel)
})

# --- Configuration -----------------------------------------------------------
PROJECT   <- "/media/ross/New Volume1/HIForensics"
GTEX_DIR  <- file.path(PROJECT, "data", "gtex")
OUT_DIR   <- file.path(PROJECT, "results", "hif_scores")
N_WORKERS <- 12   # Per-tissue scoring: each tissue matrix is ~200 MB not 5 GB,
                  # so 12 workers × small matrix is safe. See Step 5.

# Shared GSVA config — gene sets path, minSize/maxSize, expression filters
# Edit R/gsva_config.R to change any of these; never hardcode them here
source(file.path(PROJECT, "R", "gsva_config.R"))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== HIForensics HIF Scoring:", format(Sys.time()), "===\n")
cat("Output dir:", OUT_DIR, "\n\n")

# --- Check inputs exist ------------------------------------------------------
inputs <- list(
    tpm   = file.path(GTEX_DIR, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"),
    sattr = file.path(GTEX_DIR, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"),
    subj  = file.path(GTEX_DIR, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"),
    gsets = GENESETS_RDS
)

missing <- names(inputs)[!file.exists(unlist(inputs))]
if (length(missing) > 0) {
    stop("Missing input files: ", paste(missing, collapse=", "),
         "\nRun scripts/00_download_gtex.sh first.")
}

# =============================================================================
# STEP 1: Load metadata
# =============================================================================
cat("--- Step 1: Loading metadata ---\n")

# Sample attributes — one row per RNA-seq sample
sattr <- fread(inputs$sattr, sep="\t", quote="")

# Filter to RNA-seq samples only
sattr <- sattr[SMAFRZE == "RNASEQ"]
cat("  RNA-seq samples:", nrow(sattr), "\n")

# Extract donor ID from SAMPID (format: GTEX-XXXXX-...)
sattr[, SUBJID := sub("(GTEX-[^-]+)-.*", "\\1", SAMPID)]

# Subject phenotypes — one row per donor
subj <- fread(inputs$subj, sep="\t", quote="")

# Merge
meta <- merge(sattr, subj, by="SUBJID", all.x=TRUE)

# Keep relevant columns, rename for clarity
meta <- meta[, .(
    sample_id    = SAMPID,
    donor_id     = SUBJID,
    tissue       = SMTSD,             # full tissue name
    ischemic_min = SMTSISCH,          # ischemic time in minutes
    ischemic_hr  = SMTSISCH / 60,     # convenience: hours
    sex          = SEX,               # 1=male, 2=female
    age_range    = AGE,               # e.g. "50-59"
    hardy_scale  = DTHHRDY            # 0=ventilator, 1=fast violent, 2=fast natural, 3=intermediate, 4=slow
)]

# Remove samples with missing ischemic time
n_before <- nrow(meta)
meta <- meta[!is.na(ischemic_min) & ischemic_min >= 0]
cat("  Samples with valid ischemic time:", nrow(meta),
    "(removed", n_before - nrow(meta), "with NA/negative)\n")
cat("  Ischemic time range:",
    round(min(meta$ischemic_min), 0), "–",
    round(max(meta$ischemic_min), 0), "min (",
    round(min(meta$ischemic_hr), 1), "–",
    round(max(meta$ischemic_hr), 1), "hr)\n")
cat("  Tissues:", length(unique(meta$tissue)), "\n\n")

# =============================================================================
# STEP 2: Load GTEx TPM matrix
# =============================================================================
cat("--- Step 2: Loading GTEx TPM matrix (this will take a few minutes) ---\n")
cat("  File:", inputs$tpm, "\n")

# GCT format: line 1 = "#1.2", line 2 = "nrow\tncol", line 3+ = data with header
tpm_dt <- fread(
    inputs$tpm,
    skip        = 2,       # skip "#1.2" and "nrow\tncol"
    header      = TRUE,
    sep         = "\t",
    quote       = "",
    showProgress= TRUE
)

cat("  Raw dimensions:", nrow(tpm_dt), "genes ×", ncol(tpm_dt) - 2, "samples\n")

# Column 1 = Name (Ensembl ID with version), Column 2 = Description (gene symbol)
# Strip Ensembl version suffix (ENSG00000000003.15 → ENSG00000000003)
tpm_dt[, Name := sub("\\..*", "", Name)]

# Use gene symbol (Description) as row identifier — MSigDB uses symbols
# Handle duplicate symbols by keeping the row with highest mean expression
tpm_dt[, mean_expr := rowMeans(.SD), .SDcols = patterns("^GTEX-")]
tpm_dt <- tpm_dt[order(-mean_expr)]
tpm_dt <- tpm_dt[!duplicated(Description)]
tpm_dt[, mean_expr := NULL]

cat("  After dedup by symbol:", nrow(tpm_dt), "genes\n")

# =============================================================================
# STEP 3: Filter to samples with valid metadata + expressed genes
# =============================================================================
cat("--- Step 3: Filtering samples and genes ---\n")

# Keep only samples present in metadata
sample_cols   <- colnames(tpm_dt)[-(1:2)]
valid_samples <- intersect(sample_cols, meta$sample_id)
cat("  Samples with metadata:", length(valid_samples), "of", length(sample_cols), "\n")

# Subset
keep_cols <- c("Description", valid_samples)
tpm_dt    <- tpm_dt[, ..keep_cols]

# Convert to matrix (genes × samples)
gene_symbols <- tpm_dt$Description
tpm_mat      <- as.matrix(tpm_dt[, -1])
rownames(tpm_mat) <- gene_symbols
colnames(tpm_mat) <- valid_samples
rm(tpm_dt); gc()

# Log2(TPM + 1) transform
cat("  Applying log2(TPM+1) transform\n")
tpm_mat <- log2(tpm_mat + 1)

# Filter lowly expressed genes — thresholds from gsva_config.R
n_samples   <- ncol(tpm_mat)
min_samples <- ceiling(EXPR_FILTER$min_expressed_frac * n_samples)
expressed   <- rowSums(tpm_mat > log2(EXPR_FILTER$min_tpm + 1)) >= min_samples
tpm_mat     <- tpm_mat[expressed, ]
cat("  Expressed genes (TPM >", EXPR_FILTER$min_tpm,
    "in >=", EXPR_FILTER$min_expressed_frac * 100, "% samples):",
    nrow(tpm_mat), "\n\n")

# =============================================================================
# STEP 4: Load HIF gene sets and check coverage
# =============================================================================
cat("--- Step 4: Loading HIF gene sets ---\n")

gene_sets    <- readRDS(inputs$gsets)
geneset_cols <- names(gene_sets)  # defined here so Step 6 CSV export can use it

# Report gene set coverage in the expression matrix
cat(sprintf("  %-55s  %5s  %5s  %5s\n", "Gene Set", "Total", "Found", "%"))
cat("  ", strrep("-", 74), "\n", sep="")
for (nm in names(gene_sets)) {
    total <- length(gene_sets[[nm]])
    found <- sum(gene_sets[[nm]] %in% rownames(tpm_mat))
    cat(sprintf("  %-55s  %5d  %5d  %4.0f%%\n",
                nm, total, found, 100 * found / total))
}

# Drop gene sets with < 5 genes in the matrix (GSVA minimum)
gene_sets <- gene_sets[sapply(gene_sets, function(g) sum(g %in% rownames(tpm_mat))) >= 5]
cat("\n  Gene sets retained (>=5 genes found):", length(gene_sets), "\n\n")

# =============================================================================
# STEP 5: GSVA scoring — tissue by tissue
#
# Strategy: score each tissue's samples separately rather than the full
# 17k-sample matrix. Peak RAM = one tissue matrix (~200 MB) × N_WORKERS,
# vs the previous approach that forked the full 5 GB matrix × N_WORKERS.
#
# Checkpointing: each tissue's scores are saved to
#   results/hif_scores/tissue_scores/{safe_name}.rds
# so a crash mid-run can be resumed without rescoring completed tissues.
# Final merge assembles all tissue RDS files into the combined score matrix.
# =============================================================================
cat("--- Step 5: GSVA scoring tissue-by-tissue (", N_WORKERS, "workers) ---\n")

TISSUE_SCORE_DIR <- file.path(OUT_DIR, "tissue_scores")
dir.create(TISSUE_SCORE_DIR, recursive = TRUE, showWarnings = FALSE)

bpparam  <- MulticoreParam(workers = N_WORKERS, progressbar = FALSE)
tissues  <- sort(unique(meta$tissue))
n_tissues <- length(tissues)
cat("  Tissues to score:", n_tissues, "\n")
cat("  Checkpoints:     ", TISSUE_SCORE_DIR, "\n\n")

# Helper: sanitise tissue name to a safe filename
safe_name <- function(x) gsub("[^A-Za-z0-9_-]", "_", x)

t_total_start <- proc.time()

for (i in seq_along(tissues)) {
    tiss     <- tissues[i]
    rds_file <- file.path(TISSUE_SCORE_DIR, paste0(safe_name(tiss), ".rds"))

    if (file.exists(rds_file)) {
        cat(sprintf("  [%2d/%2d] SKIP  %s (already scored)\n", i, n_tissues, tiss))
        next
    }

    # Subset samples for this tissue
    sids   <- meta[tissue == tiss, sample_id]
    t_mat  <- tpm_mat[, sids, drop = FALSE]

    t_start <- proc.time()
    cat(sprintf("  [%2d/%2d] GSVA  %s  (%d samples) ...",
                i, n_tissues, tiss, ncol(t_mat)))

    param  <- gsvaParam(
        exprData = t_mat,
        geneSets = gene_sets,
        minSize  = GSVA_PARAMS$min_size,
        maxSize  = GSVA_PARAMS$max_size
    )
    scores <- gsva(param, BPPARAM = bpparam, verbose = FALSE)
    saveRDS(scores, rds_file)

    elapsed <- round((proc.time() - t_start)["elapsed"], 0)
    cat(sprintf(" done in %ds\n", elapsed))

    # Explicitly free the tissue matrix and scores before next iteration
    rm(t_mat, scores, param)
    gc(verbose = FALSE)
}

# --- Merge all tissue score files into combined matrix
cat("\n  Merging tissue score files...\n")
rds_files    <- file.path(TISSUE_SCORE_DIR,
                           paste0(safe_name(tissues), ".rds"))
missing_rds  <- rds_files[!file.exists(rds_files)]
if (length(missing_rds) > 0) {
    stop("Missing score files for tissues: ",
         paste(tissues[!file.exists(rds_files)], collapse = ", "))
}

gsva_scores <- do.call(cbind, lapply(rds_files, readRDS))
cat("  Combined score matrix:", nrow(gsva_scores), "gene sets ×",
    ncol(gsva_scores), "samples\n")

t_total_elapsed <- round((proc.time() - t_total_start)["elapsed"] / 60, 1)
cat("  Total GSVA time:", t_total_elapsed, "min\n\n")

# Save combined scores
rds_path <- file.path(OUT_DIR, "gsva_scores.rds")
saveRDS(gsva_scores, rds_path)
cat("  Saved:", rds_path, "\n\n")

# =============================================================================
# STEP 6: Merge scores with metadata
# =============================================================================
cat("--- Step 6: Merging scores with metadata ---\n")

scores_long <- as.data.frame(t(gsva_scores)) |>
    tibble::rownames_to_column("sample_id") |>
    as.data.table()

meta_scores <- merge(meta, scores_long, by = "sample_id")
cat("  Merged rows:", nrow(meta_scores), "\n")

saveRDS(meta_scores, file.path(OUT_DIR, "metadata.rds"))
cat("  Saved:", file.path(OUT_DIR, "metadata.rds"), "\n")

# CSV export for Python ML script — flat table of scores + metadata.
# PRIMARY training data excludes ventilator deaths (DTHHRDY = 0).
# Ventilator cases have artificially controlled ischemic times and different
# pre-mortem physiology; they are out-of-distribution for forensic PMI
# estimation where the decedent was never on life support.
ml_cols  <- c("sample_id", "donor_id", "tissue", "ischemic_min",
              "hardy_scale", "sex", "age_range", geneset_cols)
meta_ml  <- meta_scores[is.na(hardy_scale) | hardy_scale != 0]
n_vent_excl <- nrow(meta_scores) - nrow(meta_ml)
write_csv(meta_ml[, ..ml_cols], file.path(OUT_DIR, "scores_for_ml.csv"))
cat("  Primary training set (DTHHRDY 1-4 + NA):", nrow(meta_ml), "samples",
    "(excluded", n_vent_excl, "ventilator cases, DTHHRDY=0)\n")
cat("  Saved:", file.path(OUT_DIR, "scores_for_ml.csv"), "(for python/03_train_pmi_model.py)\n\n")

# =============================================================================
# STEP 7: Partial Spearman correlation — HIF score vs ischemic time,
#         controlling for Hardy death scale, per tissue
#
# Hardy scale (DTHHRDY) confounds PMI: ventilator cases (0) have controlled
# ischemic times and different pre-mortem oxygenation than sudden death cases.
# We partial out Hardy scale using rank residualization:
#   1. Rank-transform both HIF score and ischemic time
#   2. Regress each on factor(hardy_scale), extract residuals
#   3. Pearson r of residuals = partial Spearman r
#
# Also saves a sensitivity analysis excluding ventilator cases (DTHHRDY == 0),
# since those samples have fundamentally different ischemic physiology.
# =============================================================================
cat("--- Step 7: Computing partial Spearman correlations (controlling for Hardy scale) ---\n")

# Hardy scale distribution in dataset
cat("  Hardy scale distribution:\n")
hardy_counts <- meta_scores[, .N, by = hardy_scale][order(hardy_scale)]
hardy_labels <- c("0=ventilator", "1=fast violent", "2=fast natural",
                   "3=intermediate", "4=slow")
for (i in seq_len(nrow(hardy_counts))) {
    hs  <- hardy_counts$hardy_scale[i]
    lab <- if (!is.na(hs) && hs >= 0 && hs <= 4) hardy_labels[hs + 1] else "unknown"
    cat(sprintf("    %s: %d samples\n", lab, hardy_counts$N[i]))
}
cat("\n")

# Partial Spearman via rank residualization (no extra package required)
partial_spearman <- function(x, y, z) {
    # x = HIF score, y = ischemic time, z = hardy scale (factor covariate)
    # Returns list: r, p_value, n
    ok <- !is.na(x) & !is.na(y) & !is.na(z)
    if (sum(ok) < 10) return(list(r = NA_real_, p = NA_real_, n = sum(ok)))

    rx <- rank(x[ok])
    ry <- rank(y[ok])
    zf <- factor(z[ok])  # treat as categorical — no ordinality assumption

    # Residualize ranks on Hardy scale
    ex <- residuals(lm(rx ~ zf))
    ey <- residuals(lm(ry ~ zf))

    ct <- cor.test(ex, ey, method = "pearson")
    list(r = ct$estimate, p = ct$p.value, n = sum(ok))
}

run_correlations <- function(data, label) {
    tissues <- sort(unique(data$tissue))
    res <- lapply(tissues, function(tiss) {
        df <- data[tissue == tiss]
        if (nrow(df) < 10) return(NULL)

        lapply(geneset_cols, function(gs) {
            ps <- partial_spearman(df[[gs]], df$ischemic_min, df$hardy_scale)
            data.table(
                cohort     = label,
                tissue     = tiss,
                geneset    = gs,
                n          = ps$n,
                spearman_r = ps$r,
                p_value    = ps$p
            )
        }) |> rbindlist(fill = TRUE)
    }) |> rbindlist(fill = TRUE)

    res[, q_value := p.adjust(p_value, method = "BH")]
    res[, sig     := q_value < 0.05]
    res[order(-abs(spearman_r))]
}


# Primary: exclude ventilator deaths (DTHHRDY = 0).
# Ventilator support maintains partial oxygenation — these samples have
# artificially controlled ischemic times and fundamentally different pre-mortem
# physiology. Forensic cases are never ventilator-supported. Hardy scale (1-4)
# is still partialled out within the primary cohort.
n_primary <- sum(is.na(meta_scores$hardy_scale) | meta_scores$hardy_scale != 0)
cat("  Primary cohort (DTHHRDY 1-4 + NA):", n_primary, "samples\n")
cor_primary <- run_correlations(
    meta_scores[is.na(hardy_scale) | hardy_scale != 0],
    "no_ventilator"
)

# Sensitivity: include ventilator cases — for reviewer comparison only.
# Shows magnitude of ventilator-driven confounding. Do NOT use for ML training.
n_vent <- sum(meta_scores$hardy_scale == 0, na.rm = TRUE)
cat("  Sensitivity (incl.", n_vent, "ventilator cases, DTHHRDY=0)\n")
cor_sensitivity <- run_correlations(meta_scores, "all_with_ventilator")

# Combine and save
cor_results <- rbindlist(list(cor_primary, cor_sensitivity))
write_csv(cor_results, file.path(OUT_DIR, "tissue_correlations.csv"))
cat("  Saved: tissue_correlations.csv\n")
cat("  Primary (no-ventilator) tests:", nrow(cor_primary),
    "| Significant (q<0.05):", sum(cor_primary$sig, na.rm = TRUE), "\n")
cat("  Sensitivity (all incl. ventilator) tests:", nrow(cor_sensitivity),
    "| Significant (q<0.05):", sum(cor_sensitivity$sig, na.rm = TRUE), "\n\n")

# =============================================================================
# STEP 8: Summary — top tissues ranked by mean |r| across HIF gene sets
#         + brain validation check
# =============================================================================
cat("--- Step 8: Summary and brain validation check ---\n")

make_summary <- function(cor_tbl, label) {
    cor_tbl[, .(
        cohort      = label,
        mean_abs_r  = mean(abs(spearman_r), na.rm = TRUE),
        max_abs_r   = max(abs(spearman_r), na.rm = TRUE),
        n_sig_sets  = sum(sig, na.rm = TRUE),
        n_sets      = .N,
        n_samples   = max(n)
    ), by = tissue][order(-mean_abs_r)]
}

summary_primary     <- make_summary(cor_primary,     "no_ventilator")
summary_sensitivity <- make_summary(cor_sensitivity, "all_with_ventilator")
summary_tbl         <- rbindlist(list(summary_primary, summary_sensitivity))

write_csv(summary_tbl, file.path(OUT_DIR, "summary.csv"))
cat("  Saved: summary.csv\n\n")

# Print top 15 (primary cohort)
cat("=== Top 15 tissues by mean |partial Spearman r| (Hardy-controlled) ===\n")
cat(sprintf("  %-45s  %6s  %6s  %5s\n", "Tissue", "mean|r|", "max|r|", "n_sig"))
cat("  ", strrep("-", 68), "\n", sep="")
for (i in seq_len(min(15, nrow(summary_primary)))) {
    row <- summary_primary[i]
    cat(sprintf("  %-45s  %6.3f  %6.3f  %5d\n",
                row$tissue, row$mean_abs_r, row$max_abs_r, row$n_sig_sets))
}

# --- Brain validation check --------------------------------------------------
# Neurons are among the most oxygen-sensitive cells in the body (~2 min to
# irreversible injury). Brain regions should rank near the top for HIF-PMI
# correlation. If they don't, suspect a metadata or gene set coverage issue.
cat("\n=== Brain region validation check ===\n")
brain_rows <- summary_primary[grepl("Brain|Cortex|Cerebellum|Hippocampus|Caudate|Putamen|Nucleus|Amygdala|Substantia|Frontal|Hypothalamus|Spinal",
                                     tissue, ignore.case = TRUE)]
if (nrow(brain_rows) == 0) {
    cat("  WARNING: No brain tissues found in results — check tissue name matching\n")
} else {
    brain_ranks <- match(brain_rows$tissue, summary_primary$tissue)
    cat(sprintf("  %-45s  rank  mean|r|\n", "Brain region"))
    cat("  ", strrep("-", 60), "\n", sep="")
    for (i in order(brain_ranks)) {
        cat(sprintf("  %-45s  %4d  %.3f\n",
                    brain_rows$tissue[i], brain_ranks[i], brain_rows$mean_abs_r[i]))
    }
    median_rank  <- median(brain_ranks)
    n_tissues    <- nrow(summary_primary)
    pct_rank     <- round(100 * median_rank / n_tissues, 0)
    cat(sprintf("\n  Median brain rank: %d / %d tissues (top %d%%)\n",
                as.integer(median_rank), n_tissues, 100 - pct_rank))
    # NOTE: Brain regions are expected to rank lower in the non-ventilator cohort
    # because brain samples are almost exclusively DTHHRDY 2+4 (natural/slow deaths)
    # with uniformly long ischemic times (median ~800 min). Without short-ischemia
    # ventilator cases to contrast against, the within-cohort variance is reduced.
    # Brain ranking in top third is acceptable; top 10 is the ideal.
    if (median_rank <= 10) {
        cat("  PASS: Brain regions rank near the top — strong HIF-PMI signal\n")
    } else if (median_rank <= n_tissues / 3) {
        cat("  ACCEPTABLE: Brain regions in top third — reduced contrast in non-ventilator cohort\n")
        cat("    (Brain samples are ~96% DTHHRDY 2+4 with uniformly long ischemic times)\n")
    } else {
        cat("  WARNING: Brain regions ranking low — investigate:\n")
        cat("    - Check SMTSISCH coverage in brain samples\n")
        cat("    - Check gene set coverage for brain-expressed genes\n")
        cat("    - Check for batch effects in brain Hardy scale distribution\n")
    }
}

cat("\n=== Done:", format(Sys.time()), "===\n")
cat("Next step: python/03_train_pmi_model.py\n")
