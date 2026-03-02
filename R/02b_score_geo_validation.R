# =============================================================================
# 02b_score_geo_validation.R
# HIForensics — GSVA scoring of GEO validation datasets
#
# MUST use the same scoring method (GSVA) and gene sets as 02_score_hif.R.
# Outputs a CSV per dataset in the same column format as scores_for_ml.csv
# so that 04_validate.py can apply the trained model directly.
#
# Run AFTER: scripts/02_download_geo_validation.sh (downloads GEO data)
# Run BEFORE: python/04_validate.py
#
# Output: data/validation/{accession}/scores_for_validation.csv
# =============================================================================

suppressPackageStartupMessages({
    library(GSVA)
    library(GEOquery)
    library(data.table)
    library(dplyr)
    library(readr)
    library(BiocParallel)
    library(limma)        # for normalizeBetweenArrays on array data
    library(org.Hs.eg.db) # Entrez ID → gene symbol fallback for Affymetrix arrays
    library(AnnotationDbi)
})

PROJECT   <- "/media/ross/New Volume1/HIForensics"
GEO_DIR   <- file.path(PROJECT, "data", "validation")
N_WORKERS <- 4    # MulticoreParam forks full parent; >4 workers OOMs on 64 GB

# Shared GSVA config — gene sets path, minSize/maxSize, expression filters
# Edit R/gsva_config.R to change any of these; never hardcode them here
source(file.path(PROJECT, "R", "gsva_config.R"))

cat("=== GEO validation scoring:", format(Sys.time()), "===\n\n")

gene_sets <- readRDS(GENESETS_RDS)
bpparam   <- MulticoreParam(workers = N_WORKERS, progressbar = FALSE)

# =============================================================================
# Dataset registry
# Each entry defines how to extract expression, PMI, and tissue from the GEO
# series matrix. Adjust column names if GEO updates the metadata fields.
# =============================================================================
DATASETS <- list(
    # PRIMARY validation dataset — bulk RNA-seq, PMI in minutes, frontal cortex
    GSE216281 = list(
        pmi_col    = "post mortem delay:ch1",  # values already in minutes (150–1245)
        tissue_col = "tissue:ch1",
        geo_tissue = "brain",
        notes      = "postmortem frontal cortex, Parkinson's + controls, NovaSeq RNA-seq"
    ),
    # SECONDARY — microarray (platform mismatch expected; kept for comparison)
    GSE53987 = list(
        pmi_col    = "pmi:ch1",
        tissue_col = "tissue:ch1",
        geo_tissue = "brain",
        notes      = "postmortem human brain, well-characterised PMI (Affymetrix microarray)"
    )
    # GSE45642 excluded: PMI not in GEO series matrix (Stanley Brain Bank supplementary table).
    # GSE116754 excluded: GPL13534 = 450K methylation array, hESC source (wrong dataset).
)

# =============================================================================
# Helper: probe → gene symbol mapping
# For array platforms: collapse multiple probes per gene by max mean expression
# =============================================================================
probes_to_genes <- function(expr_mat, feature_data) {
    # Primary: look for a gene symbol column by name pattern
    sym_col <- grep("gene.symbol|gene_symbol|symbol",
                    colnames(feature_data), ignore.case=TRUE, value=TRUE)[1]

    # Fallback: SPOT_ID contains Entrez gene IDs on some Affymetrix platforms
    # (e.g. GSE45642 GPL570 — probe IDs like "10000_at", SPOT_ID = 10000 = AKT3)
    if (is.na(sym_col) && "SPOT_ID" %in% colnames(feature_data)) {
        cat("  NOTE: No symbol column; mapping SPOT_ID (Entrez IDs) → gene symbols\n")
        entrez_ids <- as.character(feature_data[["SPOT_ID"]])
        mapped     <- suppressMessages(
            mapIds(org.Hs.eg.db, keys = entrez_ids,
                   column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
        )
        feature_data[["_symbol"]] <- mapped[entrez_ids]
        sym_col <- "_symbol"
        n_mapped <- sum(!is.na(feature_data[["_symbol"]]))
        cat("  Mapped", n_mapped, "/", nrow(feature_data), "probes to gene symbols\n")
    }

    if (is.na(sym_col)) {
        cat("  WARNING: no gene symbol column found in feature data\n")
        cat("  Available:", paste(colnames(feature_data)[1:10], collapse=", "), "\n")
        return(NULL)
    }

    symbols <- feature_data[[sym_col]]
    # Remove probes with empty/ambiguous symbols
    keep    <- nchar(trimws(symbols)) > 0 & !grepl("///", symbols)
    expr_mat <- expr_mat[keep, ]
    symbols  <- trimws(symbols[keep])

    # Collapse by max mean expression per symbol
    dt <- as.data.table(expr_mat)
    dt[, gene := symbols]
    dt <- dt[, lapply(.SD, mean), by = gene, .SDcols = colnames(expr_mat)]
    mat <- as.matrix(dt[, -1])
    rownames(mat) <- dt$gene
    mat
}

# =============================================================================
# Score one dataset
# =============================================================================
score_dataset <- function(acc, cfg) {
    cat(sprintf("--- %s: %s ---\n", acc, cfg$notes))
    out_dir <- file.path(GEO_DIR, acc)
    rds_path <- file.path(out_dir, paste0(acc, "_gse.rds"))
    out_csv  <- file.path(out_dir, "scores_for_validation.csv")

    if (file.exists(out_csv)) {
        cat("  [SKIP] scores_for_validation.csv already exists\n\n")
        return(invisible(NULL))
    }

    if (!file.exists(rds_path)) {
        cat("  [SKIP] GEO data not downloaded yet — run scripts/02_download_geo_validation.sh\n\n")
        return(invisible(NULL))
    }

    gse <- readRDS(rds_path)
    if (is.list(gse)) gse <- gse[[1]]  # getGEO returns a list

    # --- Expression matrix
    expr_mat <- exprs(gse)
    cat("  Raw matrix:", nrow(expr_mat), "features ×", ncol(expr_mat), "samples\n")

    # --- RNA-seq GEO datasets: expression data is in a supplementary file, not
    # the series matrix. Detect empty matrix and load supplementary counts.
    if (nrow(expr_mat) == 0) {
        cat("  Series matrix empty — looking for supplementary count file\n")
        supp_files <- list.files(out_dir,
                                 pattern = "counts.*\\.txt\\.gz$|tpm.*\\.txt\\.gz$",
                                 ignore.case = TRUE, full.names = TRUE)
        if (length(supp_files) == 0) {
            cat("  ERROR: no supplementary count file found in", out_dir, "\n\n")
            return(invisible(NULL))
        }
        supp_path <- supp_files[1]
        cat("  Loading:", basename(supp_path), "\n")
        counts_dt  <- fread(supp_path)
        gene_ids   <- counts_dt[[1]]          # first column = Ensembl IDs
        count_mat  <- as.matrix(counts_dt[, -1])

        # Map Ensembl IDs → gene symbols (strip version suffix first)
        gene_ids_clean <- sub("\\..*", "", gene_ids)
        symbols <- suppressMessages(
            mapIds(org.Hs.eg.db, keys = gene_ids_clean,
                   column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
        )
        mapped   <- !is.na(symbols)
        count_mat <- count_mat[mapped, ]
        symbols  <- symbols[mapped]
        n_mapped <- sum(mapped)
        cat("  Ensembl → symbol: mapped", n_mapped, "/", length(gene_ids), "genes\n")

        # Collapse duplicate symbols by highest mean expression
        mean_expr  <- rowMeans(count_mat)
        keep_idx   <- !duplicated(symbols[order(-mean_expr)])
        ord        <- order(-mean_expr)
        count_mat  <- count_mat[ord, ][keep_idx, ]
        symbols    <- symbols[ord][keep_idx]
        rownames(count_mat) <- symbols

        # Assign sample IDs positionally from pData (column order assumed to
        # match the series matrix sample order — standard GEO convention)
        pd_ids <- rownames(pData(gse))
        if (ncol(count_mat) == length(pd_ids)) {
            colnames(count_mat) <- pd_ids
        } else {
            cat("  WARNING: count matrix columns (", ncol(count_mat),
                ") != pData rows (", length(pd_ids), ") — cannot align\n\n")
            return(invisible(NULL))
        }
        expr_mat <- log2(count_mat + 1)
        cat("  Platform type: RNA-seq (supplementary file) — applied log2(count+1)\n")
        cat("  Final matrix:", nrow(expr_mat), "genes ×", ncol(expr_mat), "samples\n")

    # --- Standard path: series matrix contains expression data
    } else {
        # Determine if array or RNA-seq (RNA-seq has integer-like values)
        is_rnaseq <- median(expr_mat, na.rm=TRUE) > 100 &&
                     all(expr_mat == floor(expr_mat), na.rm=TRUE)

        if (!is_rnaseq) {
            # Array: normalize between arrays, then map probes → genes
            cat("  Platform type: microarray — normalising and mapping probes\n")
            expr_mat <- normalizeBetweenArrays(expr_mat, method="quantile")
            fd       <- fData(gse)
            expr_mat <- probes_to_genes(expr_mat, fd)
            if (is.null(expr_mat)) {
                cat("  ERROR: probe→gene mapping failed\n\n"); return(invisible(NULL))
            }
            cat("  After gene collapse:", nrow(expr_mat), "genes\n")
        } else {
            cat("  Platform type: RNA-seq — applying log2(count+1)\n")
            expr_mat <- log2(expr_mat + 1)
        }
    }

    # --- Metadata: extract PMI and tissue
    pd      <- pData(gse)
    pmi_raw <- pd[[cfg$pmi_col]]

    if (is.null(pmi_raw)) {
        # Try to find PMI column by pattern match
        pmi_col_found <- grep("pmi|postmortem.interval|ischemic",
                              colnames(pd), ignore.case=TRUE, value=TRUE)[1]
        if (!is.na(pmi_col_found)) {
            cat("  NOTE: PMI column '", cfg$pmi_col, "' not found; using '",
                pmi_col_found, "'\n", sep="")
            pmi_raw <- pd[[pmi_col_found]]
        } else {
            cat("  ERROR: PMI column not found. Available columns:\n")
            cat("  ", paste(colnames(pd)[1:20], collapse="\n  "), "\n\n")
            return(invisible(NULL))
        }
    }

    # Parse PMI — may be "12 hours", "720 min", "12.5", etc.
    pmi_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", pmi_raw)))

    # Print raw values before any unit conversion so you can sanity-check.
    # GSE45642 in particular has some very short PMIs that can pull the median
    # down and trigger a false hours→minutes conversion.
    cat("  Raw PMI values (pre-conversion):\n")
    cat("    n =", length(pmi_numeric), " | NAs =", sum(is.na(pmi_numeric)), "\n")
    cat("    min =", round(min(pmi_numeric, na.rm=TRUE), 2),
        " median =", round(median(pmi_numeric, na.rm=TRUE), 2),
        " max =", round(max(pmi_numeric, na.rm=TRUE), 2), "\n")
    cat("    first 10 raw values:", paste(round(head(pmi_numeric, 10), 2), collapse=", "), "\n")
    cat("    unique values (up to 20):",
        paste(sort(unique(round(pmi_numeric, 1)))[1:min(20, length(unique(pmi_numeric)))],
              collapse=", "), "\n")

    # Detect unit: if median < 50 assume hours, convert to minutes.
    # Verify this looks right from the raw print above before trusting it.
    pmi_median <- median(pmi_numeric, na.rm=TRUE)
    if (pmi_median < 50) {
        cat("  PMI unit: median =", round(pmi_median, 2),
            "< 50 → inferred HOURS, converting to minutes (×60)\n")
        cat("  If this is wrong (values are already minutes), set cfg$pmi_unit = 'min'\n")
        pmi_numeric <- pmi_numeric * 60
    } else {
        cat("  PMI unit: median =", round(pmi_median, 2),
            ">= 50 → inferred MINUTES, no conversion\n")
    }
    cat("  PMI range after conversion:", round(min(pmi_numeric, na.rm=TRUE), 0), "–",
        round(max(pmi_numeric, na.rm=TRUE), 0), "min\n")

    # Tissue
    tissue_raw <- if (!is.null(pd[[cfg$tissue_col]])) pd[[cfg$tissue_col]] else cfg$geo_tissue

    # --- Death-related metadata
    # GEO datasets don't use Hardy scale, but often record cause/manner of death.
    # Capture whatever is available so 04_validate.py can compare distributions
    # with GTEx training data (Hardy 0-4). Columns vary by study.
    death_patterns  <- "cause.of.death|manner.of.death|death|autopsy|hardy|cod"
    death_exclude   <- "contact_zip|postal_code|age.at.death|age_at_death"  # false positives
    death_cols      <- grep(death_patterns, colnames(pd), ignore.case=TRUE, value=TRUE)
    death_cols      <- death_cols[!grepl(death_exclude, death_cols, ignore.case=TRUE)]

    if (length(death_cols) > 0) {
        cat("  Death-related metadata columns found:", paste(death_cols, collapse=", "), "\n")
        death_meta <- pd[, death_cols, drop=FALSE]
        colnames(death_meta) <- death_cols   # keep original names
    } else {
        cat("  No death-related metadata columns found in pData\n")
        cat("  Available pData columns:", paste(colnames(pd)[1:15], collapse=", "), "\n")
        death_meta <- data.frame(row.names=rownames(pd))
    }

    # --- Clinical/disease metadata
    # Capture disease staging, diagnosis, and condition columns.
    # Critical for subsetting validation datasets (e.g. Braak-0 controls in GSE216281).
    clin_patterns <- "braak|lewy|parkinson|disease|diagnosis|condition|stage|status|phenotype|case|control|group|disorder|dementia|alzheimer"
    clin_exclude  <- "data_processing|extract_protocol|library|instrument|contact|submission|update"
    clin_cols     <- grep(clin_patterns, colnames(pd), ignore.case=TRUE, value=TRUE)
    clin_cols     <- clin_cols[!grepl(clin_exclude, clin_cols, ignore.case=TRUE)]
    # Avoid double-capturing columns already in death_meta
    clin_cols     <- setdiff(clin_cols, death_cols)

    if (length(clin_cols) > 0) {
        cat("  Clinical metadata columns found:", paste(clin_cols, collapse=", "), "\n")
        clin_meta <- pd[, clin_cols, drop=FALSE]
    } else {
        cat("  No clinical metadata columns found\n")
        clin_meta <- data.frame(row.names=rownames(pd))
    }

    meta_df <- data.table(
        sample_id    = colnames(expr_mat),
        pmi_min      = pmi_numeric,
        geo_tissue   = as.character(tissue_raw),
        gtex_tissue  = cfg$geo_tissue,   # coarse label for Python tissue mapper
        accession    = acc
    )

    # Merge death metadata (row-aligned to samples)
    if (ncol(death_meta) > 0) {
        death_dt <- as.data.table(death_meta, keep.rownames="sample_id")
        meta_df  <- merge(meta_df, death_dt, by="sample_id", all.x=TRUE)
    }

    # Merge clinical/disease metadata
    if (ncol(clin_meta) > 0) {
        clin_dt <- as.data.table(clin_meta, keep.rownames="sample_id")
        meta_df <- merge(meta_df, clin_dt, by="sample_id", all.x=TRUE)
    }
    meta_df <- meta_df[!is.na(pmi_min)]
    expr_mat <- expr_mat[, meta_df$sample_id, drop=FALSE]
    cat("  Samples with valid PMI:", nrow(meta_df), "\n")

    # --- GSVA scoring (same method as training)
    cat("  Running GSVA...\n")
    # Filter gene sets to genes present in this expression matrix
    gs_filtered <- lapply(gene_sets, function(g) intersect(g, rownames(expr_mat)))
    gs_filtered <- gs_filtered[lengths(gs_filtered) >= 5]
    cat("  Gene sets with >=5 genes in matrix:", length(gs_filtered), "\n")

    # Parameters from gsva_config.R — must match 02_score_hif.R exactly
    param  <- gsvaParam(exprData = expr_mat, geneSets = gs_filtered,
                        minSize = GSVA_PARAMS$min_size,
                        maxSize = GSVA_PARAMS$max_size)
    scores <- gsva(param, BPPARAM = bpparam, verbose = FALSE)

    # --- Merge and save
    scores_dt <- as.data.table(t(scores), keep.rownames="sample_id")
    out_dt    <- merge(meta_df, scores_dt, by="sample_id")

    write_csv(out_dt, out_csv)
    cat("  Saved:", out_csv, "\n")
    cat("  Columns:", ncol(out_dt), "| Rows:", nrow(out_dt), "\n\n")
}

for (acc in names(DATASETS)) {
    tryCatch(
        score_dataset(acc, DATASETS[[acc]]),
        error = function(e) cat("  ERROR in", acc, ":", conditionMessage(e), "\n\n")
    )
}

cat("=== Done:", format(Sys.time()), "===\n")
cat("Next step: python/04_validate.py\n")
