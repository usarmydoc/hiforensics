# =============================================================================
# 01_download_msigdb.R
# HIForensics Phase 1 — Curate HIF/hypoxia gene sets from MSigDB
#
# Uses msigdbr package (no account or download required).
# Saves curated gene sets as:
#   data/msigdb/hif_genesets.rds   — named list of character vectors
#   data/msigdb/hif_genesets.gmt   — GMT format (for GSVA / fgsea)
#   data/msigdb/hif_genesets_long.csv — tidy format (geneset, gene_symbol)
# =============================================================================

library(msigdbr)
library(dplyr)
library(readr)

PROJECT <- "/media/ross/New Volume1/HIForensics"
OUT_DIR <- file.path(PROJECT, "data", "msigdb")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== MSigDB gene set download:", format(Sys.time()), "===\n")
cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n\n")

# --- Define target gene sets -------------------------------------------------
# Each entry: collection (gs_cat) + optional subcollection (gs_subcat) + name pattern
TARGET_SETS <- list(
    # Hallmark — curated, low-redundancy
    list(cat = "H",  subcat = "",           pattern = "HALLMARK_HYPOXIA"),

    # C2: curated — Reactome HIF pathways
    list(cat = "C2", subcat = "CP:REACTOME", pattern = "HIF"),
    list(cat = "C2", subcat = "CP:REACTOME", pattern = "OXYGEN"),
    list(cat = "C2", subcat = "CP:REACTOME", pattern = "HYPOXIA"),

    # C2: curated — KEGG
    list(cat = "C2", subcat = "CP:KEGG_MEDICUS", pattern = "HIF"),

    # C2: curated — BioCarta
    list(cat = "C2", subcat = "CP:BIOCARTA", pattern = "HIF"),

    # C5: ontology — GO Biological Process
    list(cat = "C5", subcat = "GO:BP",       pattern = "HYPOXIA"),
    list(cat = "C5", subcat = "GO:BP",       pattern = "HIF"),
    list(cat = "C5", subcat = "GO:BP",       pattern = "PROLYL_HYDROXYLASE"),  # PHD enzymes

    # C6: oncogenic signatures (HIF in cancer context)
    list(cat = "C6", subcat = "",            pattern = "HIF")
)

# --- Fetch and filter --------------------------------------------------------
# Exclude sets that match the pattern string but are biologically off-target
EXCLUDE_SETS <- c(
    "GOBP_VIRAL_TRANSLATIONAL_FRAMESHIFTING",                          # false positive on "HIF" substring
    "REACTOME_ERYTHROCYTES_TAKE_UP_CARBON_DIOXIDE_AND_RELEASE_OXYGEN", # O2 transport, not HIF signalling
    "REACTOME_ERYTHROCYTES_TAKE_UP_OXYGEN_AND_RELEASE_CARBON_DIOXIDE"  # same
)

fetch_sets <- function(cat, subcat, pattern) {
    if (nchar(subcat) > 0) {
        df <- msigdbr(species = "Homo sapiens", collection = cat, subcollection = subcat)
    } else {
        df <- msigdbr(species = "Homo sapiens", collection = cat)
    }
    df |> filter(grepl(pattern, gs_name, ignore.case = FALSE),
                 !gs_name %in% EXCLUDE_SETS)
}

all_sets <- lapply(TARGET_SETS, function(x) {
    tryCatch(
        fetch_sets(x$cat, x$subcat, x$pattern),
        error = function(e) {
            cat("  WARNING:", x$cat, x$subcat, x$pattern, "->", conditionMessage(e), "\n")
            NULL
        }
    )
})

combined <- bind_rows(Filter(Negate(is.null), all_sets)) |>
    select(gs_name, gs_collection, gs_subcollection, gene_symbol, ensembl_gene) |>
    distinct()

cat("Gene sets retrieved:\n")
combined |>
    count(gs_name, gs_collection, gs_subcollection) |>
    arrange(gs_collection, gs_name) |>
    print(n = Inf)
cat("\n")

# --- Also add PHD enzyme proxies manually ------------------------------------
# PHD1/2/3 (EGLN1/2/3) activity: their targets reflect HIF stability
# VHL substrate: HIF1A, HIF2A (EPAS1) hydroxylation marks
phd_proxy_genes <- c("EGLN1", "EGLN2", "EGLN3",   # PHD1/2/3
                      "VHL",                          # E3 ligase
                      "HIF1A", "EPAS1",              # substrates
                      "EP300", "CREBBP",             # coactivators
                      "ARNT")                         # HIF1B dimerisation partner

phd_df <- tibble(
    gs_name          = "CUSTOM_PHD_VHL_AXIS",
    gs_collection    = "CUSTOM",
    gs_subcollection = "",
    gene_symbol      = phd_proxy_genes,
    ensembl_gene     = NA_character_
)

combined <- bind_rows(combined, phd_df)

# --- Save: tidy CSV ----------------------------------------------------------
tidy_path <- file.path(OUT_DIR, "hif_genesets_long.csv")
write_csv(combined, tidy_path)
cat("Saved tidy CSV:", tidy_path, "\n")

# --- Save: named list RDS ----------------------------------------------------
geneset_list <- combined |>
    group_by(gs_name) |>
    summarise(genes = list(unique(gene_symbol)), .groups = "drop") |>
    tibble::deframe()

rds_path <- file.path(OUT_DIR, "hif_genesets.rds")
saveRDS(geneset_list, rds_path)
cat("Saved RDS:     ", rds_path, "\n")
cat("  Gene sets:   ", length(geneset_list), "\n")
cat("  Size range:  ", min(lengths(geneset_list)), "–",
    max(lengths(geneset_list)), "genes\n\n")

# --- Save: GMT format --------------------------------------------------------
# GMT: one gene set per line — name <tab> description <tab> gene1 <tab> gene2 ...
gmt_path <- file.path(OUT_DIR, "hif_genesets.gmt")
gmt_lines <- lapply(names(geneset_list), function(gs) {
    genes <- geneset_list[[gs]]
    paste(c(gs, "HIForensics_curated", genes), collapse = "\t")
})
writeLines(unlist(gmt_lines), gmt_path)
cat("Saved GMT:     ", gmt_path, "\n")

# --- Summary -----------------------------------------------------------------
cat("\n=== Gene set summary ===\n")
cat(sprintf("  %-55s  %s\n", "Gene Set", "n genes"))
cat(strrep("-", 65), "\n")
for (nm in names(geneset_list)) {
    cat(sprintf("  %-55s  %d\n", nm, length(geneset_list[[nm]])))
}

cat("\n=== Done:", format(Sys.time()), "===\n")
cat("Next step: R/02_score_hif.R\n")
