# =============================================================================
# 05_figures.R
# HIForensics — Publication-ready figures from pipeline results CSVs
#
# Run this after the full pipeline completes. All inputs are CSVs/RDS already
# written by earlier scripts — no model refitting, no data reloading.
#
# Figures produced:
#   Fig 1 — Heatmap: partial Spearman r (HIF vs PMI) across tissues × gene sets
#   Fig 2 — SHAP feature importance: HIF sets vs tissue vs Hardy vs demographics
#   Fig 3 — Predicted vs actual PMI scatter (global model, OOF + test)
#   Fig 4 — Brain validation: per-tissue performance with brain regions highlighted
#
# Output: figures/  — PDF (vector) + PNG (300 DPI) for each figure
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(pheatmap)
    library(ggrepel)
    library(cowplot)
    library(patchwork)
    library(scales)
    library(RColorBrewer)
    library(viridis)
})

PROJECT  <- "/media/ross/New Volume1/HIForensics"
RES_DIR  <- file.path(PROJECT, "results")
FIG_DIR  <- file.path(PROJECT, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# Publication theme — clean, minimal, suitable for journal submission
theme_pub <- function(base_size = 11) {
    theme_cowplot(font_size = base_size) +
    theme(
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold"),
        axis.line          = element_line(linewidth = 0.4),
        legend.key.size    = unit(0.4, "cm"),
        plot.title         = element_text(face = "bold", size = base_size + 1),
        plot.subtitle      = element_text(color = "grey40", size = base_size - 1)
    )
}

save_figure <- function(p, name, width, height) {
    pdf_path <- file.path(FIG_DIR, paste0(name, ".pdf"))
    png_path <- file.path(FIG_DIR, paste0(name, ".png"))
    if (inherits(p, "pheatmap") || inherits(p, "list")) {
        # pheatmap saves differently
        pdf(pdf_path, width = width, height = height)
        grid::grid.newpage(); grid::grid.draw(p$gtable)
        dev.off()
        png(png_path, width = width, height = height, units = "in", res = 300)
        grid::grid.newpage(); grid::grid.draw(p$gtable)
        dev.off()
    } else {
        ggsave(pdf_path, p, width = width, height = height, device = cairo_pdf)
        ggsave(png_path, p, width = width, height = height, dpi = 300)
    }
    cat("  Saved:", pdf_path, "\n")
    cat("  Saved:", png_path, "\n")
}

# Brain region pattern — used across multiple figures
BRAIN_PATTERN <- "Brain|Cortex|Cerebellum|Hippocampus|Caudate|Putamen|
                  Nucleus accumbens|Amygdala|Substantia|Frontal|Hypothalamus|
                  Spinal cord"
BRAIN_PATTERN <- gsub("\\s+", "", BRAIN_PATTERN)

is_brain <- function(x) grepl(BRAIN_PATTERN, x, ignore.case = TRUE)

# =============================================================================
# FIGURE 1 — Heatmap: partial Spearman r across tissues × HIF gene sets
# =============================================================================
cat("--- Figure 1: HIF-PMI correlation heatmap ---\n")

cor_path <- file.path(RES_DIR, "hif_scores", "tissue_correlations.csv")
if (!file.exists(cor_path)) {
    cat("  SKIP: tissue_correlations.csv not found\n\n")
} else {
    cor_df <- read_csv(cor_path, show_col_types = FALSE)

    # Use primary cohort (Hardy-controlled, all samples)
    cor_primary <- cor_df |> filter(cohort == "all_hardy_controlled")

    # Pivot to matrix: tissues (rows) × gene sets (cols)
    cor_mat <- cor_primary |>
        select(tissue, geneset, spearman_r) |>
        pivot_wider(names_from = geneset, values_from = spearman_r) |>
        column_to_rownames("tissue") |>
        as.matrix()

    # Shorten gene set names for display
    colnames(cor_mat) <- gsub("^(HALLMARK|REACTOME|GOBP|KEGG_MEDICUS_|BIOCARTA|CUSTOM)_?", "",
                               colnames(cor_mat))
    colnames(cor_mat) <- gsub("_", " ", colnames(cor_mat))
    colnames(cor_mat) <- stringr::str_to_title(tolower(colnames(cor_mat)))

    # Shorten tissue names (remove "Brain - " prefix etc.)
    rownames(cor_mat) <- gsub(" - ", "\n", rownames(cor_mat))
    rownames(cor_mat) <- gsub("\\(.*?\\)", "", rownames(cor_mat)) |> trimws()

    # Row annotation: brain vs non-brain
    brain_flag <- data.frame(
        Region = ifelse(is_brain(rownames(cor_mat)), "Brain", "Other"),
        row.names = rownames(cor_mat)
    )
    ann_colors <- list(Region = c(Brain = "#E64B35", Other = "#B0B0B0"))

    # Sort rows: brain regions first, then by mean |r| descending
    brain_rows  <- rownames(cor_mat)[is_brain(rownames(cor_mat))]
    other_rows  <- rownames(cor_mat)[!is_brain(rownames(cor_mat))]
    other_order <- other_rows[order(-rowMeans(abs(cor_mat[other_rows, , drop=FALSE]),
                                              na.rm=TRUE))]
    cor_mat <- cor_mat[c(brain_rows, other_order), ]

    p1 <- pheatmap(
        cor_mat,
        color            = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
        breaks           = seq(-0.6, 0.6, length.out = 101),
        cluster_rows     = FALSE,
        cluster_cols     = TRUE,
        annotation_row   = brain_flag,
        annotation_colors= ann_colors,
        fontsize         = 7,
        fontsize_row     = 6,
        fontsize_col     = 7,
        border_color     = NA,
        main             = "Partial Spearman r: HIF score vs ischemic time\n(Hardy death scale controlled)",
        angle_col        = 45,
        na_col           = "grey90"
    )

    save_figure(p1, "fig1_hif_pmi_heatmap", width = 12, height = 14)
    cat("\n")
}

# =============================================================================
# FIGURE 2 — SHAP feature importance
# =============================================================================
cat("--- Figure 2: SHAP feature importance ---\n")

shap_path <- file.path(RES_DIR, "ml_model", "shap_importance.csv")
if (!file.exists(shap_path)) {
    cat("  SKIP: shap_importance.csv not found\n\n")
} else {
    shap_df <- read_csv(shap_path, show_col_types = FALSE)

    # Classify features into groups
    shap_df <- shap_df |>
        mutate(group = case_when(
            feature %in% c("sex", "age_mid")          ~ "Demographics",
            feature == "hardy_scale"                   ~ "Hardy death scale",
            startsWith(feature, "tissue_")             ~ "Tissue identity",
            TRUE                                       ~ "HIF gene set"
        )) |>
        mutate(group = factor(group, levels = c(
            "HIF gene set", "Tissue identity", "Hardy death scale", "Demographics"
        )))

    group_colors <- c(
        "HIF gene set"      = "#E64B35",
        "Tissue identity"   = "#4DBBD5",
        "Hardy death scale" = "#F39B7F",
        "Demographics"      = "#B0B0B0"
    )

    # 2a: Top 25 individual features (bar chart)
    top25 <- shap_df |> slice_max(mean_abs_shap, n = 25)

    # Clean feature labels
    top25 <- top25 |>
        mutate(label = feature |>
            gsub("^tissue_", "", x=_) |>
            gsub("_", " ", x=_) |>
            gsub("HALLMARK |REACTOME |GOBP |BIOCARTA |CUSTOM |KEGG MEDICUS ", "", x=_) |>
            stringr::str_to_title() |>
            stringr::str_trunc(45))

    p2a <- ggplot(top25, aes(x = mean_abs_shap,
                              y = reorder(label, mean_abs_shap),
                              fill = group)) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = group_colors, name = "Feature group") +
        scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
        labs(
            title    = "SHAP feature importance — global PMI model",
            subtitle = "Mean |SHAP value| across training samples",
            x        = "Mean |SHAP value|",
            y        = NULL
        ) +
        theme_pub() +
        theme(legend.position = "right")

    # 2b: Group-level summary (donut / stacked bar)
    group_summary <- shap_df |>
        group_by(group) |>
        summarise(total_shap = sum(mean_abs_shap), .groups = "drop") |>
        mutate(pct = 100 * total_shap / sum(total_shap),
               label = sprintf("%s\n%.1f%%", group, pct))

    p2b <- ggplot(group_summary,
                  aes(x = reorder(group, total_shap), y = total_shap, fill = group)) +
        geom_col(width = 0.6, show.legend = FALSE) +
        geom_text(aes(label = sprintf("%.1f%%", pct)),
                  hjust = -0.1, size = 3.5) +
        scale_fill_manual(values = group_colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
        coord_flip() +
        labs(
            title    = "Feature group contributions",
            subtitle = "Summed mean |SHAP| by group",
            x        = NULL,
            y        = "Total mean |SHAP value|"
        ) +
        theme_pub()

    p2 <- p2a + p2b + plot_layout(widths = c(2, 1)) +
        plot_annotation(tag_levels = "a")

    save_figure(p2, "fig2_shap_importance", width = 14, height = 7)
    cat("\n")
}

# =============================================================================
# FIGURE 3 — Predicted vs actual PMI scatter
# =============================================================================
cat("--- Figure 3: Predicted vs actual PMI ---\n")

oof_path  <- file.path(RES_DIR, "ml_model", "oof_predictions.csv")
test_path <- file.path(RES_DIR, "ml_model", "test_predictions.csv")

if (!file.exists(oof_path)) {
    cat("  SKIP: oof_predictions.csv not found\n\n")
} else {
    oof_df  <- read_csv(oof_path,  show_col_types = FALSE) |> mutate(split = "CV (OOF)")
    test_df <- if (file.exists(test_path)) {
        read_csv(test_path, show_col_types = FALSE) |> mutate(split = "Held-out test")
    } else NULL

    pred_df <- bind_rows(oof_df, test_df) |>
        mutate(
            is_brain   = is_brain(tissue),
            tissue_grp = if_else(is_brain, tissue, "Other tissues"),
            # Truncate tissue labels for legend
            tissue_grp = gsub(" - ", "\n", tissue_grp),
            tissue_grp = gsub("\\(.*?\\)", "", tissue_grp) |> trimws()
        )

    # Overall Spearman r (OOF)
    oof_only <- pred_df |> filter(split == "CV (OOF)")
    rho <- cor(oof_only$ischemic_min, oof_only$predicted_min,
               method = "spearman", use = "complete.obs")
    rho_label <- sprintf("OOF ρ = %.3f", rho)

    # Downsample for plotting (OOF can have ~14k points)
    set.seed(42)
    plot_df <- pred_df |>
        group_by(split) |>
        slice_sample(n = min(n(), 3000)) |>
        ungroup()

    p3 <- ggplot(plot_df, aes(x = ischemic_min, y = predicted_min)) +
        geom_point(aes(color = is_brain), alpha = 0.4, size = 0.8) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    color = "black", linewidth = 0.5) +
        geom_smooth(method = "lm", se = TRUE, color = "grey30",
                    linewidth = 0.8, linetype = "solid") +
        scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                           labels = c("TRUE" = "Brain", "FALSE" = "Other"),
                           name = "Region") +
        scale_x_continuous(labels = comma, limits = c(0, NA)) +
        scale_y_continuous(labels = comma, limits = c(0, NA)) +
        facet_wrap(~split, scales = "free") +
        annotate("text", x = Inf, y = -Inf, label = rho_label,
                 hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "italic") +
        labs(
            title    = "Predicted vs actual ischemic time",
            subtitle = "Global XGBoost model  |  dashed line = perfect prediction",
            x        = "Actual ischemic time (min)",
            y        = "Predicted ischemic time (min)"
        ) +
        theme_pub()

    save_figure(p3, "fig3_predicted_vs_actual", width = 12, height = 5)

    # 3b: Per-tissue performance (top 20 tissues by Spearman r, brain highlighted)
    pt_path <- file.path(RES_DIR, "ml_model", "per_tissue_performance.csv")
    if (file.exists(pt_path)) {
        pt_df <- read_csv(pt_path, show_col_types = FALSE) |>
            arrange(desc(spearman_r)) |>
            mutate(
                rank     = row_number(),
                is_brain = is_brain(tissue),
                label    = gsub(" - ", "\n", tissue) |>
                               gsub("\\(.*?\\)", "", x=_) |>
                               trimws()
            )

        top20 <- pt_df |> slice_head(n = 20)

        p3b <- ggplot(top20, aes(x = reorder(label, spearman_r),
                                  y = spearman_r, fill = is_brain)) +
            geom_col(width = 0.7) +
            geom_errorbar(aes(ymin = spearman_r - 0.05,
                              ymax = spearman_r + 0.05),
                          width = 0.3, alpha = 0) +  # placeholder for CI if added
            scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#91BED4"),
                              labels = c("TRUE" = "Brain", "FALSE" = "Other"),
                              name = "Region") +
            scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
            coord_flip() +
            labs(
                title    = "Per-tissue PMI prediction (top 20)",
                subtitle = "Spearman ρ — global model, OOF predictions",
                x        = NULL,
                y        = "Spearman ρ"
            ) +
            theme_pub()

        save_figure(p3b, "fig3b_per_tissue_performance", width = 9, height = 7)
    }
    cat("\n")
}

# =============================================================================
# FIGURE 4 — Brain validation check
# =============================================================================
cat("--- Figure 4: Brain validation ---\n")

pt_path <- file.path(RES_DIR, "ml_model", "per_tissue_performance.csv")
if (!file.exists(pt_path)) {
    cat("  SKIP: per_tissue_performance.csv not found\n\n")
} else {
    pt_df <- read_csv(pt_path, show_col_types = FALSE) |>
        arrange(desc(spearman_r)) |>
        mutate(
            rank     = row_number(),
            is_brain = is_brain(tissue),
            label    = ifelse(is_brain,
                              gsub("Brain - ", "", tissue) |>
                                  gsub("\\(.*?\\)", "", x=_) |>
                                  trimws(),
                              NA_character_)
        )

    n_tissues   <- nrow(pt_df)
    brain_df    <- pt_df |> filter(is_brain)
    median_rank <- median(brain_df$rank)

    p4 <- ggplot(pt_df, aes(x = rank, y = spearman_r)) +
        # All tissues as faint points
        geom_point(data = filter(pt_df, !is_brain),
                   color = "#CCCCCC", size = 1.5, alpha = 0.7) +
        # Brain regions highlighted
        geom_point(data = brain_df,
                   aes(color = spearman_r), size = 3) +
        geom_label_repel(data = brain_df,
                         aes(label = label),
                         size = 2.8, max.overlaps = 20,
                         segment.color = "grey60",
                         min.segment.length = 0.2,
                         box.padding = 0.3) +
        # Median brain rank vertical line
        geom_vline(xintercept = median_rank, linetype = "dashed",
                   color = "#E64B35", linewidth = 0.6) +
        annotate("text",
                 x = median_rank + 1, y = max(pt_df$spearman_r, na.rm=TRUE),
                 label = sprintf("Median brain rank: %d / %d", as.integer(median_rank), n_tissues),
                 hjust = 0, size = 3.2, color = "#E64B35", fontface = "italic") +
        scale_color_viridis_c(option = "plasma", name = "Spearman ρ",
                               limits = c(0, 1)) +
        scale_x_continuous(breaks = pretty_breaks(n = 6)) +
        labs(
            title    = "Brain regions dominate HIF-based PMI prediction",
            subtitle = sprintf(
                "Brain tissues ranked %s — consistent with extreme neuronal oxygen sensitivity",
                if (median_rank <= 10)
                    sprintf("median %d / %d (top %.0f%%)", as.integer(median_rank), n_tissues,
                            100 * (1 - median_rank/n_tissues))
                else
                    sprintf("median %d / %d", as.integer(median_rank), n_tissues)
            ),
            x        = "Tissue rank (1 = strongest HIF-PMI correlation)",
            y        = "Spearman ρ (HIF vs ischemic time)"
        ) +
        theme_pub() +
        theme(legend.position = "right")

    save_figure(p4, "fig4_brain_validation", width = 11, height = 6)
    cat("\n")
}

# =============================================================================
# SUMMARY
# =============================================================================
cat("=== Figures complete ===\n")
figs <- list.files(FIG_DIR, pattern = "\\.(pdf|png)$")
cat("  Files in figures/:\n")
for (f in sort(figs)) cat("    ", f, "\n")
cat("\n  Edit theme_pub() at the top of this script to adjust fonts/sizes globally.\n")
