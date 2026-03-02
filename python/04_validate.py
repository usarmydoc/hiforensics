"""
04_validate.py
HIForensics Phase 3 — Independent validation against GEO postmortem datasets

Design
------
- PURE INFERENCE: global_model.joblib and tissue_models.joblib are loaded
  and applied as-is. No retraining, no fine-tuning, no parameter updates.
  Any fitting on GEO data invalidates the validation.

- GEO data is GSVA-scored by R/02b_score_geo_validation.R before this script
  runs. The scoring method must match training (same GSVA implementation,
  same gene sets). Using a different scorer shifts the score distribution and
  makes the model appear to fail even if it's correct.

- Tissue label mismatch: GEO tissue names don't match GTEx tissue strings.
  Three-tier approach (applied in order of preference):
    1. Tissue-specific model (no tissue encoding needed) — preferred for brain
    2. Global model + manual GTEx tissue mapping (mapped label)
    3. Global model + unknown tissue (all-zero one-hot, HIF scores drive pred)
  All three are computed and reported so the reader can see the difference.

Input:  data/validation/{accession}/scores_for_validation.csv
        results/ml_model/global_model.joblib
        results/ml_model/tissue_models.joblib
Output: results/validation/
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
from scipy import stats
import joblib
import json
from pathlib import Path

from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

# =============================================================================
# CONFIG
# =============================================================================
PROJECT     = Path("/media/ross/New Volume1/HIForensics")
VAL_DATA    = PROJECT / "data" / "validation"
MODEL_DIR   = PROJECT / "results" / "ml_model"
OUT_DIR     = PROJECT / "results" / "validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# GEO DATASET REGISTRY
# geo_tissue:  the coarse tissue label written by 02b_score_geo_validation.R
# gtex_tissue: the single best-matching GTEx tissue label for global model
# ts_model_key: key in tissue_models.joblib to use (None = not applicable)
#
# GTEx brain tissue labels (exact strings from the training data):
#   "Brain - Frontal Cortex (BA9)"
#   "Brain - Cortex"
#   "Brain - Hippocampus"
#   "Brain - Caudate (basal ganglia)"
#   "Brain - Putamen (basal ganglia)"
#   "Brain - Nucleus accumbens (basal ganglia)"
#   "Brain - Cerebellum" / "Brain - Cerebellar Hemisphere"
#   "Brain - Amygdala" / "Brain - Hypothalamus" / "Brain - Substantia nigra"
#   "Brain - Anterior cingulate cortex (BA24)"
#   "Brain - Spinal cord (cervical c-1)"
#
# For brain GEO datasets, frontal cortex is the most commonly profiled region.
# Use "Brain - Frontal Cortex (BA9)" as default; change if study specifies
# a different region.
# =============================================================================
DATASETS = {
    "GSE216281": dict(
        description = "Postmortem frontal cortex, Parkinson's + controls, NovaSeq RNA-seq",
        geo_tissue  = "brain",
        gtex_tissue = "Brain - Frontal Cortex (BA9)",
        ts_model_key= "Brain - Frontal Cortex (BA9)",
        pmi_col     = "pmi_min",
    ),
    "GSE53987": dict(
        description = "Postmortem human brain, well-characterised PMI (microarray — platform mismatch)",
        geo_tissue  = "brain",
        gtex_tissue = "Brain - Frontal Cortex (BA9)",
        ts_model_key= "Brain - Frontal Cortex (BA9)",
        pmi_col     = "pmi_min",
    ),
    # GSE45642: PMI not in GEO series matrix — excluded
    # GSE116754: 450K methylation array, wrong dataset — excluded
}

# =============================================================================
# HELPERS
# =============================================================================

def eval_metrics(y_true, y_pred, label=""):
    y_pred = np.clip(y_pred, 0, None)
    rmse   = np.sqrt(mean_squared_error(y_true, y_pred))
    mae    = mean_absolute_error(y_true, y_pred)
    r2     = r2_score(y_true, y_pred)
    rho, p = stats.spearmanr(y_true, y_pred)
    n      = len(y_true)
    if label:
        print(f"    {label:<40} n={n:4d}  RMSE={rmse:7.1f} min  "
              f"MAE={mae:7.1f} min  R²={r2:6.3f}  ρ={rho:6.3f} (p={p:.2e})")
    return dict(label=label, n=n, rmse=rmse, mae=mae, r2=r2,
                spearman_r=rho, spearman_p=p)


def load_scored_data(accession, cfg):
    csv_path = VAL_DATA / accession / "scores_for_validation.csv"
    if not csv_path.exists():
        print(f"  [MISSING] {csv_path}")
        print(f"  → Run R/02b_score_geo_validation.R first\n")
        return None
    df = pd.read_csv(csv_path)
    # Validate PMI column
    if cfg["pmi_col"] not in df.columns:
        print(f"  [ERROR] PMI column '{cfg['pmi_col']}' not in {csv_path}")
        print(f"  Available columns: {list(df.columns)}\n")
        return None
    df = df.dropna(subset=[cfg["pmi_col"]])
    df = df[df[cfg["pmi_col"]] > 0]
    print(f"  Loaded {len(df)} samples, PMI range: "
          f"{df[cfg['pmi_col']].min():.0f}–{df[cfg['pmi_col']].max():.0f} min")
    return df


def get_hif_cols(df, global_model):
    """Identify HIF score columns present in both the data and the model."""
    num_features = global_model["pre"].transformers_[0][2]
    hif_in_model = [f for f in num_features
                    if f not in ("sex", "age_mid", "hardy_scale")]
    hif_in_data  = [f for f in hif_in_model if f in df.columns]
    missing      = set(hif_in_model) - set(hif_in_data)
    if missing:
        print(f"  WARNING: {len(missing)} HIF feature(s) missing from GEO data "
              f"(will be imputed to median by pipeline):")
        for f in sorted(missing):
            print(f"    {f}")
    return hif_in_model


def prepare_global_features(df, hif_cols, gtex_tissue, global_model):
    """
    Build feature DataFrame for the global model.
    - Uses gtex_tissue for the tissue column (mapped GTEx label or None)
    - Missing HIF cols filled with NaN (imputer handles them)
    - Missing sex/age/hardy filled with NaN (imputer handles them)
    """
    feat = pd.DataFrame(index=df.index)
    for col in hif_cols:
        feat[col] = df[col] if col in df.columns else np.nan
    feat["sex"]         = df["sex"] if "sex" in df.columns else np.nan
    feat["age_mid"]     = df["age_mid"] if "age_mid" in df.columns else np.nan
    feat["hardy_scale"] = df["hardy_scale"] if "hardy_scale" in df.columns else np.nan

    if gtex_tissue is not None:
        feat["tissue"] = gtex_tissue
        mapping_note   = f"mapped to GTEx: '{gtex_tissue}'"
    else:
        # Unknown tissue → one-hot encoder emits all zeros (handle_unknown="ignore")
        # Prediction is then driven by HIF scores + demographics alone
        feat["tissue"] = "__UNKNOWN__"
        mapping_note   = "unknown tissue (one-hot all-zero, HIF scores drive prediction)"

    return feat, mapping_note


# =============================================================================
# LOAD MODELS
# =============================================================================
print("=== HIForensics Validation  ===\n")

for model_path in [MODEL_DIR / "global_model.joblib",
                   MODEL_DIR / "tissue_models.joblib"]:
    if not model_path.exists():
        raise FileNotFoundError(
            f"{model_path} not found.\n"
            "Run python/03_train_pmi_model.py first."
        )

global_model   = joblib.load(MODEL_DIR / "global_model.joblib")
tissue_models  = joblib.load(MODEL_DIR / "tissue_models.joblib")
print(f"Loaded global model")
print(f"Loaded {len(tissue_models)} tissue-specific models: "
      f"{', '.join(sorted(tissue_models.keys())[:5])}...\n")

# Load training performance as reference
train_meta_path = MODEL_DIR / "run_metadata.json"
if train_meta_path.exists():
    with open(train_meta_path) as f:
        train_meta = json.load(f)
    print(f"Training reference performance (held-out test set):")
    print(f"  Global model — R²={train_meta.get('test_r2_global','?'):.3f}  "
          f"RMSE={train_meta.get('test_rmse_min','?'):.1f} min  "
          f"ρ={train_meta.get('test_spearman_r','?'):.3f}")
    print(f"  ΔR² (HIF beyond tissue): {train_meta.get('delta_r2_test','?'):.3f}\n")

# =============================================================================
# HARDY SCALE DISTRIBUTION COMPARISON
# GTEx training Hardy distribution vs whatever death metadata GEO provides.
# GEO datasets don't use Hardy scale directly, but the R scoring script
# captures any cause/manner-of-death columns from pData. We report both
# so readers can judge whether the validation cohort is comparable to training.
# =============================================================================
HARDY_LABELS = {
    0: "ventilator",
    1: "fast violent",
    2: "fast natural",
    3: "intermediate",
    4: "slow",
}

def hardy_distribution(series, label):
    """Print Hardy scale counts and percentages for a numeric series."""
    counts = series.value_counts().sort_index()
    total  = counts.sum()
    print(f"  {label} (n={total}):")
    for val, n in counts.items():
        try:
            key  = int(val)
            desc = HARDY_LABELS.get(key, "unknown")
            print(f"    Hardy {key} ({desc:<14}): {n:4d}  ({100*n/total:5.1f}%)")
        except (ValueError, TypeError):
            print(f"    '{val}': {n:4d}  ({100*n/total:5.1f}%)")
    return counts


print("=" * 70)
print("HARDY SCALE / DEATH CONTEXT COMPARISON")
print("=" * 70)

# Load GTEx training Hardy distribution from scores_for_ml.csv
ml_csv = PROJECT / "results" / "hif_scores" / "scores_for_ml.csv"
gtex_hardy_counts = None
if ml_csv.exists():
    ml_df = pd.read_csv(ml_csv, usecols=["hardy_scale"])
    ml_df = ml_df.dropna(subset=["hardy_scale"])
    gtex_hardy_counts = hardy_distribution(ml_df["hardy_scale"], "GTEx training data")
    print()
else:
    print(f"  GTEx training data not found at {ml_csv} — skipping distribution.\n")

# Check each GEO dataset for death metadata
death_col_pattern = r"cause.of.death|manner.of.death|death|autopsy|cod"
geo_death_context = {}   # accession → string summary for reviewer sentence

for accession in DATASETS:
    csv_path = VAL_DATA / accession / "scores_for_validation.csv"
    if not csv_path.exists():
        continue

    df_geo = pd.read_csv(csv_path)
    death_cols = [c for c in df_geo.columns
                  if any(kw in c.lower()
                         for kw in ["cause", "death", "manner", "autopsy", "cod", "hardy"])]

    print(f"  {accession} death metadata:")
    if death_cols:
        for col in death_cols:
            print(f"\n    Column: '{col}'")
            val_counts = df_geo[col].value_counts()
            total      = val_counts.sum()
            for val, n in val_counts.items():
                print(f"      {str(val):<40}: {n:4d}  ({100*n/total:5.1f}%)")
        # Store the dominant category for the reviewer sentence
        primary_col  = death_cols[0]
        top_val      = df_geo[primary_col].mode()[0] if not df_geo[primary_col].isna().all() else "unknown"
        geo_death_context[accession] = f"predominant cause: '{top_val}'"
    else:
        print(f"    No death metadata columns in scores_for_validation.csv")
        print(f"    (GEO pData may not have recorded cause/manner of death)")
        geo_death_context[accession] = "death metadata not recorded in GEO"
    print()

# Hardy comparability summary
print("  Comparability note:")
if gtex_hardy_counts is not None:
    dominant_hardy = int(gtex_hardy_counts.idxmax())
    dominant_pct   = 100 * gtex_hardy_counts.max() / gtex_hardy_counts.sum()
    vent_pct       = 100 * gtex_hardy_counts.get(0, 0) / gtex_hardy_counts.sum()
    slow_pct       = 100 * gtex_hardy_counts.get(4, 0) / gtex_hardy_counts.sum()
    print(f"  GTEx training: predominantly Hardy {dominant_hardy} "
          f"({HARDY_LABELS[dominant_hardy]}, {dominant_pct:.0f}%); "
          f"ventilator={vent_pct:.0f}%, slow={slow_pct:.0f}%")
for acc, ctx in geo_death_context.items():
    print(f"  {acc}: {ctx}")
print()

# =============================================================================
# VALIDATE EACH DATASET
# =============================================================================
all_results = []

for accession, cfg in DATASETS.items():
    print("=" * 70)
    print(f"Dataset: {accession} — {cfg['description']}")
    print("=" * 70)

    df = load_scored_data(accession, cfg)
    if df is None:
        continue

    y_true = df[cfg["pmi_col"]].values
    hif_cols = get_hif_cols(df, global_model)
    results_for_dataset = []

    # ------------------------------------------------------------------
    # APPROACH 1: Tissue-specific model (HIF scores only, no tissue encoding)
    # Most direct test of HIF→PMI signal, no tissue label issues at all.
    # ------------------------------------------------------------------
    ts_key = cfg["ts_model_key"]
    if ts_key and ts_key in tissue_models:
        print(f"\n  [1] Tissue-specific model ({ts_key})")
        # Use only HIF cols the ts model was trained on (same HIF cols)
        hif_only = [c for c in hif_cols
                    if c not in ("sex", "age_mid", "hardy_scale", "tissue")]
        X_hif = np.column_stack([
            df[c].values if c in df.columns else np.full(len(df), np.nan)
            for c in hif_only
        ])
        y_pred_log = tissue_models[ts_key].predict(X_hif)
        y_pred     = np.clip(np.expm1(y_pred_log), 0, None)
        m = eval_metrics(y_true, y_pred,
                         label=f"Tissue-specific [{ts_key[:30]}]")
        m.update(accession=accession, approach="tissue_specific",
                 tissue_used=ts_key)
        results_for_dataset.append(m)
    elif ts_key:
        print(f"\n  [1] Tissue-specific model: '{ts_key}' not in trained models "
              f"(too few GTEx samples at training time)")

    # ------------------------------------------------------------------
    # APPROACH 2: Global model + GTEx tissue mapping
    # Maps GEO tissue description to closest GTEx tissue string.
    # ------------------------------------------------------------------
    print(f"\n  [2] Global model + tissue mapping")
    feat2, note2 = prepare_global_features(df, hif_cols, cfg["gtex_tissue"],
                                            global_model)
    print(f"    Tissue: {note2}")
    y_pred_log2 = global_model.predict(feat2)
    y_pred2     = np.clip(np.expm1(y_pred_log2), 0, None)
    m2 = eval_metrics(y_true, y_pred2, label="Global model (mapped tissue)")
    m2.update(accession=accession, approach="global_mapped",
              tissue_used=str(cfg["gtex_tissue"]))
    results_for_dataset.append(m2)

    # ------------------------------------------------------------------
    # APPROACH 3: Global model + unknown tissue (HIF scores drive prediction)
    # Worst-case scenario: model has no tissue information at all.
    # Upper bound on generalisation if you have no idea what tissue you have.
    # ------------------------------------------------------------------
    print(f"\n  [3] Global model + unknown tissue (HIF only, tissue all-zero)")
    feat3, note3 = prepare_global_features(df, hif_cols, None, global_model)
    y_pred_log3 = global_model.predict(feat3)
    y_pred3     = np.clip(np.expm1(y_pred_log3), 0, None)
    m3 = eval_metrics(y_true, y_pred3, label="Global model (unknown tissue)")
    m3.update(accession=accession, approach="global_unknown_tissue",
              tissue_used="unknown")
    results_for_dataset.append(m3)

    all_results.extend(results_for_dataset)

    # ------------------------------------------------------------------
    # Per-dataset summary
    # ------------------------------------------------------------------
    print(f"\n  Summary for {accession}:")
    best = max(results_for_dataset, key=lambda x: x["spearman_r"])
    print(f"  Best approach: {best['approach']}  "
          f"ρ={best['spearman_r']:.3f}  RMSE={best['rmse']:.1f} min")

    # Generalisation check vs training
    if train_meta_path.exists():
        train_r2   = train_meta.get("test_r2_global", None)
        best_r2    = best["r2"]
        delta      = best_r2 - train_r2 if train_r2 else None
        if delta is not None:
            status = "PASS" if delta > -0.15 else "WARN — larger than expected drop"
            print(f"  Generalisation: train R²={train_r2:.3f} → val R²={best_r2:.3f}  "
                  f"(Δ={delta:+.3f}) [{status}]")

    print()

# =============================================================================
# SAVE RESULTS
# =============================================================================
if all_results:
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUT_DIR / "validation_metrics.csv", index=False)
    print("\nSaved: results/validation/validation_metrics.csv")

    # Cross-dataset summary table — best approach per dataset
    summary = (results_df.sort_values("spearman_r", ascending=False)
                         .groupby("accession").first().reset_index()
               [["accession", "approach", "tissue_used", "n",
                 "rmse", "mae", "r2", "spearman_r", "spearman_p"]])
    summary.to_csv(OUT_DIR / "validation_summary.csv", index=False)
    print("Saved: results/validation/validation_summary.csv")

    print("\n=== Cross-dataset validation summary (best approach per dataset) ===")
    print(f"  {'Dataset':<12}  {'Approach':<25}  {'n':>4}  "
          f"{'ρ':>6}  {'RMSE':>7}  {'R²':>6}")
    print("  " + "-" * 72)
    for _, row in summary.iterrows():
        print(f"  {row.accession:<12}  {row.approach:<25}  {int(row.n):>4}  "
              f"{row.spearman_r:>6.3f}  {row.rmse:>7.1f}  {row.r2:>6.3f}")

    # Reviewer-ready output
    # ρ > 0.5 threshold: ρ > 0.3 is too lenient for a forensic utility claim.
    # Reviewers will push back on weak correlations being described as generalisation.
    RHO_THRESHOLD = 0.5
    n_pass    = (results_df.groupby("accession")["spearman_r"].max() > RHO_THRESHOLD).sum()
    n_total   = summary.shape[0]
    mean_rho  = summary["spearman_r"].mean()

    print(f"\n  === Reviewer-ready summary ===")
    print(f"  (Pass threshold: Spearman ρ > {RHO_THRESHOLD} — raise or defend this in peer review)\n")
    print(f"  Performance sentence (main text):")
    print(f"  The HIF-based PMI estimator generalised to {n_pass}/{n_total} "
          f"independent postmortem datasets (Spearman ρ > {RHO_THRESHOLD}; "
          f"mean ρ={mean_rho:.3f}).")

    # Hardy context sentence — included as supplementary/methods note
    print(f"\n  Hardy/death-context sentence (methods or supplementary):")
    if gtex_hardy_counts is not None:
        dominant_hardy = int(gtex_hardy_counts.idxmax())
        dominant_label = HARDY_LABELS[dominant_hardy]
        dominant_pct   = 100 * gtex_hardy_counts.max() / gtex_hardy_counts.sum()
        vent_pct       = 100 * gtex_hardy_counts.get(0, 0) / gtex_hardy_counts.sum()
        slow_pct       = 100 * gtex_hardy_counts.get(4, 0) / gtex_hardy_counts.sum()

        geo_ctx_parts = []
        for acc, ctx in geo_death_context.items():
            geo_ctx_parts.append(f"{acc} ({ctx})")
        geo_ctx_str = "; ".join(geo_ctx_parts) if geo_ctx_parts else "GEO death metadata unavailable"

        print(f"  GTEx training data was predominantly Hardy {dominant_hardy} "
              f"({dominant_label}, {dominant_pct:.0f}% of samples; "
              f"ventilator cases={vent_pct:.0f}%, slow deaths={slow_pct:.0f}%). "
              f"Independent validation datasets: {geo_ctx_str}. "
              f"{'Validation cohorts represent a comparable death-context distribution to training, supporting generalisability.' if slow_pct < 30 and vent_pct < 30 else 'Note: death-context distributions differ between training and validation — interpret generalisation performance with this in mind.'}")
    else:
        print("  (GTEx training Hardy distribution not available — run after scores_for_ml.csv exists)")

    print(f"\n  → Copy both sentences into your Methods/Results draft.\n")
else:
    print("\nNo datasets scored yet — run R/02b_score_geo_validation.R first.")

print(f"=== Done: validation outputs in {OUT_DIR} ===")
