"""
03_train_pmi_model.py
HIForensics Phase 2 — PMI prediction from HIF gene set scores

Design principles
-----------------
1. Log1p target: ischemic_min is right-skewed (long tail past 1000 min).
   Predict log1p(ischemic_min), evaluate predictions on original-scale minutes.

2. Two model types, different questions:
   a. Tissue-specific XGBoost (primary analysis):
      Features = HIF scores only. Trained within each tissue separately.
      Answers: "Do HIF scores predict PMI independent of tissue identity?"
   b. Global XGBoost (deployable model):
      Features = HIF scores + tissue (one-hot). Works on any input tissue.
      A tissue-only baseline is also trained to confirm HIF adds signal
      beyond just knowing the tissue.

3. Donor-grouped CV throughout (GroupKFold, group = donor_id).
   A donor's samples appear in exactly one fold. This is critical because
   the same donor contributes samples from multiple tissues — random splitting
   leaks the shared ischemic time across folds and inflates performance.

4. Held-out test set: 20% of donors reserved before any fitting.
   Stratified to ensure all tissues are represented in both splits.

Input:  results/hif_scores/scores_for_ml.csv
Output: results/ml_model/
"""

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import numpy as np
import pandas as pd
from scipy import stats
import joblib
import json
import os
from pathlib import Path

from sklearn.model_selection import GroupKFold, GroupShuffleSplit
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
import xgboost as xgb
import shap

# =============================================================================
# CONFIG
# =============================================================================
PROJECT  = Path("/media/ross/New Volume1/HIForensics")
IN_FILE  = PROJECT / "results" / "hif_scores" / "scores_for_ml.csv"
OUT_DIR  = PROJECT / "results" / "ml_model"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_FOLDS       = 5      # GroupKFold splits
TEST_FRAC     = 0.20   # donor fraction held out as final test set
MIN_TISSUE_N  = 30     # minimum samples for a tissue-specific model
RANDOM_STATE  = 42
N_JOBS        = 20

# XGBoost hyperparameters — reasonable defaults for tabular biomedical data
# Tune with Optuna in a later pass if needed
XGB_PARAMS = dict(
    n_estimators      = 500,
    max_depth         = 4,       # shallow trees reduce overfitting on small N
    learning_rate     = 0.05,
    subsample         = 0.8,
    colsample_bytree  = 0.8,
    min_child_weight  = 5,
    reg_alpha         = 0.1,
    reg_lambda        = 1.0,
    tree_method       = "hist",  # "hist" + device="cuda" is correct for XGBoost >=2.0
    device            = "cuda",  # RTX 4070 Ti — gpu_hist is deprecated in XGBoost 3.x
    random_state      = RANDOM_STATE,
    n_jobs            = N_JOBS,
    verbosity         = 0,
)

# =============================================================================
# HELPERS
# =============================================================================

def eval_metrics(y_true_log, y_pred_log, label=""):
    """Evaluate on original scale (ischemic minutes) after back-transforming."""
    y_true = np.expm1(y_true_log)
    y_pred = np.expm1(y_pred_log)
    y_pred = np.clip(y_pred, 0, None)   # predictions can't be negative minutes

    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    mae  = mean_absolute_error(y_true, y_pred)
    r2   = r2_score(y_true, y_pred)
    rho, p = stats.spearmanr(y_true, y_pred)

    if label:
        print(f"  {label}: RMSE={rmse:.1f} min  MAE={mae:.1f} min  "
              f"R²={r2:.3f}  Spearman r={rho:.3f} (p={p:.2e})")
    return dict(rmse=rmse, mae=mae, r2=r2, spearman_r=rho, spearman_p=p)


def age_to_midpoint(age_str):
    """Convert GTEx age range string (e.g. '50-59') to numeric midpoint."""
    try:
        lo, hi = map(int, str(age_str).split("-"))
        return (lo + hi) / 2
    except Exception:
        return np.nan


def donor_stratified_split(df, test_frac, random_state):
    """
    Split donors into train/test ensuring every tissue appears in both sets.
    Falls back to pure GroupShuffleSplit if a tissue has too few donors to
    guarantee representation.
    """
    donors = df[["donor_id", "tissue"]].drop_duplicates()
    # Tissues with only one donor can't be stratified — flag them
    single_donor_tissues = (
        donors.groupby("tissue")["donor_id"].nunique()
        .pipe(lambda s: s[s == 1].index.tolist())
    )
    if single_donor_tissues:
        print(f"  Note: {len(single_donor_tissues)} tissue(s) have only 1 donor "
              f"— cannot guarantee test representation for these.")

    gss = GroupShuffleSplit(n_splits=1, test_size=test_frac,
                            random_state=random_state)
    idx_train, idx_test = next(gss.split(df, groups=df["donor_id"]))
    return df.iloc[idx_train].copy(), df.iloc[idx_test].copy()


# =============================================================================
# 1. LOAD AND VALIDATE
# =============================================================================
print(f"=== HIForensics PMI Model Training  ===")
print(f"Input: {IN_FILE}\n")

if not IN_FILE.exists():
    raise FileNotFoundError(
        f"{IN_FILE} not found.\n"
        "Run R/02_score_hif.R first to generate scores_for_ml.csv."
    )

df = pd.read_csv(IN_FILE)
print(f"Loaded: {len(df):,} samples × {df.shape[1]} columns")

# Identify HIF gene set score columns (everything that isn't a metadata column)
META_COLS = ["sample_id", "donor_id", "tissue", "ischemic_min",
             "hardy_scale", "sex", "age_range"]
HIF_COLS  = [c for c in df.columns if c not in META_COLS]
print(f"HIF gene set columns ({len(HIF_COLS)}): {HIF_COLS}\n")

# --- Clean and engineer features
df["age_mid"]     = df["age_range"].apply(age_to_midpoint)
df["target_log"]  = np.log1p(df["ischemic_min"])

# Drop samples with missing target
n_before = len(df)
df = df.dropna(subset=["ischemic_min", "donor_id", "tissue"])
print(f"After dropping missing target/donor/tissue: {len(df):,} "
      f"(removed {n_before - len(df)})")

# Ischemic time distribution
print(f"\nIschemic time (minutes):")
print(f"  range:  {df.ischemic_min.min():.0f} – {df.ischemic_min.max():.0f}")
print(f"  median: {df.ischemic_min.median():.0f}  mean: {df.ischemic_min.mean():.0f}")
print(f"  >1000 min: {(df.ischemic_min > 1000).sum()} samples "
      f"({100*(df.ischemic_min>1000).mean():.1f}% — long-tail confirmed, log transform correct)")
print(f"\nTissues: {df.tissue.nunique()}")
print(f"Donors:  {df.donor_id.nunique()}\n")

# =============================================================================
# 2. HELD-OUT TEST SET (20% of donors, set aside before any fitting)
# =============================================================================
print("--- Splitting held-out test set (20% of donors) ---")
df_train, df_test = donor_stratified_split(df, TEST_FRAC, RANDOM_STATE)

train_donors = df_train.donor_id.nunique()
test_donors  = df_test.donor_id.nunique()
print(f"  Train: {len(df_train):,} samples ({train_donors} donors)")
print(f"  Test:  {len(df_test):,} samples ({test_donors} donors)")
print(f"  Test tissues covered: {df_test.tissue.nunique()} / {df.tissue.nunique()}\n")

# =============================================================================
# 3. TISSUE-SPECIFIC MODELS (primary analysis)
#    Features: HIF scores only. No tissue feature — answers whether HIF alone
#    predicts PMI within a tissue, removing the tissue-identity confounder.
# =============================================================================
print("=" * 70)
print("ANALYSIS 1: Tissue-specific XGBoost models (HIF scores only)")
print("=" * 70)

tissue_cv_records = []
tissue_test_records = []
tissue_models = {}

tissues = sorted(df_train.tissue.unique())
tissues_eligible = [t for t in tissues
                    if (df_train.tissue == t).sum() >= MIN_TISSUE_N]
print(f"Tissues with ≥{MIN_TISSUE_N} train samples: {len(tissues_eligible)} / {len(tissues)}\n")

for tissue in tissues_eligible:
    tr = df_train[df_train.tissue == tissue].copy()
    te = df_test[df_test.tissue == tissue].copy() if tissue in df_test.tissue.values else None

    X_tr = tr[HIF_COLS].values
    y_tr = tr["target_log"].values
    groups_tr = tr["donor_id"].values

    n_donors = tr.donor_id.nunique()
    n_folds  = min(N_FOLDS, n_donors)

    if n_donors < 2:
        continue

    # --- Cross-validation
    gkf  = GroupKFold(n_splits=n_folds)
    oof_pred  = np.full(len(tr), np.nan)
    fold_metrics = []

    for fold, (idx_fit, idx_val) in enumerate(gkf.split(X_tr, y_tr, groups_tr)):
        model = xgb.XGBRegressor(**XGB_PARAMS)
        model.fit(X_tr[idx_fit], y_tr[idx_fit],
                  eval_set=[(X_tr[idx_val], y_tr[idx_val])],
                  verbose=False)
        oof_pred[idx_val] = model.predict(X_tr[idx_val])
        m = eval_metrics(y_tr[idx_val], oof_pred[idx_val])
        m["fold"] = fold
        fold_metrics.append(m)

    oof_m = eval_metrics(y_tr, oof_pred)
    tissue_cv_records.append(dict(tissue=tissue, n_train=len(tr),
                                  n_donors=n_donors, **oof_m))

    # --- Train final tissue model on all training samples
    final_model = xgb.XGBRegressor(**XGB_PARAMS)
    final_model.fit(X_tr, y_tr, verbose=False)
    tissue_models[tissue] = final_model

    # --- Test set evaluation
    if te is not None and len(te) >= 5:
        X_te = te[HIF_COLS].values
        y_te = te["target_log"].values
        y_pred_te = final_model.predict(X_te)
        te_m = eval_metrics(y_te, y_pred_te)
        tissue_test_records.append(dict(tissue=tissue, n_test=len(te), **te_m))

# Print tissue-specific CV summary
ts_cv_df = pd.DataFrame(tissue_cv_records).sort_values("spearman_r", ascending=False)
print(f"\nTissue-specific CV results (top 20 by Spearman r):")
print(f"  {'Tissue':<45}  {'n':>5}  {'r':>6}  {'RMSE':>7}  {'MAE':>7}")
print("  " + "-" * 72)
for _, row in ts_cv_df.head(20).iterrows():
    print(f"  {row.tissue:<45}  {int(row.n_train):>5}  "
          f"{row.spearman_r:>6.3f}  {row.rmse:>7.1f}  {row.mae:>7.1f}")

ts_cv_df.to_csv(OUT_DIR / "tissue_specific_cv.csv", index=False)
pd.DataFrame(tissue_test_records).to_csv(OUT_DIR / "tissue_specific_test.csv", index=False)
print(f"\nSaved: tissue_specific_cv.csv, tissue_specific_test.csv")

# Save all tissue models
joblib.dump(tissue_models, OUT_DIR / "tissue_models.joblib")
print(f"Saved: tissue_models.joblib ({len(tissue_models)} tissue models)\n")

# =============================================================================
# 4. GLOBAL MODEL — tissue as feature (deployable one-stop-shop)
#    Includes: HIF scores + tissue (one-hot) + sex + age_mid
#    Hardy scale intentionally EXCLUDED: not available in prospective forensics
# =============================================================================
print("=" * 70)
print("ANALYSIS 2: Global XGBoost model (HIF scores + tissue as feature)")
print("=" * 70)
print("Note: Hardy scale excluded — not available prospectively in casework\n")

# hardy_scale included as a covariate: slow deaths (Hardy 4) have systematically
# different HIF signatures and without it the model risks learning Hardy scale
# implicitly through the HIF features. Include it explicitly so the model can
# account for it cleanly. hardy_scale can be NA for some GTEx donors, so it
# goes through SimpleImputer(median) before StandardScaler.
# Note: a prospective-only variant without hardy_scale can be derived later.
GLOBAL_FEAT_COLS = HIF_COLS + ["sex", "age_mid", "hardy_scale"]
CAT_COLS         = ["tissue"]
NUM_COLS         = GLOBAL_FEAT_COLS

num_transformer = Pipeline([
    ("impute", SimpleImputer(strategy="median")),  # handles NA in hardy_scale/age
    ("scale",  StandardScaler()),
])

preprocessor = ColumnTransformer(transformers=[
    ("num", num_transformer,                                              NUM_COLS),
    ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), CAT_COLS),
], remainder="drop")

def build_global_pipeline(model):
    return Pipeline([("pre", preprocessor), ("model", model)])

# Prepare full feature frame
def make_X(data):
    return data[NUM_COLS + CAT_COLS].copy()

X_tr_full = make_X(df_train)
y_tr_full  = df_train["target_log"].values
g_tr_full  = df_train["donor_id"].values

X_te_full  = make_X(df_test)
y_te_full   = df_test["target_log"].values

# --- 4a. Tissue-only baseline (to quantify how much HIF adds)
print("--- Baseline: tissue-only model ---")
base_preprocessor = ColumnTransformer(transformers=[
    ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), ["tissue"]),
], remainder="drop")
baseline_pipe = Pipeline([
    ("pre", base_preprocessor),
    ("model", xgb.XGBRegressor(**XGB_PARAMS))
])

gkf = GroupKFold(n_splits=N_FOLDS)
base_oof = np.full(len(df_train), np.nan)
for idx_fit, idx_val in gkf.split(X_tr_full, y_tr_full, g_tr_full):
    baseline_pipe.fit(
        X_tr_full.iloc[idx_fit][["tissue"]],
        y_tr_full[idx_fit]
    )
    base_oof[idx_val] = baseline_pipe.predict(X_tr_full.iloc[idx_val][["tissue"]])

base_cv_m = eval_metrics(y_tr_full, base_oof, label="Tissue-only CV")

# Retrain on all train data, evaluate on test
baseline_pipe.fit(df_train[["tissue"]], y_tr_full)
base_test_pred = baseline_pipe.predict(df_test[["tissue"]])
base_test_m = eval_metrics(y_te_full, base_test_pred, label="Tissue-only TEST")

# --- 4b. Global HIF + tissue model
print("\n--- Global model: HIF scores + tissue ---")
global_pipe = build_global_pipeline(xgb.XGBRegressor(**XGB_PARAMS))

global_oof = np.full(len(df_train), np.nan)
for idx_fit, idx_val in gkf.split(X_tr_full, y_tr_full, g_tr_full):
    global_pipe.fit(X_tr_full.iloc[idx_fit], y_tr_full[idx_fit])
    global_oof[idx_val] = global_pipe.predict(X_tr_full.iloc[idx_val])

global_cv_m = eval_metrics(y_tr_full, global_oof, label="Global (HIF+tissue) CV")

# Retrain on all train data, evaluate on test
global_pipe.fit(X_tr_full, y_tr_full)
global_test_pred = global_pipe.predict(X_te_full)
global_test_m = eval_metrics(y_te_full, global_test_pred, label="Global (HIF+tissue) TEST")

# Delta R² — how much HIF adds on top of tissue alone
delta_r2_cv   = global_cv_m["r2"]   - base_cv_m["r2"]
delta_r2_test = global_test_m["r2"] - base_test_m["r2"]
print(f"\n  ΔR² (HIF contribution beyond tissue identity): "
      f"CV={delta_r2_cv:+.3f}  Test={delta_r2_test:+.3f}")

# --- 4c. Random Forest (comparison)
print("\n--- Comparison: Random Forest ---")
rf_pipe = build_global_pipeline(
    RandomForestRegressor(n_estimators=300, max_depth=6,
                          min_samples_leaf=5, n_jobs=N_JOBS,
                          random_state=RANDOM_STATE)
)
rf_oof = np.full(len(df_train), np.nan)
for idx_fit, idx_val in gkf.split(X_tr_full, y_tr_full, g_tr_full):
    rf_pipe.fit(X_tr_full.iloc[idx_fit], y_tr_full[idx_fit])
    rf_oof[idx_val] = rf_pipe.predict(X_tr_full.iloc[idx_val])
eval_metrics(y_tr_full, rf_oof, label="Random Forest CV")

# =============================================================================
# 5. SAVE GLOBAL MODEL + OOF PREDICTIONS
# =============================================================================
print("\n--- Saving global model and predictions ---")
joblib.dump(global_pipe, OUT_DIR / "global_model.joblib")
print("Saved: global_model.joblib")

# Out-of-fold predictions frame
oof_df = df_train[["sample_id", "donor_id", "tissue", "ischemic_min",
                    "hardy_scale"]].copy()
oof_df["predicted_min"] = np.clip(np.expm1(global_oof), 0, None)
oof_df["residual_min"]  = oof_df["predicted_min"] - oof_df["ischemic_min"]
oof_df.to_csv(OUT_DIR / "oof_predictions.csv", index=False)
print("Saved: oof_predictions.csv")

# Test set predictions
test_df = df_test[["sample_id", "donor_id", "tissue", "ischemic_min",
                    "hardy_scale"]].copy()
test_df["predicted_min"] = np.clip(np.expm1(global_test_pred), 0, None)
test_df["residual_min"]  = test_df["predicted_min"] - test_df["ischemic_min"]
test_df.to_csv(OUT_DIR / "test_predictions.csv", index=False)
print("Saved: test_predictions.csv")

# Model comparison summary
comparison = pd.DataFrame([
    dict(model="tissue_only",        split="cv",   **base_cv_m),
    dict(model="tissue_only",        split="test", **base_test_m),
    dict(model="global_hif+tissue",  split="cv",   **global_cv_m),
    dict(model="global_hif+tissue",  split="test", **global_test_m),
])
comparison.to_csv(OUT_DIR / "model_comparison.csv", index=False)
print("Saved: model_comparison.csv")

# =============================================================================
# 6. SHAP FEATURE IMPORTANCE (global model)
# =============================================================================
print("\n--- Computing SHAP values (global model) ---")

# SHAP on XGBoost directly (faster than through Pipeline)
# Get the transformed feature matrix for the training set
X_tr_transformed = global_pipe["pre"].transform(X_tr_full)
feature_names = (
    NUM_COLS +
    list(global_pipe["pre"].named_transformers_["cat"]
         .get_feature_names_out(CAT_COLS))
)

explainer    = shap.TreeExplainer(global_pipe["model"])
shap_values  = explainer.shap_values(X_tr_transformed)
mean_abs_shap = np.abs(shap_values).mean(axis=0)

shap_df = (
    pd.DataFrame({"feature": feature_names, "mean_abs_shap": mean_abs_shap})
    .sort_values("mean_abs_shap", ascending=False)
    .reset_index(drop=True)
)

# Summarise: HIF sets vs tissue vs Hardy scale vs demographics
# Hardy scale is tracked separately — it's a confounder, not a scientific signal.
# If its SHAP is large that confirms it was right to include it explicitly.
hif_shap    = shap_df[shap_df.feature.isin(HIF_COLS)].mean_abs_shap.sum()
tissue_shap = shap_df[shap_df.feature.str.startswith("tissue_")].mean_abs_shap.sum()
hardy_shap  = shap_df[shap_df.feature == "hardy_scale"].mean_abs_shap.sum()
demo_shap   = shap_df[shap_df.feature.isin(["sex", "age_mid"])].mean_abs_shap.sum()
total_shap  = hif_shap + tissue_shap + hardy_shap + demo_shap

print(f"\n  Feature group contributions (mean |SHAP|, summed):")
print(f"    HIF gene sets:       {hif_shap:.4f}  ({100*hif_shap/total_shap:.1f}%)")
print(f"    Tissue identity:     {tissue_shap:.4f}  ({100*tissue_shap/total_shap:.1f}%)")
print(f"    Hardy death scale:   {hardy_shap:.4f}  ({100*hardy_shap/total_shap:.1f}%)")
print(f"    Demographics:        {demo_shap:.4f}  ({100*demo_shap/total_shap:.1f}%)")
if hardy_shap / total_shap > 0.15:
    print(f"    → Hardy scale is a major predictor: good call including it explicitly.")
print(f"\n  Top 10 individual features:")
for _, row in shap_df.head(10).iterrows():
    print(f"    {row.feature:<50}  {row.mean_abs_shap:.4f}")

shap_df.to_csv(OUT_DIR / "shap_importance.csv", index=False)
print("\nSaved: shap_importance.csv")

# =============================================================================
# 7. PER-TISSUE PERFORMANCE (global model, OOF predictions)
# =============================================================================
print("\n--- Per-tissue performance (global model, OOF) ---")
per_tissue = []
for tissue, grp in oof_df.groupby("tissue"):
    if len(grp) < 10:
        continue
    m = eval_metrics(
        np.log1p(grp.ischemic_min.values),
        np.log1p(grp.predicted_min.values)
    )
    per_tissue.append(dict(tissue=tissue, n=len(grp), **m))

per_tissue_df = (pd.DataFrame(per_tissue)
                 .sort_values("spearman_r", ascending=False))
per_tissue_df.to_csv(OUT_DIR / "per_tissue_performance.csv", index=False)
print("Saved: per_tissue_performance.csv")

# Brain validation check
brain_mask = per_tissue_df.tissue.str.contains(
    "Brain|Cortex|Cerebellum|Hippocampus|Caudate|Putamen|Nucleus accumbens|"
    "Amygdala|Substantia|Frontal|Hypothalamus|Spinal", case=False, regex=True
)
brain_df   = per_tissue_df[brain_mask]
brain_ranks = per_tissue_df.index[brain_mask].tolist()
overall_n   = len(per_tissue_df)

print(f"\n  Brain region performance (global model):")
print(f"  {'Tissue':<45}  {'rank':>5}  {'r':>6}  {'RMSE':>7}")
for i, (_, row) in enumerate(brain_df.iterrows()):
    rank = per_tissue_df.index.get_loc(row.name) + 1
    print(f"  {row.tissue:<45}  {rank:>5}  {row.spearman_r:>6.3f}  {row.rmse:>7.1f}")

if len(brain_df):
    med_rank = int(np.median([per_tissue_df.index.get_loc(i) + 1
                               for i in brain_df.index]))
    print(f"\n  Median brain rank: {med_rank} / {overall_n}  "
          f"({'PASS' if med_rank <= 10 else 'WARNING — investigate'})")

# =============================================================================
# 8. SAVE METADATA
# =============================================================================
meta = {
    "n_samples_total":  int(len(df)),
    "n_donors_total":   int(df.donor_id.nunique()),
    "n_tissues":        int(df.tissue.nunique()),
    "n_train_samples":  int(len(df_train)),
    "n_test_samples":   int(len(df_test)),
    "hif_cols":         HIF_COLS,
    "n_folds":          N_FOLDS,
    "xgb_params":       XGB_PARAMS,
    "cv_r2_global":         round(global_cv_m["r2"],             4),
    "cv_r2_baseline":       round(base_cv_m["r2"],               4),
    "delta_r2_cv":          round(delta_r2_cv,                   4),
    "test_r2_global":       round(global_test_m["r2"],           4),
    "delta_r2_test":        round(delta_r2_test,                 4),
    "test_rmse_min":        round(global_test_m["rmse"],         2),
    "test_mae_min":         round(global_test_m["mae"],          2),
    "test_spearman_r":      round(global_test_m["spearman_r"],   4),
    "shap_pct_hif":         round(100 * hif_shap / total_shap,   1),
    "shap_pct_tissue":      round(100 * tissue_shap / total_shap, 1),
    "shap_pct_hardy":       round(100 * hardy_shap / total_shap,  1),
    "shap_pct_demographics":round(100 * demo_shap / total_shap,   1),
}
with open(OUT_DIR / "run_metadata.json", "w") as f:
    json.dump(meta, f, indent=2, default=float)
print("\nSaved: run_metadata.json")

print(f"\n=== Done. All outputs in {OUT_DIR} ===")
print("Next step: python/04_validate.py (independent GEO datasets)")
