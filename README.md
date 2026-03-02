# HIForensics

**A cell-type-resolved hypoxic transcriptional decay atlas for postmortem interval estimation**

---

## Biological motivation

When a person dies, tissues become hypoxic and the HIF (Hypoxia-Inducible Factor) pathway activates within minutes. As ischemic time increases, HIF target gene expression rises, plateaus, and eventually decays as cells lyse and RNA degrades. This creates a transcriptional clock that begins ticking at the moment of death.

The GTEx project has collected RNA-seq from 50+ tissues across ~900 donors, each with a recorded **ischemic time** (time from death to tissue preservation). Nobody has systematically characterised how HIF pathway activation and subsequent decay varies by tissue and cell type across that ischemic window — that is the gap this project addresses.

**Key scientific contributions:**
- First systematic multi-tissue, cell-type-aware HIF-based PMI estimator
- Trained on the largest postmortem expression dataset in existence (GTEx V8, ~17,000 samples)
- Mechanistically grounded: HIF biology explains *why* the predictor works, not just that it does
- Validated against three independent GEO postmortem datasets

---

## Pipeline architecture

```
Phase 1 — Scoring
  scripts/00_download_gtex.sh          GTEx V8 expression + metadata (~14 GB)
  scripts/01_download_msigdb.R         HIF/hypoxia gene sets via msigdbr
  R/02_score_hif.R                     GSVA scoring across GTEx tissues
                                        Partial Spearman vs ischemic time
                                        (Hardy death scale controlled)

Phase 2 — Machine learning
  python/03_train_pmi_model.py         XGBoost PMI predictor
                                        Tissue-specific models (HIF only)
                                        Global model (HIF + tissue)
                                        Donor-grouped cross-validation
                                        SHAP feature importance

Phase 3 — Validation
  scripts/02_download_geo_validation.sh  GSE45642 · GSE53987 · GSE116754
  R/02b_score_geo_validation.R           GSVA scoring of GEO datasets
  python/04_validate.py                  Pure inference, no retraining
                                          Hardy distribution comparison

Figures
  R/05_figures.R                       Publication-ready figures from CSVs
```

**Shared config:** `R/gsva_config.R` — gene set path, GSVA parameters, expression
filters. Edit here only; both scoring scripts source it automatically.

---

## Reproducing the pipeline

### Dependencies

**R (≥ 4.3)**
```r
install.packages(c("data.table", "dplyr", "tidyr", "readr",
                   "ggplot2", "cowplot", "pheatmap", "ggrepel",
                   "scales", "viridis", "patchwork", "RColorBrewer"))

if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("GSVA", "BiocParallel", "GEOquery",
                       "scDblFinder", "limma"))

install.packages("msigdbr")
```

**Python (≥ 3.10)**
```bash
pip install -r python/requirements.txt
# pandas==2.3.3  numpy==2.4.0  scipy==1.17.0  scikit-learn==1.8.0
# xgboost==3.2.0  joblib==1.5.3  shap==0.50.0
```

### Step-by-step reproduction

```bash
# 1. Gene sets — fast, run first
Rscript scripts/01_download_msigdb.R

# 2. GTEx V8 — ~14 GB, run overnight
bash scripts/00_download_gtex.sh

# 3. HIF scoring across GTEx (~30–60 min, 20 cores, ~25 GB RAM)
Rscript R/02_score_hif.R

# 4. Train PMI model (~2–3 hr, 20 cores)
python3 python/03_train_pmi_model.py

# 5. Download and score GEO validation datasets
bash scripts/02_download_geo_validation.sh
Rscript R/02b_score_geo_validation.R

# 6. Validate (pure inference)
python3 python/04_validate.py

# 7. Generate figures
Rscript R/05_figures.R
```

### Dependency order

```
01_download_msigdb.R
        │
        ▼
00_download_gtex.sh ──► 02_score_hif.R ──► 03_train_pmi_model.py
                                                      │
02_download_geo_validation.sh                         │
        │                                             ▼
        ▼                               02b_score_geo_validation.R
02b_score_geo_validation.R  ──────────► 04_validate.py

All results CSVs ──► 05_figures.R  (can be written before pipeline finishes)
```

---

## Data sources

| Dataset | Description | Access |
|---|---|---|
| [GTEx V8](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression) | ~17,000 RNA-seq samples, 54 tissues, ischemic time recorded | Open — no dbGaP required |
| [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/) | HIF/hypoxia gene sets (Hallmark, Reactome, GO, KEGG) | Open via `msigdbr` |
| [GSE45642](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45642) | Postmortem brain, multiple PMI timepoints | Open |
| [GSE53987](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53987) | Postmortem brain, well-characterised PMI | Open |
| [GSE116754](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116754) | Postmortem RNA stability, multi-tissue | Open |

All data is fully open access. No controlled-access application required.

---

## Key outputs

| File | Description |
|---|---|
| `results/hif_scores/tissue_correlations.csv` | Partial Spearman r (HIF vs PMI) per tissue × gene set, Hardy-controlled |
| `results/hif_scores/summary.csv` | Tissues ranked by mean \|r\| — brain regions should rank first |
| `results/hif_scores/scores_for_ml.csv` | Flat table of GSVA scores + metadata for ML |
| `results/ml_model/global_model.joblib` | Deployable XGBoost PMI predictor |
| `results/ml_model/tissue_models.joblib` | Per-tissue models (HIF scores only) |
| `results/ml_model/shap_importance.csv` | SHAP feature group breakdown (HIF vs tissue vs Hardy) |
| `results/ml_model/model_comparison.csv` | ΔR² — HIF contribution beyond tissue identity alone |
| `results/validation/validation_summary.csv` | Generalisation to independent GEO datasets |

---

## Design decisions worth noting

**Why GTEx V8, not V10?** V8 is fully open access with no dbGaP requirement. Full reproducibility — anyone can run this pipeline.

**Why partial Spearman, not Pearson?** HIF-PMI relationships are nonlinear (activation plateau followed by decay). Rank-based correlation is more robust.

**Why Hardy scale is controlled but not excluded from ML features?** In the correlation analysis (Phase 1), Hardy scale is a confounder — we partial it out. In the ML model (Phase 2), it is included explicitly as a covariate so the model accounts for it rather than learning it implicitly through the HIF features. A prospective-use model without Hardy scale can be derived from the same training run.

**Why donor-grouped cross-validation?** The same donor contributes samples from multiple tissues. Random splitting leaks the shared ischemic time across folds and inflates performance estimates.

---

## Citation

Manuscript in preparation.
