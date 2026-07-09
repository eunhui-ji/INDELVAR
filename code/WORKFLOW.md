# INDELVAR pipeline

Install the dependencies (`pip install -r requirements.txt`; `renv::restore()` in R) and set `INDELVAR_ROOT` to the repository root (scripts resolve all paths from it). The trained model and the ClinVar training set are in the repository. The genome-wide score database is on [Zenodo](https://doi.org/10.5281/zenodo.21285601); download it into `data/`. The held-out benchmark test set is under HGMD license and is supplied by the user.

Pipeline order: cohorts, features, train, calibrate, database, then score/lookup; benchmark, tables, and figures are downstream.

## 1. Cohorts: `build_cohort.py`

MANE Select in-frame-indel cohorts (`train`, `test`, `vus`, or `benign`), built with Ensembl VEP.

- In: [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) and [gnomAD v4](https://gnomad.broadinstitute.org/downloads) VCFs, [Ensembl VEP](https://github.com/Ensembl/ensembl-vep)
- Out: per-cohort variant table

## 2. Features: `compute_features.py`

The 39 model features for a cohort, computed from per-protein and per-gene reference data.

- In: a cohort table, plus [AlphaFold](https://alphafold.ebi.ac.uk/) structures, [UniProt](https://www.uniprot.org/) 2024_06, [gnomAD v4.1 constraint](https://gnomad.broadinstitute.org/data#v4-constraint), [phyloP100way](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/) and [GERP](https://ftp.ensembl.org/pub/release-111/compara/conservation_scores/) bigWigs, and the [MANE Select](https://www.ncbi.nlm.nih.gov/refseq/MANE/) v1.5 GFF
- Out: feature table (tsv/parquet)

## 3. Train: `train.R`

Random forest (caret, ntree=1000, 5×3 CV); score = raw RF pathogenic probability.

- In: `dataset/train_data/train_set.tsv`
- Out: `model/{indelvar.rds, train_oof_predictions.tsv, feature_medians.rds}`

## 4. Calibrate: `calibrate.R` (+ `calib_params.R`)

Per-type ACMG cutoffs from the OOF predictions and the gnomAD benign background.

- In: `dataset/train_data/{train_oof_predictions,train_set_coords}.tsv`, `dataset/calibration/gnomad_benign_scored.tsv`
- Out: `model/{evidence_cutoffs.tsv, cutoffs_local_1sided95.tsv}`

## 5. Database: `build_database.py`

Assembles the genome-wide Zenodo tables (features, scores, or both) from the per-type chunks produced in stage 2.

- In: `data/precomputed/{features_*,precomputed_*,*_genomic}.parquet`
- Out: `data/precomputed/{release_features.parquet, release_scores.parquet, release_scores.tsv.gz}`

## 6. Score / look up

- `score_variants.R <features.{tsv,parquet}> <out.tsv> [del|ins]`: score a feature table (score, ACMG tier and points, SHAP groups).
- `indelvar_lookup.py --gene CFTR --hgvsp p.Leu127dup [--features]`: query the database by MANE gene and HGVS p. (reads `data/release_scores.parquet` from Zenodo).

## 7. Benchmark: `benchmark.R`

INDELVAR against external tools on the held-out test set, which is under HGMD license and user-supplied (place it under `dataset/test_data/`). Feeds Tables S4 and S5.

## 8. Tables (`code/tables/`)

- `gen_tables.R`: Tables 1, S4, S5, S7.
- `gen_tables.py`: Table S6.

## 9. Figures (`code/figures/`)

| script | figure |
|---|---|
| `plot_landscape.R` | Figure 2 |
| `plot_features.R` | Figure 3 |
| `plot_importance.R` | Figure 4A |
| `plot_benchmark_roc.R` | Figure 4B,C * |
| `plot_calibration.R` | Figure 5 |
| `plot_application.R` | Figure 6 |
| `plot_hyperparam.R` | Figure S1 |
| `plot_feature_importance.R` | Figure S2 |
| `plot_feature_correlation.R` | Figure S3 |
| `plot_pr_curves.R` | Figure S4 * |
| `plot_subgroup_auroc.R` | Figure S5 |
| `plot_calibration_decision.R` | Figure S6 * |

`*` uses the HGMD-licensed test set (step 7). Figure 1 is not scripted.

## How to run

- **Reproduce a figure or table**: run any `code/figures/plot_*.R` or `code/tables/gen_tables.{R,py}`; the figures marked `*` and the test-set tables (S4, S5, S6) use the HGMD data.
- **Rebuild from scratch**: run stages 1 to 5 in order, providing the reference data listed in stages 1 and 2.

To score your own variants, see Usage in the [README](../README.md).
