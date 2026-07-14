# INDELVAR pipeline

Install the dependencies (`pip install -r requirements.txt`; `renv::restore()` in R) and set `INDELVAR_ROOT` to the repository root (scripts resolve all paths from it). The trained model and the ClinVar training set are in the repository. The genome-wide score database is on [Zenodo](https://doi.org/10.5281/zenodo.21285600); download it into `data/`.

Pipeline order: cohorts, features, train, calibrate, database, then score / look up.

## 1. Cohorts: `build_cohort.py`

MANE Select in-frame-indel cohorts (`train`, `test`, `vus`, or `benign`), built with Ensembl VEP. The `train` cohort applies the ClinVar review-status (>=1 star) inclusion gate.

- In: [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) and [gnomAD v4](https://gnomad.broadinstitute.org/downloads) VCFs, [Ensembl VEP](https://github.com/Ensembl/ensembl-vep)
- Out: per-cohort variant table

## 2. Features: `compute_features.py`

The 39 model features for a cohort, computed from per-protein and per-gene reference data.

- In: a cohort table, plus [AlphaFold](https://alphafold.ebi.ac.uk/) structures, [UniProt](https://www.uniprot.org/) 2024_06, [gnomAD v4.1 constraint](https://gnomad.broadinstitute.org/data#v4-constraint), [phyloP100way](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/) and [GERP](https://ftp.ensembl.org/pub/release-111/compara/conservation_scores/) bigWigs, and the [MANE Select](https://www.ncbi.nlm.nih.gov/refseq/MANE/) v1.5 GFF
- Out: feature table (tsv/parquet)

## 3. Train: `train.R`

Random forest (caret, ntree=1000, 5×3 CV); score = raw RF pathogenic probability.

- In: `dataset/train_data/train_set_features.tsv`
- Out: `model/{indelvar.rds, train_oof_predictions.tsv, feature_medians.rds}`

## 4. Calibrate: `calibrate.R` (+ `calib_params.R`)

Per-type ACMG cutoffs from the OOF predictions and the gnomAD benign background.

- In: `dataset/train_data/{train_oof_predictions,train_set_coords}.tsv`, `dataset/calibration/gnomad_calibration_scored.tsv`
- Out: `model/{evidence_cutoffs.tsv, evidence_strength_targets.tsv}`

## 5. Database: `build_database.py`

Assembles the genome-wide Zenodo tables (features, scores, or both) from the per-type chunks produced in stage 2.

- In: `data/precomputed/{features_*,precomputed_*,*_genomic}.parquet`
- Out: `data/precomputed/{release_features.parquet, release_scores.parquet, release_scores.tsv.gz}`

## 6. Score / look up

- `score_variants.R <features.{tsv,parquet}> <out.tsv> [del|ins]`: score a feature table (score, ACMG tier and points, SHAP groups).
- `indelvar_lookup.py --gene CFTR --hgvsp p.Leu127dup [--features]`: query the database by MANE gene and HGVS p. (reads `data/release_scores.parquet` from Zenodo).

## How to run

Run stages 1 to 5 in order, providing the reference data listed in stages 1 and 2, to rebuild the model and the precomputed database from scratch. To score your own variants without rebuilding, see Usage in the [README](../README.md).
