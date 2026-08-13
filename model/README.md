# Model files

INDELVAR uses a random forest trained on 4,927 in-frame indels with 37 features. The model returns a score from 0 to 1.

| Setting | Value |
|---|---|
| Algorithm | Random forest |
| Implementation | R `randomForest` 4.7.1.2 |
| Trees | 1,000 |
| Variables per split (`mtry`) | 7 |
| Terminal node size | 1 |
| Sampling fraction | 1.0 with replacement |
| Maximum nodes | Unlimited |

Missing features are replaced with training-set medians during scoring. PP3/BP4 evidence thresholds are calibrated separately for deletions and insertions.

## Files

| File | Description |
|---|---|
| `indelvar_rf.rds` | Trained random-forest model used to calculate INDELVAR scores |
| `features.tsv` | The 37 model features in order: what each one measures, its unit, when it is missing, and where it comes from |
| `imputation_medians.tsv` | Training-set median used for each feature when it is missing |
| `evidence_cutoffs.tsv` | Score thresholds at which each PP3/BP4 evidence strength is reached |
| `evidence_strength_targets.tsv` | Likelihood-ratio targets those thresholds were calibrated to meet |
| `model_manifest.json` | Model version, hyperparameters, and SHA-256 checksums |
