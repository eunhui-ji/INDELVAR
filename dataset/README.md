# Datasets

All datasets use GRCh38 and MANE Select v1.5 and contain in-frame insertions and deletions spanning 1–10 amino acids.

| Directory | Description |
|---|---|
| `training_set/` | ClinVar-derived model development set (n = 4,927) |
| `test_set/` | Independent ClinVar test set (n = 1,262) |
| `vus_set/` | ClinVar variants of uncertain significance (n = 20,639) |
| `calibration_reference/` | Population reference used for calibration (n = 23,253) |

Each directory contains variant annotations, model features, and/or INDELVAR predictions as applicable. The 37 model features are defined in [`../model/features.tsv`](../model/features.tsv). Missing values are written as `NA`.

