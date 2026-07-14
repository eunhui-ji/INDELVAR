#!/usr/bin/env Rscript
## Train the INDELVAR random forest on the ClinVar in-frame-indel training set.
suppressPackageStartupMessages({
  library(data.table); library(caret); library(randomForest); library(pROC)
})
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
SEED <- 42; set.seed(SEED)
SRC  <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- Sys.getenv("INDELVAR_MODEL_OUT", file.path(SRC, "model"))
dir.create(OUT, showWarnings = FALSE)
say  <- function(...) cat(format(Sys.time(), "%H:%M:%S "), sprintf(...), "\n", sep = "")

full  <- fread(file.path(SRC, "dataset", "train_data", "train_set_features.tsv"))
feats <- readRDS(file.path(OUT, "col_definitions.rds"))$feature_cols
say("training set = %d (P=%d B=%d; del=%d ins=%d); features = %d", nrow(full),
    sum(full$label == "Pathogenic"), sum(full$label == "Benign"),
    sum(full$indel_type == "deletion"), sum(full$indel_type == "insertion"), length(feats))

df <- as.data.frame(full)
df$label <- factor(as.character(df$label), levels = c("Benign", "Pathogenic"))
if ("indel_type" %in% feats)
  df$indel_type <- factor(as.character(df$indel_type), levels = c("deletion", "insertion"))
saveRDS(df, file.path(OUT, "train_set.rds"))

ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 3, classProbs = TRUE,
                     summaryFunction = twoClassSummary, savePredictions = "final")
p <- length(feats)
grid <- expand.grid(mtry = sort(unique(c(max(1, floor(sqrt(p))), max(1, floor(p / 3)),
                                         max(1, floor(p / 2)), p))))
say("mtry grid: %s", paste(grid$mtry, collapse = ", "))

t0 <- Sys.time()
rf_fit <- train(x = df[, feats], y = df$label, method = "rf", trControl = ctrl,
                tuneGrid = grid, ntree = 1000, importance = TRUE, metric = "ROC")
say("RF trained in %.1f min; bestTune mtry = %d", as.numeric(Sys.time() - t0, units = "mins"),
    rf_fit$bestTune$mtry)
saveRDS(rf_fit, file.path(OUT, "indelvar.rds"))

oof     <- as.data.table(rf_fit$pred)[mtry == rf_fit$bestTune$mtry]
oof_avg <- oof[, .(indelvar = mean(Pathogenic), label_oof = first(obs)), by = rowIndex]
oof_avg[, variant_id := df$variant_id[rowIndex]]
oof_avg[, y := as.integer(label_oof == "Pathogenic")]
fwrite(oof_avg, file.path(OUT, "train_oof_predictions.tsv"), sep = "\t")
say("OOF: N = %d (P = %d B = %d); AUROC(raw) = %.4f", nrow(oof_avg),
    sum(oof_avg$y), sum(1 - oof_avg$y),
    as.numeric(pROC::auc(pROC::roc(oof_avg$y, oof_avg$indelvar, quiet = TRUE,
                                   levels = c(0, 1), direction = "<"))))

med <- sapply(feats, function(f) {
  v <- suppressWarnings(as.numeric(full[[f]])); if (all(is.na(v))) 0 else median(v, na.rm = TRUE) })
saveRDS(med, file.path(OUT, "feature_medians.rds"))
say("DONE training. artifacts -> %s", OUT)
