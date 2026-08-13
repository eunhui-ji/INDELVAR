#!/usr/bin/env Rscript
# Score a 37-feature table with the trained random forest.

args <- commandArgs(trailingOnly = TRUE)
get <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing required argument ", flag)
    return(default)
  }
  args[[i + 1L]]
}

FEATURES_P <- get("--features")
OUT_P      <- get("--out")
MODEL_P    <- get("--model")
CONTRACT_P <- get("--contract")
MEDIANS_P  <- get("--medians")
RUNTIME_P  <- get("--runtime-versions")
PROV_P     <- paste0(OUT_P, ".provenance.tsv")

for (p in c(OUT_P, PROV_P))
  if (file.exists(p)) stop("OUTPUT_EXISTS: ", basename(p))

for (p in c(MODEL_P, CONTRACT_P, MEDIANS_P, RUNTIME_P))
  if (!file.exists(p)) stop("MISSING_ARTIFACT: ", basename(p))

# Publish both outputs only after they are complete.
TMP_OUT  <- tempfile("indelvar", dirname(OUT_P), ".part")
TMP_PROV <- tempfile("indelvar", dirname(PROV_P), ".part")
created <- c(TMP_OUT, TMP_PROV)
finished <- FALSE
cleanup <- function() {
  if (!finished) for (f in created) if (file.exists(f)) unlink(f)
}
options(error = function() { cleanup(); quit(save = "no", status = 1) })

runtime <- read.delim(RUNTIME_P, sep = "\t", header = TRUE, quote = "",
                      colClasses = "character", check.names = FALSE,
                      comment.char = "#")
if (!all(c("package", "version", "consumed_by") %in% names(runtime)))
  stop("MALFORMED_RUNTIME_TABLE: ", basename(RUNTIME_P))
need <- runtime[runtime$consumed_by == "model scoring", , drop = FALSE]
if (nrow(need) == 0L) stop("NO_R_PACKAGE_DECLARED: ", basename(RUNTIME_P))
for (i in seq_len(nrow(need))) {
  pkg <- need$package[i]
  if (!requireNamespace(pkg, quietly = TRUE)) stop("MISSING_R_PACKAGE: ", pkg)
  if (!identical(as.character(packageVersion(pkg)),
                 as.character(package_version(need$version[i]))))
    stop("R_PACKAGE_VERSION_MISMATCH: ", pkg)
}

suppressPackageStartupMessages({
  library(randomForest)
})

sha256 <- function(p) digest::digest(file = p, algo = "sha256")
DIGESTS <- c(model = sha256(MODEL_P), contract = sha256(CONTRACT_P),
             medians = sha256(MEDIANS_P))
if (any(is.na(DIGESTS)) || !all(nchar(DIGESTS) == 64L))
  stop("DIGEST_FAILED: ", paste(names(DIGESTS)[nchar(DIGESTS) != 64L],
                                collapse = ", "))

contract <- read.delim(CONTRACT_P, sep = "\t", header = TRUE, quote = "",
                       colClasses = "character", check.names = FALSE)
FEATS <- contract$feature_id
stopifnot(length(FEATS) == 37L, !anyDuplicated(FEATS))
# `indel_type` enters the model matrix as one binary indicator.
MM_INDEL <- "indel_type_insertion"
COLS <- ifelse(FEATS == "indel_type", MM_INDEL, FEATS)

rf <- readRDS(MODEL_P)
stopifnot("model artifact has no preprocessing" = !is.null(rf$prep),
          "model artifact has no forest" = !is.null(rf$model))
if (!identical(rf$feature_order_raw, FEATS))
  stop("the model artifact's feature order differs from the contract")
if (!identical(rf$design_matrix_columns, COLS))
  stop("the model artifact's model matrix differs from the contract encoding")

# Check the exported medians against the model's.
prep <- rf$prep
{
  md <- read.delim(MEDIANS_P, sep = "\t", header = TRUE, quote = "",
                   colClasses = "character", check.names = FALSE)
  exported <- suppressWarnings(setNames(as.numeric(md$training_median), md$feature))
  if (!setequal(names(exported), COLS))
    stop("imputation_medians.tsv does not cover the model matrix columns")
  for (f in COLS) {
    s <- prep$median[[f]]
    if (is.na(s)) next
    e <- exported[[f]]
    if (is.na(e) || !isTRUE(all.equal(s, e, tolerance = 1e-14)))
      stop(sprintf("median for %s disagrees with the export: %.17g vs %.17g",
                   f, s, e))
  }
}

d <- read.delim(FEATURES_P, sep = "\t", header = TRUE, quote = "",
                colClasses = "character", check.names = FALSE,
                na.strings = character(0))
if (!"variant_id" %in% names(d)) stop("input has no variant_id column")
if (anyDuplicated(names(d))) stop("input has duplicate column names: ",
                                  paste(names(d)[duplicated(names(d))], collapse = ", "))
missing <- setdiff(FEATS, names(d))
if (length(missing)) stop("input is missing contract features: ",
                          paste(missing, collapse = ", "))
if ("label" %in% names(d)) stop("input carries a label column")
if (nrow(d) == 0L) stop("input has no rows")

MISSING_TOKENS <- c("", "NA")
x <- data.frame(row.names = seq_len(nrow(d)))
for (f in FEATS) {
  v <- d[[f]]
  v[v %in% MISSING_TOKENS] <- NA_character_
  if (f == "indel_type") {
    if (any(is.na(v)))
      stop("indel_type is missing")
    bad <- setdiff(unique(v), c("deletion", "insertion"))
    if (length(bad))
      stop("unrecognised indel_type value: ", paste(bad, collapse = ", "))
    x[[MM_INDEL]] <- as.numeric(v == "insertion")
  } else {
    n <- suppressWarnings(as.numeric(v))
    if (any(is.na(n) & !is.na(v)))
      stop("non-numeric value in ", f, ": ",
           paste(utils::head(unique(v[is.na(n) & !is.na(v)]), 3), collapse = ", "))
    x[[f]] <- n
  }
}
stopifnot(identical(names(x), COLS))

# The imputed-feature report names the raw feature.
was_na <- vapply(COLS, function(f) is.na(x[[f]]), logical(nrow(x)))
if (is.null(dim(was_na))) was_na <- matrix(was_na, nrow = nrow(x),
                                           dimnames = list(NULL, COLS))
colnames(was_na) <- FEATS
# Impute missing values using the model medians.
stopifnot(identical(names(x), prep$features))
for (f in prep$features) {
  m <- prep$median[[f]]
  if (is.na(m)) stop("no median for ", f)
  x[[f]][is.na(x[[f]])] <- m
}
xm <- as.matrix(x)
stopifnot(!any(is.na(xm)), all(is.finite(xm)))

score <- predict(rf$model, newdata = xm, type = "prob")[, "pathogenic"]
stopifnot(all(is.finite(score)), all(score >= 0), all(score <= 1))

imputed_names <- apply(was_na, 1L, function(r) paste(FEATS[r], collapse = ","))
out <- data.frame(
  variant_id = d$variant_id,
  indelvar_score = sprintf("%.17g", score),
  n_features_imputed = rowSums(was_na),
  imputed_features = imputed_names,
  stringsAsFactors = FALSE
)
write.table(out, TMP_OUT, sep = "\t", quote = FALSE, row.names = FALSE,
            na = "", fileEncoding = "UTF-8")

prov <- c(
  sprintf("model_sha256\t%s", DIGESTS[["model"]]),
  sprintf("feature_contract_sha256\t%s", DIGESTS[["contract"]]),
  sprintf("imputation_medians_sha256\t%s", DIGESTS[["medians"]]),
  sprintf("imputation_source\tmodel_artifact"),
  sprintf("n_rows\t%d", nrow(out)),
  sprintf("r_version\t%s", R.version.string),
  sprintf("randomForest_version\t%s", as.character(packageVersion("randomForest")))
)
writeLines(prov, TMP_PROV)

for (p in c(OUT_P, PROV_P))
  if (file.exists(p)) stop("OUTPUT_EXISTS: ", basename(p))
if (!file.rename(TMP_OUT, OUT_P)) stop("WRITE_FAILED: ", basename(OUT_P))
created <- c(created, OUT_P)
if (!file.rename(TMP_PROV, PROV_P)) stop("WRITE_FAILED: ", basename(PROV_P))
created <- c(created, PROV_P)

finished <- TRUE
cat(sprintf("scored %d rows -> %s\n", nrow(out), OUT_P))
