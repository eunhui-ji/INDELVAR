#!/usr/bin/env Rscript
## Score in-frame indels: raw RF Pathogenic probability, ACMG/AMP stratum, and SHAP
## feature-group attribution. Input needs the 39 features plus an indel_type column.
## Usage: Rscript score_variants.R <in.{parquet,tsv,csv}> <out.tsv> [del|ins] [--nsim N] [--per_feature]
suppressWarnings(suppressPackageStartupMessages({
  library(data.table); library(randomForest); library(caret); library(arrow)
}))

.indelvar_env <- new.env(parent = emptyenv())

.INDELVAR_GRP <- list(
  Structure    = c("af_plddt_mean","af_helix_frac","af_sheet_frac","af_rsa_mean",
                   "low_plddt_region","wcn_mean","contact_density","bridging_contacts",
                   "ss_element_break","ss_element_frac_affected","helix_register_preserved",
                   "hbond_density","hydrophobic_exposure_risk","local_ss_class"),
  Conservation = c("phylop100_mean","gerp_mean"),
  SeqContext   = c("local_entropy","dist_to_splice","indel_type","indel_size_aa",
                   "indel_size_relative"),
  Constraint   = c("pli","loeuf","mis_z","cds_length"),
  Annotation   = c("in_domain","n_domains","in_repeat","ptm_nearby","ptm_count_near",
                   "disulfide_nearby","in_active_site","subcell_Cytoplasm","subcell_Membrane",
                   "subcell_Mitochondrion","subcell_Nucleus","subcell_Secreted",
                   "subcell_Other","subcell_Unknown"))
.INDELVAR_GRP_ORDER <- c("Structure","Conservation","SeqContext","Constraint","Annotation")

indelvar_load <- function(model_dir = Sys.getenv("INDELVAR_MODEL",
                                               file.path(Sys.getenv("INDELVAR_ROOT","."), "model")),
                        n_background = 0L) {
  e <- .indelvar_env
  e$rf    <- readRDS(file.path(model_dir, "indelvar.rds"))
  e$feats <- readRDS(file.path(model_dir, "col_definitions.rds"))$feature_cols
  e$cut   <- fread(file.path(model_dir, "evidence_cutoffs.tsv"))
  e$grp   <- .INDELVAR_GRP
  tr <- fread(file.path(Sys.getenv("INDELVAR_ROOT","."), "dataset", "train_data", "train_set.tsv"))
  mf <- file.path(model_dir, "feature_medians.rds")
  if (file.exists(mf)) e$med <- readRDS(mf) else {
    e$med <- sapply(e$feats, function(f){
      v <- suppressWarnings(as.numeric(tr[[f]])); if (all(is.na(v))) 0 else median(v, na.rm = TRUE) })
    saveRDS(e$med, mf)
  }
  bg <- .prep_X(tr)                                    # SHAP background = training distribution
  if (n_background > 0L && n_background < nrow(bg)) {
    set.seed(1L); bg <- bg[sample.int(nrow(bg), n_background), , drop = FALSE]
  }
  e$bg <- bg
  invisible(e)
}

.getco <- function(gg) {
  c <- .indelvar_env$cut[subset == gg]
  g <- function(ev){ v <- c[evidence == ev]$cutoff; if (length(v) == 0) NA_real_ else v }
  list(PSt = g("PP3_Strong"), PMo3 = g("PP3_Moderate3"), PMo = g("PP3_Moderate"),
       PSu = g("PP3_Supporting"), BSu = g("BP4_Supporting"), BMo = g("BP4_Moderate"),
       BMo3 = g("BP4_Moderate3"), BSt = g("BP4_Strong"))
}

# 8-tier assignment, weakest first so the strongest applicable tier wins; NA cutoffs skipped.
# PP3_Moderate spans +2 and +3 (BP4_Moderate -2/-3); use acmg_points to tell them apart.
.catv <- function(score, gg) {
  k <- .getco(gg); n <- length(score)
  cat <- rep("Uncertain", n); pts <- rep(0L, n); ok <- !is.na(score)
  if (!is.na(k$BSu))  { i <- ok & score <= k$BSu;  cat[i] <- "BP4_Supporting"; pts[i] <- -1L }
  if (!is.na(k$BMo))  { i <- ok & score <= k$BMo;  cat[i] <- "BP4_Moderate";   pts[i] <- -2L }
  if (!is.na(k$BMo3)) { i <- ok & score <= k$BMo3; cat[i] <- "BP4_Moderate";   pts[i] <- -3L }
  if (!is.na(k$BSt))  { i <- ok & score <= k$BSt;  cat[i] <- "BP4_Strong";     pts[i] <- -4L }
  if (!is.na(k$PSu))  { i <- ok & score >= k$PSu;  cat[i] <- "PP3_Supporting"; pts[i] <-  1L }
  if (!is.na(k$PMo))  { i <- ok & score >= k$PMo;  cat[i] <- "PP3_Moderate";   pts[i] <-  2L }
  if (!is.na(k$PMo3)) { i <- ok & score >= k$PMo3; cat[i] <- "PP3_Moderate";   pts[i] <-  3L }
  if (!is.na(k$PSt))  { i <- ok & score >= k$PSt;  cat[i] <- "PP3_Strong";     pts[i] <-  4L }
  cat[!ok] <- NA_character_; pts[!ok] <- NA_integer_
  list(category = cat, points = pts)
}

.prep_X <- function(df) {
  e <- .indelvar_env; df <- as.data.table(df); fc <- e$feats
  for (c in fc) {
    if (!(c %in% names(df))) df[, (c) := NA_real_]
    if (c != "indel_type") {
      df[[c]] <- suppressWarnings(as.numeric(df[[c]]))
      if (any(is.na(df[[c]]))) df[is.na(get(c)), (c) := e$med[[c]]]
    }
  }
  X <- as.data.frame(df[, ..fc])
  if ("indel_type" %in% fc)
    X$indel_type <- factor(ifelse(df$indel_type == "deletion" | df$indel_type == 1,
                                  "deletion", "insertion"), levels = c("deletion", "insertion"))
  X
}

.pw <- function(X) predict(.indelvar_env$rf, X, type = "prob")[, "Pathogenic"]

# per-feature permutation Shapley: each sim reveals features in a random order from a random
# background reference, accumulating prediction deltas; sum of phi = f(x) - f(reference).
.perm_shap <- function(Xexp, nsim = 150L, seed = 1L) {
  e <- .indelvar_env; Xbg <- e$bg; feats <- e$feats
  set.seed(seed)
  n <- nrow(Xexp); p <- length(feats); m <- nrow(Xbg)
  phi <- matrix(0, n, p, dimnames = list(NULL, feats))
  for (s in seq_len(nsim)) {
    o   <- sample.int(p)
    cur <- Xbg[sample.int(m, n, replace = TRUE), , drop = FALSE]
    prev <- .pw(cur)
    for (j in o) {
      cur[[j]] <- Xexp[[j]]
      pj <- .pw(cur)
      phi[, j] <- phi[, j] + (pj - prev)
      prev <- pj
    }
  }
  phi / nsim
}

indelvar_explain <- function(df, subset = c("del", "ins"), with_shap = TRUE,
                           nsim = 150L, per_feature = FALSE) {
  subset <- match.arg(subset); e <- .indelvar_env
  if (is.null(e$rf)) stop("call indelvar_load() first")
  if (subset == "ins") { df <- as.data.table(df); df[, helix_register_preserved := NA_real_] }  # undefined for inserted residues
  Xexp  <- .prep_X(df)
  score <- round(.pw(Xexp), 6)
  cc <- .catv(score, subset)
  out <- data.table(indelvar_score = score, indelvar_category = cc$category, acmg_points = cc$points)
  if (with_shap) {
    phi <- .perm_shap(Xexp, nsim = nsim)
    gsum <- function(g) {
      cols <- intersect(e$grp[[g]], colnames(phi))
      if (!length(cols)) return(rep(0, nrow(phi)))
      rowSums(phi[, cols, drop = FALSE])
    }
    out[, `:=`(shap_structure    = round(gsum("Structure"), 6),
               shap_conservation = round(gsum("Conservation"), 6),
               shap_seqcontext   = round(gsum("SeqContext"), 6),
               shap_constraint   = round(gsum("Constraint"), 6),
               shap_annotation   = round(gsum("Annotation"), 6))]
    gm <- abs(as.matrix(out[, .(shap_structure, shap_conservation, shap_seqcontext,
                                shap_constraint, shap_annotation)]))
    out[, top_driver_group := .INDELVAR_GRP_ORDER[max.col(gm, ties.method = "first")]]
    out[, top3_features := vapply(seq_len(nrow(phi)), function(r) {
      v <- phi[r, ]; o <- order(abs(v), decreasing = TRUE)[1:3]
      paste(sprintf("%s(%+.3f)", names(v)[o], v[o]), collapse = "; ")
    }, character(1))]
    if (per_feature) {
      pf <- as.data.table(round(phi, 6)); setnames(pf, paste0("shap_", colnames(phi)))
      out <- cbind(out, pf)
    }
  }
  out[]
}

# ---- CLI ----
args <- commandArgs(trailingOnly = TRUE)
nsim <- 150L; per_feature <- FALSE; keep <- rep(TRUE, length(args))
i <- which(args == "--nsim")
if (length(i)) { nsim <- as.integer(args[i + 1]); keep[c(i, i + 1)] <- FALSE }
i <- which(args == "--per_feature")
if (length(i)) { per_feature <- TRUE; keep[i] <- FALSE }
pos <- args[keep]
if (length(pos) < 2)
  stop("usage: Rscript score_variants.R <input> <output.tsv> [del|ins] [--nsim N] [--per_feature]")
infile <- pos[1]; outfile <- pos[2]
sub <- if (length(pos) >= 3) pos[3] else NA_character_

read_any <- function(p) {
  if (grepl("\\.parquet$", p)) as.data.table(read_parquet(p))
  else as.data.table(fread(p))
}
d <- read_any(infile)
if (!"indel_type" %in% names(d)) {
  if (is.na(sub)) stop("input needs an 'indel_type' column, or pass 'del' or 'ins' as the 3rd argument")
  d[, indel_type := if (sub == "del") "deletion" else "insertion"]
}
norm_sub <- function(x) ifelse(x == "deletion" | x == 1, "del", "ins")
d[, .subset := if (!is.na(sub)) sub else norm_sub(indel_type)]

cat(sprintf("INDELVAR: scoring %s variants (%s) ...\n",
            format(nrow(d), big.mark = ","),
            paste(sprintf("%s=%d", names(table(d$.subset)), table(d$.subset)), collapse = ", ")))
indelvar_load()
res <- rbindlist(lapply(split(seq_len(nrow(d)), d$.subset), function(idx) {
  s <- d$.subset[idx[1]]
  cbind(data.table(.row = idx),
        indelvar_explain(d[idx], subset = s, with_shap = TRUE, nsim = nsim, per_feature = per_feature))
}))
setorder(res, .row); res[, .row := NULL]
d[, .subset := NULL]
out <- cbind(d, res)
fwrite(out, outfile, sep = "\t")
cat(sprintf("wrote %s (%s rows, %d cols)\n", outfile, format(nrow(out), big.mark = ","), ncol(out)))
