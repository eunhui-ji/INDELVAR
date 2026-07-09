#!/usr/bin/env Rscript
## benchmark.R: test-set per-tool AUROC and INDELVAR-vs-FATHMM leakage-free DeLong.
## Usage: Rscript benchmark.R
suppressPackageStartupMessages({library(data.table); library(pROC)})
S <- Sys.getenv("INDELVAR_ROOT", ".")
CAN <- file.path(S, "results/figtable_data")
need <- function(p) { if (!file.exists(p)) stop(sprintf(
  "missing input: %s\nbenchmark.R reads development-repo comparator scores and test-set tables (results/, tools/, data/) that are not shipped in the release; see README 'Reproduce'.", p), call. = FALSE); p }

## (1) per-tool AUROC benchmark
invisible((function() {
cx <- fread(need(file.path(CAN, "crosssource_benchmark_scored.tsv")))
m <- unique(cx[, .(variant_id, label, indel_type,
                   INDELVAR = indelvar, CADD = cadd, `ESM-1b` = esm1b)], by = "variant_id")
con <- unique(fread(file.path(S, "data/processed/hgmd_test_conservation_subcell.tsv")),
              by = "variant_id")
m <- merge(m, con[, .(variant_id, phyloP100 = phylop100_mean, `GERP++` = gerp_mean)],
           by = "variant_id", all.x = TRUE)
fa <- fread(file.path(CAN, "external_tool_scores/fathmm_hgmd_test.tsv"))
m <- merge(m, unique(fa[, .(variant_id, `FATHMM-indel` = score_raw)], by = "variant_id"),
           by = "variant_id", all.x = TRUE)
pg <- fread(file.path(CAN, "external_tool_scores/progen2_hgmd_test.tsv"))
m <- merge(m, unique(pg[, .(variant_id, ProGen2 = score_raw)], by = "variant_id"),
           by = "variant_id", all.x = TRUE)

## SHINE (deletion-only)
g2u <- as.data.table(readRDS(file.path(S, "data/processed/gene_to_uniprot_expanded.rds")))
u2g <- unique(g2u[!is.na(uniprot) & !is.na(gene_symbol), .(uniprot, gene_symbol, ensg)], by = "uniprot")
shmap <- fread(file.path(S, "tools/SHINE/MSA/ensembl_gene_protein_id.txt"))
g2e <- setNames(shmap[Gene.name != "", Protein.stable.ID], shmap[Gene.name != "", Gene.name])
eg2e <- setNames(shmap[Gene.stable.ID != "", Protein.stable.ID], shmap[Gene.stable.ID != "", Gene.stable.ID])
PRE <- file.path(S, "tools/SHINE/precomputed/deletion")
cache <- new.env()
load_ensp <- function(e) {
  if (is.null(e) || is.na(e) || e == "") return(NULL)
  if (exists(e, cache)) return(get(e, cache))
  p <- file.path(PRE, paste0(e, ".txt")); d <- NULL
  if (file.exists(p)) { x <- fread(p, header = TRUE); d <- setNames(x[[3]], x[[1]]) }
  assign(e, d, cache); d
}
td <- unique(pg[indel_type == "deletion",
               .(variant_id, uniprot, ps = protein_pos_start, pe = protein_pos_end)],
             by = "variant_id")
td <- merge(td, u2g, by = "uniprot", all.x = TRUE)
td[, ensp := g2e[gene_symbol]]
td[is.na(ensp), ensp := eg2e[ensg]]
shine_score <- function(e, ps, pe) {
  d <- load_ensp(e); if (is.null(d)) return(NA_real_)
  v <- d[as.character(ps:pe)]; v <- v[!is.na(v)]
  if (!length(v)) return(NA_real_); mean(v)
}
td[, SHINE := mapply(shine_score, ensp, ps, pe)]
m <- merge(m, td[, .(variant_id, SHINE)], by = "variant_id", all.x = TRUE)
## PON-Del (deletion-only)
pd <- fread(file.path(CAN, "external_tool_scores/pondel_hgmd_test.tsv"))
m <- merge(m, unique(pd[, .(variant_id, `PON-Del` = score_raw)], by = "variant_id"),
           by = "variant_id", all.x = TRUE)
## IndeLLM
il <- fread(file.path(S, "data/comparators/indellm_scores.txt"))
m <- merge(m, unique(il[, .(variant_id = id, IndeLLM = score)], by = "variant_id"),
           by = "variant_id", all.x = TRUE)
m[, y := as.integer(label == "Pathogenic")]
cat(sprintf("assembled n=%d (P=%d B=%d; del=%d ins=%d); SHINE on %d deletions\n",
            nrow(m), sum(m$y), sum(1 - m$y), sum(m$indel_type == "deletion"),
            sum(m$indel_type == "insertion"), sum(!is.na(m$SHINE))))

TOOLS <- c("INDELVAR", "IndeLLM", "FATHMM-indel", "CADD", "ESM-1b", "ProGen2", "phyloP100", "GERP++", "SHINE", "PON-Del")
# TRUE = tool has training overlap with this test set
LEAK <- c(INDELVAR = FALSE, IndeLLM = TRUE, `FATHMM-indel` = FALSE, CADD = FALSE,
          `ESM-1b` = FALSE, ProGen2 = FALSE, phyloP100 = FALSE, `GERP++` = FALSE,
          SHINE = TRUE, `PON-Del` = TRUE)
au <- function(d, t) {
  s <- suppressWarnings(as.numeric(d[[t]])); k <- !is.na(s) & !is.na(d$y); y <- d$y[k]; s <- s[k]
  if (length(unique(y)) < 2 || sum(y) < 3 || sum(1 - y) < 3) return(NULL)
  r <- roc(y, s, quiet = TRUE, levels = c(0, 1), direction = "<"); a <- as.numeric(auc(r))
  if (a < 0.5) { r <- roc(y, s, quiet = TRUE, levels = c(0, 1), direction = ">"); a <- as.numeric(auc(r)) }
  ci <- as.numeric(ci.auc(r))
  data.table(AUROC = a, lo = ci[1], hi = ci[3], n = length(y), nP = sum(y), nB = sum(1 - y))
}
summ <- list()
for (sub in c("all", "deletion", "insertion")) {
  d <- if (sub == "all") m else m[indel_type == sub]
  for (t in TOOLS) { a <- au(d, t)
    if (!is.null(a)) summ[[length(summ) + 1]] <- cbind(subset = sub, tool = t, leak = LEAK[[t]],
                                                       cov = round(100 * a$n / nrow(d)), a) }
}
S2 <- rbindlist(summ)
fwrite(m, file.path(CAN, "test_set_tool_benchmark_scored.tsv"), sep = "\t")
fwrite(S2, file.path(CAN, "test_set_tool_benchmark_auroc.tsv"), sep = "\t")

cat("\nVALIDATION: ALL-set should match Fig-4C (INDELVAR .927 IndeLLM .925 FATHMM .918 phyloP .865 CADD .855 ESM .815 GERP .753 ProGen2 .733)\n")
for (sub in c("all", "deletion", "insertion")) {
  cat(sprintf("\n== %s ==\n", toupper(sub)))
  print(S2[subset == sub][order(-AUROC), .(tool, AUROC = round(AUROC, 3),
        CI = sprintf("%.3f-%.3f", lo, hi), cov = paste0(cov, "%"), nP, nB, leak)])
}
cat("\nwrote results/figtable_data/test_set_tool_benchmark_{scored,auroc}.tsv\n")
})())

## (2) FATHMM-indel leakage-free DeLong (seed=42, pROC)
invisible((function() {
set.seed(42)
st <- fread(need(file.path(CAN, "test_set_scored.tsv")))
fa <- unique(fread(need(file.path(S, "results/benchmark/external_tool_scores/fathmm_hgmd_test.tsv")))[, .(variant_id, fathmm = score_raw)], by = "variant_id")
d  <- merge(st, fa, by = "variant_id")[!is.na(fathmm) & !is.na(indelvar)]
ov <- fread(need(file.path(CAN, "leakage/fathmm_test_overlap.tsv")))$variant_id
d[, leak := variant_id %in% ov]
cat(sprintf("merged n=%d (P=%d B=%d); overlap removed=%d; pROC %s\n",
            nrow(d), sum(d$label == "Pathogenic"), sum(d$label == "Benign"),
            sum(d$leak), as.character(packageVersion("pROC"))))
delong <- function(x) {
  rs <- roc(x$label, x$indelvar, levels = c("Benign", "Pathogenic"), direction = "<", quiet = TRUE)
  rf <- roc(x$label, x$fathmm, levels = c("Benign", "Pathogenic"), direction = "<", quiet = TRUE)
  list(t = roc.test(rs, rf, method = "delong", paired = TRUE),
       as = as.numeric(auc(rs)), af = as.numeric(auc(rf)))
}
for (sub in c("all", "deletion", "insertion"))
  for (setn in c("full", "leakage-free")) {
    x <- d[(sub == "all" | indel_type == sub) & (setn == "full" | !leak)]
    if (length(unique(x$label)) < 2) next
    r <- delong(x)
    cat(sprintf("[%-9s %-12s] n=%4d  INDELVAR=%.3f  FATHMM-indel=%.3f  dAUROC=%+.3f  DeLong p=%.3f\n",
                sub, setn, nrow(x), r$as, r$af, r$as - r$af, r$t$p.value))
  }
})())
