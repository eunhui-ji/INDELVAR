#!/usr/bin/env Rscript
## Figure 4: ROC curves + AUROC bar chart, INDELVAR vs comparator tools.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(ggtext); library(pROC)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[1] else
  file.path(ROOT, "dataset/test_data/test_set_tool_benchmark_scored.tsv")
prefix <- if (length(args) >= 2) args[2] else file.path(OUT, "indelvar_benchmark")

if (!file.exists(infile)) {
  message(sprintf(paste0(
    "no scored test table at:\n  %s\n\n",
    "This figure needs a test set scored by INDELVAR and the comparator tools\n",
    "(columns: label, indelvar, <tool columns>). The cross-source test set is\n",
    "HGMD-derived and licensed, so it ships with no pathogenic variants; produce\n",
    "it with benchmark.R under your own HGMD license, then:\n",
    "  Rscript plot_benchmark_roc.R <scored.tsv>\nSee README 'Reproduce'."), infile))
  quit(save = "no", status = 0)
}

COL_INDELVAR <- "#FD6467"; COL_TEXT <- "#1F2933"; COL_NEUT <- "#9AA3AD"
theme_indelvar <- function(base = 10) theme_minimal(base_size = base) + theme(
  plot.title   = element_text(size = base + 1, face = "bold", colour = COL_TEXT,
                              hjust = 0.5, margin = margin(b = 5)),
  axis.title   = element_text(colour = COL_TEXT), axis.text = element_text(colour = COL_TEXT),
  legend.title = element_text(face = "bold", colour = COL_TEXT),
  legend.text  = element_text(colour = COL_TEXT), legend.key.size = unit(3.5, "mm"),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
  panel.border = element_rect(colour = COL_TEXT, fill = NA, linewidth = 0.4),
  plot.background = element_rect(fill = "white", colour = NA))

d <- fread(infile)
if (!"label" %in% names(d)) stop("input must contain a 'label' column", call. = FALSE)
d[, y := as.integer(label %in% c("Pathogenic", "P_LP", "P/LP") | label == 1)]
ivar <- grep("^indelvar$", names(d), ignore.case = TRUE, value = TRUE)[1]
if (is.na(ivar)) stop("input must contain an 'indelvar' score column", call. = FALSE)

## comparator columns: numeric, excluding identifiers/metadata
meta  <- c("y", "aa_start", "indel_size_aa", "pos", "acmg_points", "n_mechanisms", "cov", "leak")
tools <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], meta)
tools <- c(ivar, setdiff(tools, ivar))            # INDELVAR first

## AUROC + ROC per tool, auto-oriented score direction
roc1 <- function(col) {
  s <- suppressWarnings(as.numeric(d[[col]])); k <- is.finite(s) & !is.na(d$y)
  if (length(unique(d$y[k])) < 2 || sum(k) < 10) return(NULL)
  r <- roc(d$y[k], s[k], quiet = TRUE, direction = "<", levels = c(0, 1))
  if (as.numeric(auc(r)) < 0.5)
    r <- roc(d$y[k], s[k], quiet = TRUE, direction = ">", levels = c(0, 1))
  r
}
rocs <- Filter(Negate(is.null), setNames(lapply(tools, roc1), tools))
if (!length(rocs)) stop("no scorable columns (need >=2 classes and >=10 labelled rows)", call. = FALSE)
au  <- vapply(rocs, function(r) as.numeric(auc(r)), numeric(1))
ord <- names(sort(au, decreasing = TRUE))
cat(sprintf("[roc] %d tools on n=%d (P=%d B=%d):\n", length(rocs), sum(!is.na(d$y)),
            sum(d$y == 1, na.rm = TRUE), sum(d$y == 0, na.rm = TRUE)))
for (t in ord) cat(sprintf("   %-16s AUROC = %.3f\n", t, au[t]))

lab_of  <- setNames(sprintf("%s (%.3f)", ord, au[ord]), ord)
is_ivar <- function(t) tolower(t) == "indelvar"
## INDELVAR rose, comparators muted grey
comp    <- setdiff(ord, ord[vapply(ord, is_ivar, logical(1))])
comp_col <- setNames(colorRampPalette(c("#4A5568", "#B7C0CA"))(max(length(comp), 1)), comp)
pal     <- c(comp_col, setNames(COL_INDELVAR, ord[vapply(ord, is_ivar, logical(1))]))
names(pal) <- lab_of[names(pal)]

## ---- ROC curves ----
roc_df <- rbindlist(lapply(ord, function(t) {
  co <- coords(rocs[[t]], "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  data.table(tool = lab_of[t], fpr = 1 - co$specificity, tpr = co$sensitivity,
             ivar = is_ivar(t))
}))
roc_df[, tool := factor(tool, levels = lab_of[ord])]
pR <- ggplot(roc_df, aes(fpr, tpr, colour = tool, linewidth = ivar)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = COL_NEUT, linewidth = 0.4) +
  geom_path() +
  scale_colour_manual(values = pal, name = NULL) +
  scale_linewidth_manual(values = c(`FALSE` = 0.55, `TRUE` = 1.3), guide = "none") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  coord_equal() +
  labs(x = "1 - specificity", y = "Sensitivity", title = "Cross-source ROC") +
  theme_indelvar() + theme(legend.position = "right")
ggsave(paste0(prefix, "_roc.png"), pR, width = 175, height = 130, units = "mm", dpi = 600)

## ---- AUROC bar chart ----
bar <- data.table(tool = factor(ord, levels = rev(ord)), auroc = au[ord],
                  ivar = vapply(ord, is_ivar, logical(1)))
pB <- ggplot(bar, aes(auroc, tool, fill = ivar)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = sprintf("%.3f", auroc)), hjust = -0.15, size = 3, colour = COL_TEXT) +
  scale_fill_manual(values = c(`FALSE` = COL_NEUT, `TRUE` = COL_INDELVAR), guide = "none") +
  scale_x_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(x = "AUROC", y = NULL, title = "Cross-source AUROC") +
  theme_indelvar()
ggsave(paste0(prefix, "_auroc.png"), pB, width = 150, height = 120, units = "mm", dpi = 600)
cat(sprintf("[done] wrote %s_roc.png and %s_auroc.png\n", prefix, prefix))
