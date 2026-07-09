#!/usr/bin/env Rscript
## Figure S4: test-set precision-recall curves by indel class, INDELVAR vs comparator tools.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(PRROC)
  library(patchwork)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")
infile <- file.path(ROOT, "dataset/test_data/test_set_tool_benchmark_scored.tsv")

if (!file.exists(infile)) {
  message(sprintf(paste0(
    "no scored test table at:\n  %s\n\n",
    "This figure needs the test set scored by INDELVAR and the comparator tools.\n",
    "The cross-source test set is HGMD-derived and licensed, so it ships with no\n",
    "pathogenic variants; produce it under your own HGMD license, then rerun.\n",
    "See README 'Reproduce'."), infile))
  quit(save = "no", status = 0)
}

COL_INDELVAR <- "#FD6467"

# per-tool colours shared with Fig 4C so each tool keeps one colour across figures
TOOL_COLORS <- c(
  INDELVAR = "#FD6467", `FATHMM-indel` = "#2A9D8F", CADD = "#E0A100",
  `ESM-1b` = "#8E6FB3", ProGen2 = "#A0703C", phyloP100 = "#4F86C6",
  `GERP++` = "#E8743B", SHINE = "#7A8A99", `PON-Del` = "#5E9C76",
  IndeLLM = "#A23E8C"
)
COL_TEXT <- "#1F2933"

theme_indelvar <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Inter") +
    theme(
      plot.title       = element_text(size = base_size, face = "bold",
                                      colour = COL_TEXT),
      axis.title       = element_text(size = base_size, face = "plain",
                                      colour = COL_TEXT),
      axis.text        = element_text(size = base_size - 1, colour = COL_TEXT),
      legend.title     = element_text(size = base_size - 1, face = "bold",
                                      colour = COL_TEXT),
      legend.text      = element_text(size = base_size - 2, colour = COL_TEXT),
      legend.key.size  = unit(3, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(colour = "#1F2933", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# tool set; oriented higher = pathogenic per tool by pathogenic-vs-benign mean
TOOLS <- c("INDELVAR", "IndeLLM", "FATHMM-indel", "CADD", "ESM-1b", "ProGen2",
           "phyloP100", "GERP++", "SHINE", "PON-Del")
DISPLAY <- c(INDELVAR = "INDELVAR", IndeLLM = "IndeLLM",
             `FATHMM-indel` = "FATHMM-indel",
             CADD = "CADD v1.7", `ESM-1b` = "ESM-1b", ProGen2 = "ProGen2",
             phyloP100 = "phyloP", `GERP++` = "GERP++", SHINE = "SHINE",
             `PON-Del` = "PON-Del")
# leakage-flagged tools (trained on test sources); marked with a trailing * in the legend
LEAK_TOOLS <- c("SHINE", "PON-Del", "IndeLLM")

m <- fread(infile)
cat(sprintf("[load] test set n = %d (del=%d ins=%d)\n", nrow(m),
            sum(m$indel_type == "deletion"), sum(m$indel_type == "insertion")))

pr_for_tool <- function(d, tool) {
  s <- suppressWarnings(as.numeric(d[[tool]])); y <- d$label
  keep <- !is.na(s) & !is.na(y); s <- s[keep]; y <- y[keep]
  if (length(s) < 20 || length(unique(y)) < 2) return(NULL)
  if (sum(y == "Pathogenic") < 3 || sum(y == "Benign") < 3) return(NULL)
  if (mean(s[y == "Pathogenic"]) < mean(s[y == "Benign"])) s <- -s
  pr <- pr.curve(scores.class0 = s[y == "Pathogenic"],
                 scores.class1 = s[y == "Benign"], curve = TRUE)
  data.table(tool = tool, recall = pr$curve[, 1], precision = pr$curve[, 2],
             auprc = pr$auc.integral, n = length(s))
}

build_panel <- function(sub, title) {
  d0 <- m[indel_type == sub]
  num <- function(x) suppressWarnings(as.numeric(x))
  # AUPRC is prevalence-dependent: restrict to tools with >= 50% coverage in
  # this class and evaluate on the common subset they all cover
  shown <- TOOLS[sapply(TOOLS, function(t) mean(!is.na(num(d0[[t]]))) >= 0.5)]
  d <- d0[Reduce(`&`, lapply(shown, function(t) !is.na(num(d0[[t]]))))]
  n_coh <- nrow(d); prev <- mean(d$label == "Pathogenic")
  pr_dt <- rbindlist(lapply(shown, function(t) pr_for_tool(d, t)))
  tab <- unique(pr_dt[, .(tool, auprc, n)])[order(-auprc)]
  ordered <- tab$tool
  col_vec <- TOOL_COLORS[ordered]
  lab_vec <- sprintf("%s, %.3f", DISPLAY[ordered], tab$auprc)
  lab_vec <- setNames(paste0(lab_vec, ifelse(ordered %in% LEAK_TOOLS, " *", "")), ordered)
  pr_dt[, tool := factor(tool, levels = ordered)]
  ggplot(pr_dt, aes(x = recall, y = precision, colour = tool, linewidth = tool)) +
    geom_hline(yintercept = prev, linetype = "22", colour = "#9AA3AD",
               linewidth = 0.35) +
    annotate("text", x = 0.03, y = prev + 0.03, hjust = 0, size = 2.2,
             colour = "#586271", family = "Inter",
             label = sprintf("prevalence = %.2f", prev)) +
    geom_line() +
    scale_colour_manual(values = col_vec, labels = lab_vec,
                        name = "Tool, AUPRC") +
    scale_linewidth_manual(values = setNames(
      ifelse(ordered == "INDELVAR", 1.05, 0.5), ordered), guide = "none") +
    scale_x_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1), expand = c(0, 0)) +
    coord_fixed() +
    labs(x = "Recall", y = "Precision",
         title = sprintf("%s (n = %s)", title,
                         formatC(n_coh, big.mark = ",", format = "d"))) +
    theme_indelvar() +
    theme(legend.position = "right",
          legend.key.height = unit(3.0, "mm"),
          legend.text = element_text(size = 6.6, colour = COL_TEXT, family = "Inter"),
          legend.title = element_text(size = 7),
          plot.margin = margin(6, 8, 6, 8))
}

cat("[AUPRC] computing del / ins panels\n")
pDel <- build_panel("deletion",  "Deletions")
pIns <- build_panel("insertion", "Insertions")
figS4 <- pDel / pIns + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = COL_TEXT))

FIG_NAME <- "figS4_pr_curves"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, paste0(FIG_NAME, ".pdf")), figS4, width = 200, height = 230,
       units = "mm", device = "pdf", dpi = 600)
ggsave(file.path(OUT, paste0(FIG_NAME, ".png")), figS4, width = 200, height = 230,
       units = "mm", device = "png", dpi = 600)
cat("[done] wrote figS4_pr_curves.{pdf,png}\n")
