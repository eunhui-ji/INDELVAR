#!/usr/bin/env Rscript
## Figure S6: calibration (reliability) and decision-curve analysis on the test set.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(patchwork)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")
infile <- file.path(ROOT, "dataset/test_data/test_set_scored.tsv")

if (!file.exists(infile)) {
  message(sprintf(paste0(
    "no scored test table at:\n  %s\n\n",
    "This figure needs the test set scored by INDELVAR (columns: variant_id,\n",
    "label, indel_type, indelvar). The cross-source test set is HGMD-derived\n",
    "and licensed, so it ships with no pathogenic variants; produce it under\n",
    "your own HGMD license, then rerun. See README 'Reproduce'."), infile))
  quit(save = "no", status = 0)
}

set.seed(42)

COL_INDELVAR <- "#FD6467"; COL_PATH <- "#E47DA3"; COL_BENIGN <- "#7294D4"
COL_TEXT <- "#1F2933"; COL_NEUT <- "#9AA5B1"
theme_indelvar <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Inter") +
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold", colour = COL_TEXT, hjust = 0.5),
      axis.title       = element_text(size = base_size, colour = COL_TEXT),
      axis.text        = element_text(size = base_size - 1, colour = COL_TEXT),
      legend.title      = element_blank(),
      legend.text       = element_text(size = base_size - 1, colour = COL_TEXT),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(colour = "#1F2933", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA))
}

d <- fread(infile)
d[, y := as.integer(label == "Pathogenic")]
d[, p := indelvar]

brier <- function(dt) mean((dt$p - dt$y)^2)
BR <- c(Pooled = brier(d), Deletions = brier(d[indel_type == "deletion"]),
        Insertions = brier(d[indel_type == "insertion"]))

## ---- Panel A: reliability curve (10 fixed bins, Wilson 95% CI) ----------------
d[, bin := cut(p, breaks = seq(0, 1, 0.1), include.lowest = TRUE)]
rel <- d[, {
  n <- .N; k <- sum(y); ph <- k / n
  z <- 1.96; den <- 1 + z^2 / n
  centre <- (ph + z^2 / (2 * n)) / den
  half <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / den
  .(mean_p = mean(p), obs = ph, lo = pmax(0, centre - half),
    hi = pmin(1, centre + half), n = n)
}, by = bin][order(mean_p)]

pA <- ggplot(rel, aes(mean_p, obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = COL_NEUT, linewidth = 0.4) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.02, colour = COL_INDELVAR, linewidth = 0.4) +
  geom_line(colour = COL_INDELVAR, linewidth = 0.6) +
  geom_point(aes(size = n), colour = COL_INDELVAR) +
  scale_size_area(max_size = 3.2, guide = "none") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Calibration (reliability)",
       x = "Mean predicted probability", y = "Observed fraction pathogenic") +
  theme_indelvar()

## ---- Panel B: decision-curve analysis (net benefit) --------------------------
n <- nrow(d); prev <- mean(d$y); pt <- seq(0.01, 0.99, 0.01)
nb_model <- sapply(pt, function(t) {
  pred <- d$p >= t
  tp <- sum(pred & d$y == 1); fp <- sum(pred & d$y == 0)
  tp / n - (fp / n) * (t / (1 - t))
})
nb_all <- prev - (1 - prev) * (pt / (1 - pt))
dca <- rbind(
  data.table(pt = pt, nb = nb_model, k = "INDELVAR"),
  data.table(pt = pt, nb = nb_all,   k = "Treat all"),
  data.table(pt = pt, nb = 0,        k = "Treat none"))
dca[, k := factor(k, levels = c("INDELVAR", "Treat all", "Treat none"))]

pB <- ggplot(dca, aes(pt, nb, colour = k, linetype = k)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c("INDELVAR" = COL_INDELVAR, "Treat all" = COL_BENIGN,
                                 "Treat none" = COL_NEUT)) +
  scale_linetype_manual(values = c("INDELVAR" = 1, "Treat all" = 2, "Treat none" = 3)) +
  coord_cartesian(ylim = c(-0.05, max(nb_model) * 1.05)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(title = "Decision curve", x = "Threshold probability", y = "Net benefit") +
  theme_indelvar() +
  theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        legend.key.height = unit(3, "mm"))

fig <- pA + pB + plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", family = "Inter", colour = COL_TEXT))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT, "figS6_calibration_decision.png"), fig,
       width = 7.2, height = 3.3, dpi = 600)
ggsave(file.path(OUT, "figS6_calibration_decision.pdf"), fig,
       width = 7.2, height = 3.3)
cat(sprintf("[figS6] Brier pooled=%.4f del=%.4f ins=%.4f ; prevalence=%.3f ; n=%d\n",
            BR["Pooled"], BR["Deletions"], BR["Insertions"], prev, n))
