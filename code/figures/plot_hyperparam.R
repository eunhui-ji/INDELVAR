# Figure S1: RF hyperparameter tuning and classifier comparison.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

COL_BEST   <- "#446455"
COL_TEXT   <- "#1F2933"
COL_NEUT   <- "#9AA3AD"

theme_indelvar <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Inter") +
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      colour = COL_TEXT),
      axis.title       = element_text(size = base_size, face = "plain",
                                      colour = COL_TEXT),
      axis.text        = element_text(size = base_size - 1, colour = COL_TEXT),
      legend.title     = element_text(size = base_size - 1, face = "bold",
                                      colour = COL_TEXT),
      legend.text      = element_text(size = base_size - 1, colour = COL_TEXT),
      legend.key.size  = unit(3, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(colour = "#1F2933", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# Panel A: classifier comparison
clf <- data.table(
  classifier = c("Random Forest", "SVM (RBF)", "LightGBM", "XGBoost",
                 "Elastic Net", "k-NN", "MLP (2-layer)", "Logistic regression"),
  auroc      = c(0.970, 0.963, 0.956, 0.955, 0.951, 0.947, 0.924, 0.912),
  sd         = c(0.006, 0.007, 0.008, 0.008, 0.009, 0.010, 0.012, 0.014)
)
clf[, classifier := factor(classifier, levels = classifier)]
clf[, is_rf := classifier == "Random Forest"]
clf[, lo := auroc - sd]
clf[, hi := pmin(auroc + sd, 1)]

pA <- ggplot(clf, aes(x = classifier, y = auroc, fill = is_rf)) +
  geom_col(width = 0.65, colour = NA) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                width = 0.18, colour = COL_TEXT, linewidth = 0.4) +
  geom_text(aes(y = hi, label = sprintf("%.3f", auroc)),
            vjust = -0.7, size = 2.6, colour = COL_TEXT,
            family = "Inter", fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = COL_BEST, `FALSE` = COL_NEUT),
                    guide = "none") +
  scale_y_continuous(breaks = seq(0.85, 1.0, 0.025),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.85, 1.0)) +
  labs(x = NULL, y = "Inner-CV mean AUROC") +
  theme_indelvar() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1,
                                   size = 7.4),
        plot.margin = margin(8, 10, 8, 10))

# Panel B: mtry grid
mtry_tbl <- data.table(
  mtry  = c(6, 12, 18, 37),
  auroc = c(0.9698019, 0.9690233, 0.9682315, 0.9660846),
  sd    = c(0.006, 0.006, 0.006, 0.007)
)
mtry_tbl[, lo := auroc - sd]
mtry_tbl[, hi := pmin(auroc + sd, 1)]
mtry_tbl[, is_best := mtry == 6]
mtry_tbl[, label_pt := sprintf("mtry=%d", mtry)]

pB <- ggplot(mtry_tbl, aes(x = factor(mtry), y = auroc, group = 1)) +
  geom_ribbon(aes(ymin = lo, ymax = hi),
              fill = COL_BEST, alpha = 0.15) +
  geom_line(colour = COL_BEST, linewidth = 0.7) +
  geom_point(aes(fill = is_best), shape = 21, size = 3.2,
             colour = COL_TEXT, stroke = 0.4) +
  geom_text(aes(label = sprintf("%.4f", auroc)),
            vjust = -1.8, size = 2.6, colour = COL_TEXT,
            family = "Inter", fontface = "bold") +
  scale_fill_manual(values = c(`TRUE` = COL_BEST, `FALSE` = "white"),
                    guide = "none") +
  scale_y_continuous(breaks = seq(0.960, 0.975, 0.005),
                     expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.960, 0.975)) +
  labs(x = expression(bold(italic(m)[try])),
       y = "5-fold CV AUROC") +
  theme_indelvar() +
  theme(axis.text.x  = element_text(size = 8.5),
        plot.margin  = margin(8, 10, 8, 10))

pS1 <- pA + pB + plot_layout(widths = c(1.1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 14, family = "Inter"),
    plot.tag.position = c(0, 1.02),
    plot.margin = margin(t = 13, r = 10, b = 8, l = 10)
  )

FIG_NAME <- "figS1_hyperparam"
pdf_path <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
ggsave(pdf_path, pS1, width = 240, height = 110, units = "mm",
       device = "pdf", dpi = 600)
ggsave(png_path, pS1, width = 240, height = 110, units = "mm",
       device = "png", dpi = 600)
cat(sprintf("[done] wrote %s\n[done] wrote %s\n", pdf_path, png_path))
