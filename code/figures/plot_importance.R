#!/usr/bin/env Rscript
## Figure 4A: SHAP feature-group attribution.
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

COL_TEXT <- "#1F2933"
# Darjeeling2-style group palette, matching Fig 4 in the manuscript
GROUP_COLOURS <- c(Structure = "#046C9A", Conservation = "#D69C4E",
                   SeqContext = "#ABDDDE", Constraint = "#4D4D4D",
                   Annotation = "#ECCBAE")

theme_indelvar <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Inter") +
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      colour = COL_TEXT),
      axis.title       = element_text(size = base_size, face = "plain",
                                      colour = COL_TEXT),
      axis.text        = element_text(size = base_size - 1, colour = COL_TEXT),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
      panel.border     = element_rect(colour = COL_TEXT, fill = NA,
                                      linewidth = 0.4),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

# per-feature mean|SHAP| aggregated to the five feature groups
imp <- fread(file.path(ROOT, "dataset/table_fig_data/shap_global_feature.tsv"))
imp[, pct := 100 * mean_abs_shap / sum(mean_abs_shap)]

GRP_LABEL <- c(Structure = "Protein structure", Conservation = "Conservation",
               SeqContext = "Sequence context", Constraint = "Gene constraint",
               Annotation = "Protein annotation")
grp_sum <- imp[, .(pct = sum(pct)), by = group][order(pct)]
grp_sum[, glab := factor(GRP_LABEL[group], levels = GRP_LABEL[group])]

cat("[panel A] group %:\n")
for (i in seq_len(nrow(grp_sum))) cat(sprintf("   %-18s %.1f%%\n", grp_sum$group[i], grp_sum$pct[i]))

pA <- ggplot(grp_sum, aes(x = pct, y = glab, fill = group)) +
  geom_col(width = 0.66, colour = NA) +
  geom_text(aes(label = sprintf("%.1f%%", pct)),
            hjust = -0.15, size = 2.8, colour = COL_TEXT, family = "Inter",
            fontface = "bold") +
  scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18)),
                     name = "SHAP attribution (% of total)") +
  labs(y = NULL) +
  theme_indelvar() +
  theme(axis.text.y = element_text(size = 8.5, colour = COL_TEXT),
        axis.text.x = element_text(size = 8, colour = COL_TEXT),
        panel.grid.major.y = element_blank(),
        plot.margin  = margin(6, 12, 4, 6))

FIG_NAME <- "fig04a_importance"
pdf_path <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
cairo_ok <- tryCatch({ .tf <- tempfile(fileext = ".pdf"); cairo_pdf(.tf); dev.off()
                       ok <- file.exists(.tf) && file.info(.tf)$size > 0; unlink(.tf); ok },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (cairo_ok) ggsave(pdf_path, pA, width = 140, height = 90, units = "mm",
                     device = cairo_pdf, dpi = 600) else
  message("[fig04a] cairo_pdf unavailable (install XQuartz for vector PDF); PNG is the deliverable")
ggsave(png_path, pA, width = 140, height = 90, units = "mm",
       device = ragg::agg_png, dpi = 600)
cat(sprintf("[done] wrote %s\n", png_path))
