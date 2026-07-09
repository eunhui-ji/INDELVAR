#!/usr/bin/env Rscript
## Figure S5: test AUROC by AlphaFold pLDDT bin (structural confidence).
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")

COL_BAR  <- "#7C8794"
COL_TEXT <- "#1F2933"
COL_NEUT <- "#9AA3AD"
CATPAL   <- c("#046C9A", "#7294D4", "#D69C4E", "#9AA3AD", "#4D4D4D", "#ECCBAE")

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

pb <- fread(file.path(ROOT, "dataset/table_fig_data/plddt_bin_auroc_test.tsv"))
pb[, plddt_bin := gsub(">=", "≥", plddt_bin, fixed = TRUE)]   # render >= as the combined >= glyph
pb[, plddt_bin := factor(plddt_bin,
  levels = c("Very low (<50, IDR)", "Low (50-70)",
             "Confident (70-90)", "Very high (≥90)"))]
pB <- ggplot(pb, aes(x = plddt_bin, y = auroc, group = 1)) +
  geom_ribbon(aes(ymin = lo95, ymax = hi95), fill = "#C7CDD6", alpha = 0.5) +
  geom_line(colour = COL_BAR, linewidth = 0.8) +
  geom_point(aes(fill = plddt_bin), shape = 21, colour = "white",
             size = 3, stroke = 0.6) +
  scale_fill_manual(values = CATPAL[seq_len(nrow(pb))],
                    guide = "none") +
  geom_text(aes(label = sprintf("%.3f", auroc)), vjust = -1.1, size = 2.5,
            colour = COL_TEXT, family = "Inter", fontface = "bold") +
  geom_text(aes(y = lo95, label = sprintf("n=%d", n)), vjust = 1.8, size = 2.1,
            colour = COL_NEUT, family = "Inter") +
  scale_y_continuous(limits = c(0.5, 1.0), breaks = seq(0.5, 1.0, 0.1),
                     expand = c(0, 0)) +
  labs(x = "Mean pLDDT (deleted residues)",
       y = "Test set AUROC") +
  theme_indelvar() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 8),
        plot.margin = margin(10, 12, 8, 10))

FIG_NAME <- "figS5_subgroup_auroc"
pdf_path <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
cairo_ok <- tryCatch({ .tf <- tempfile(fileext = ".pdf"); cairo_pdf(.tf); dev.off()
                       unlink(.tf); TRUE }, error = function(e) FALSE)
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (cairo_ok) ggsave(pdf_path, pB, width = 120, height = 110, units = "mm",
                     device = cairo_pdf, dpi = 600) else
  message("[figS5] cairo_pdf unavailable (install XQuartz for vector PDF); PNG is the deliverable")
ggsave(png_path, pB, width = 120, height = 110, units = "mm",
       device = ragg::agg_png, dpi = 600)
cat(sprintf("[done] wrote %s\n", png_path))
