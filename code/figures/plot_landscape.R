# Figure 2: genome-wide indel landscape: pLDDT x RSA enrichment, size x SS heatmap, contact density, domain location.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(ggtext)
})

PROJ <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(PROJ, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Palette
COL_INDELVAR <- "#FD6467"
COL_PATH   <- "#E47DA3"
COL_BENIGN <- "#7294D4"
COL_TEXT   <- "#1F2933"
COL_NEUT   <- "#9AA3AD"
OR_LIM     <- c(-3.5, 5)   # log2-OR range, panels A & B

theme_indelvar <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Inter") +
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold",
                                      colour = COL_TEXT, hjust = 0.5,
                                      margin = margin(b = 5)),
      axis.title       = element_text(size = base_size, face = "plain",
                                      colour = COL_TEXT),
      axis.title.y     = element_text(margin = margin(r = 4)),
      axis.title.x     = element_text(margin = margin(t = 4)),
      axis.text        = element_text(size = base_size - 1, colour = COL_TEXT),
      legend.title     = element_markdown(size = base_size - 1, face = "bold",
                                          colour = COL_TEXT, lineheight = 1.2),
      legend.text      = element_text(size = base_size - 1, colour = COL_TEXT),
      legend.key.size  = unit(3, "mm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(colour = "#1F2933", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA)
    )
}

dat <- as.data.table(fread(file.path(PROJ,
  "dataset/train_data/train_set.tsv")))
dat[, label := factor(label, levels = c("Benign", "Pathogenic"))]
cat(sprintf("[load] cohort total=%d  P=%d  B=%d  genes=%d\n",
            nrow(dat),
            sum(dat$label == "Pathogenic"),
            sum(dat$label == "Benign"),
            length(unique(dat$gene_symbol))))

# Panel A: pLDDT x RSA bins
panA_dat <- copy(dat)[is.finite(af_plddt_mean) & is.finite(af_rsa_mean)]
panA_dat[, plddt_bin := cut(af_plddt_mean,
                            breaks = c(-Inf, 50, 70, 90, Inf),
                            labels = c("<50", "50-70", "70-90", "≥90"))]
panA_dat[, rsa_bin   := cut(af_rsa_mean,
                            breaks = c(-Inf, 0.2, 0.5, Inf),
                            labels = c("Buried",
                                       "Intermediate",
                                       "Exposed"))]

panA_summ <- panA_dat[, .(n = .N,
                           n_p = sum(label == "Pathogenic"),
                           n_b = sum(label == "Benign")),
                       by = .(plddt_bin, rsa_bin)]
panA_summ[, frac_p := n_p / n]
# marginal pathogenic fraction
P_global <- mean(panA_dat$label == "Pathogenic")
# odds ratio per cell, log2-scaled
glob_or <- sum(dat$label == "Pathogenic") / sum(dat$label == "Benign")
panA_summ[, or := pmax((n_p + 0.5) / (n_b + 0.5) / glob_or, 0.1)]
panA_summ[, log_or := log2(or)]
panA_summ <- panA_summ[!is.na(plddt_bin) & !is.na(rsa_bin)]
# n<10 bins -> grey "n<10" tile
panA_summ[, lab    := ifelse(n >= 10, sprintf("%.1f", or), "n<10")]
panA_summ[, log_or := ifelse(n >= 10, log_or, NA_real_)]
cat("[A] pLDDT x RSA odds-ratio table:\n"); print(panA_summ[order(-or)])

pA <- ggplot(panA_summ,
             aes(x = rsa_bin, y = plddt_bin, fill = log_or)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab),
            size = 2.6, colour = COL_TEXT, family = "Inter",
            fontface = "bold") +
  # diverging fill, mid = OR 1
  scale_fill_gradient2(low = "#4A6FE3", mid = "#E2E2E2",
                       high = "#D33F6A", midpoint = 0,
                       limits = OR_LIM, oob = scales::squish,
                       na.value = "#9AA3AD",
                       name = "log<sub>2</sub> odds ratio<br>(P/LP versus B/LB)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "AlphaFold pLDDT",
       title = "Relative solvent accessibility") +
  theme_indelvar() +
  theme(legend.position = "right",
        panel.grid = element_blank())

# Panel B: indel size x SS class
panB_dat <- copy(dat)[is.finite(indel_size_aa)]
panB_dat[, size_bin := cut(indel_size_aa,
                            breaks = c(0, 1, 2, 3, 4, 5, 10),
                            labels = c("1", "2", "3", "4", "5", "6-10"),
                            right = TRUE, include.lowest = TRUE)]
# local_ss_class: 1 helix, 2 strand; low_plddt_region flags disorder
panB_dat[, ss_class := fcase(
  local_ss_class   == 1,                           "Helix",
  local_ss_class   == 2,                           "Strand",
  low_plddt_region == 1,                           "IDR / coil",
  default                                        = "Loop"
)]
panB_summ <- panB_dat[, .(n = .N,
                           n_p = sum(label == "Pathogenic"),
                           n_b = sum(label == "Benign")),
                       by = .(size_bin, ss_class)]
panB_summ <- panB_summ[n >= 10 & !is.na(size_bin)]

# odds ratio = per-cell / global
glob_or <- sum(dat$label == "Pathogenic") / sum(dat$label == "Benign")
panB_summ[, or := pmax((n_p + 0.5) / (n_b + 0.5) / glob_or, 0.1)]
panB_summ[, log_or := log2(or)]
panB_summ[, ss_class := factor(ss_class,
  levels = c("Helix", "Strand", "Loop", "IDR / coil"))]
panB_summ[, size_bin := factor(size_bin,
  levels = c("1", "2", "3", "4", "5", "6-10"))]

pB <- ggplot(panB_summ,
             aes(x = ss_class, y = size_bin, fill = log_or)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", or)),
            size = 2.5, colour = COL_TEXT, family = "Inter",
            fontface = "bold") +
  scale_fill_gradient2(low = "#4A6FE3", mid = "#E2E2E2",   # Blue-Red2
                       high = "#D33F6A", midpoint = 0,
                       limits = OR_LIM, oob = scales::squish,
                       name = "log<sub>2</sub> odds ratio<br>(P/LP versus B/LB)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "Indel length",
       title = "Secondary structure") +
  theme_indelvar() +
  theme(legend.position = "right",
        panel.grid = element_blank())

# Panel C: contact density step
panC_dat <- copy(dat)[is.finite(contact_density)]
# contact density deciles
cd_qs <- unique(quantile(panC_dat$contact_density,
                         probs = seq(0, 1, 0.1), na.rm = TRUE))
panC_dat[, cd_decile := cut(contact_density,
                             breaks = cd_qs,
                             include.lowest = TRUE, labels = FALSE)]
panC_summ <- panC_dat[!is.na(cd_decile),
                       .(n = .N,
                         frac_p = mean(label == "Pathogenic"),
                         cd_mid = median(contact_density)),
                       by = cd_decile][order(cd_decile)]
# Wilson 95% CI for pathogenic fraction
wilson_ci <- function(k, n) {
  if (n == 0) return(c(NA, NA))
  z <- 1.96
  ph <- k / n
  denom <- 1 + z^2 / n
  cen <- (ph + z^2 / (2 * n)) / denom
  hw  <- z * sqrt(ph * (1 - ph) / n + z^2 / (4 * n^2)) / denom
  c(cen - hw, cen + hw)
}
panC_summ[, c("lo", "hi") := as.data.table(t(mapply(
  function(k, n) wilson_ci(k, n),
  round(frac_p * n), n)))]

pC <- ggplot(panC_summ, aes(x = cd_decile, y = frac_p)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#F8BBD0", alpha = 0.5) +
  geom_step(colour = COL_PATH, linewidth = 0.9, direction = "mid") +
  geom_point(colour = COL_PATH, size = 1.6) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * frac_p)),
            vjust = -0.7, size = 2.4, colour = COL_TEXT,
            family = "Inter", fontface = "bold") +
  scale_x_continuous(breaks = 1:10, expand = c(0.02, 0.02)) +
  scale_y_continuous(limits = c(0, 1.05),
                     labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = c(0, 0)) +
  labs(x = "Decile", y = "Observed P/LP fraction (%)",
       title = "Contact density") +
  theme_indelvar()

# Panel D: domain location
panD_dat <- copy(dat)[, dom_cat := fcase(
  in_domain == 1L & in_repeat == 1L, "Domain + repeat",
  in_domain == 1L,                   "Domain",
  in_repeat == 1L,                   "Repeat",
  default                          = "No annotation"
)]
panD_summ <- panD_dat[, .N, by = .(dom_cat, label)]
panD_summ[, frac := N / sum(N), by = label]
panD_summ[, dom_cat := factor(dom_cat,
  levels = c("Domain", "Domain + repeat", "Repeat", "No annotation"))]

pD <- ggplot(panD_summ,
             aes(x = label, y = frac, fill = dom_cat)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(frac >= 0.05,
                               sprintf("%.0f%%", 100 * frac), "")),
            position = position_stack(vjust = 0.5),
            size = 2.5, colour = COL_TEXT, family = "Inter",
            fontface = "bold") +
  # Darjeeling2 palette
  scale_fill_manual(values = c(`Domain`           = "#ECCBAE",
                                `Domain + repeat`  = "#046C9A",
                                `Repeat`           = "#D69C4E",
                                `No annotation`    = "#C8CDD3"),
                    name = NULL) +
  scale_x_discrete(labels = c(Benign = "B/LB", Pathogenic = "P/LP")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = c(0, 0), limits = c(0, 1.001)) +
  labs(x = NULL, y = "Within-class proportion (%)",
       title = "Domain overlap") +
  theme_indelvar() +
  theme(legend.position = "right",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7),
        panel.grid.major.x = element_blank())

# assemble: 2x2 grid
pA <- pA + labs(tag = "A")
pB <- pB + labs(tag = "B")
pC <- pC + labs(tag = "C")
pD <- pD + labs(tag = "D")

# shared legend, A+B
row1 <- (pA | pB) + plot_layout(guides = "collect")

body <- (pC | pD)

fig2 <- row1 / body +
  plot_layout(heights = c(1.3, 1.2)) +
  plot_annotation(
    theme = theme(plot.margin = margin(8, 8, 8, 8))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14, family = "Inter"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(16, 15, 14, 15)
  )

FIG_NAME  <- "fig02_landscape"
pdf_path  <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
# probe cairo_pdf (needs XQuartz); skip PDF if unavailable
cairo_ok <- tryCatch({ .tf <- tempfile(fileext = ".pdf"); cairo_pdf(.tf); dev.off()
                       ok <- file.exists(.tf) && file.info(.tf)$size > 0; unlink(.tf); ok },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (cairo_ok) ggsave(pdf_path,  fig2, width = 297, height = 180, units = "mm",
                     device = cairo_pdf,  dpi = 600) else
  message("[pdf] SKIPPED ", basename(pdf_path),
          " - cairo_pdf unavailable (install XQuartz for vector PDF); 600 dpi PNG is the deliverable")
ggsave(png_path, fig2, width = 297, height = 180, units = "mm",
       device = "png", dpi = 600)
cat(sprintf("[done] wrote %s\n[done] wrote %s\n", pdf_path, png_path))
