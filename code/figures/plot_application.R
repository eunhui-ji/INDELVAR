# Figure 6: VUS computational evidence assignment, functional validation, deployment.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(arrow)
  library(scales)
  library(ggtext)
})

PROJ <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(PROJ, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Palette
COL_INDELVAR   <- "#FD6467"
COL_PATH     <- "#E47DA3"
COL_BENIGN   <- "#7294D4"
COL_TEXT     <- "#1F2933"
COL_NEUT     <- "#9AA3AD"
# bin palette: graded B/LB blue, neutral grey, graded P/LP red
PEJAVER_COLS <- c(
  BP4_Strong     = "#3E5FA6",   # deep periwinkle
  BP4_Moderate3  = "#5878BD",   # periwinkle (-3)
  BP4_Moderate   = "#7294D4",   # periwinkle
  BP4_Supporting = "#AEC4E6",   # pale periwinkle
  No_Evidence    = "#D9DCE0",   # light grey
  PP3_Supporting = "#F2C0D3",   # pale rose
  PP3_Moderate   = "#E47DA3",   # rose
  PP3_Moderate3  = "#D66A90",   # rose (+3)
  PP3_Strong     = "#C0507A"    # deep rose
)
PEJAVER_ORDER <- names(PEJAVER_COLS)

# Theme
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

# deletion cutoffs (Table 1); del has no PP3_Strong (+4)
cutoffs <- fread(file.path(PROJ,
  "model/evidence_cutoffs.tsv"))
CUT <- setNames(cutoffs[subset == "del"]$cutoff, cutoffs[subset == "del"]$evidence)
T_BP4_SUP <- CUT["BP4_Supporting"]
T_PP3_SUP <- CUT["PP3_Supporting"]

# Panel A: diverging evidence bar
vus <- fread(file.path(PROJ, "dataset/table_fig_data/vus_set_scored.tsv"))
n_vus <- nrow(vus)
n_pp3 <- sum(grepl("^PP3", vus$acmg_call))         # P/LP-candidate (any PP3)
n_bp4 <- sum(grepl("^BP4", vus$acmg_call))         # B/LB-candidate (any BP4)
n_unc <- sum(vus$acmg_call == "Uncertain")
cat(sprintf("[6A] VUS %d: VUS-High %.1f%%, uncertain %.1f%%, VUS-Low %.1f%%; PP3/BP4-assigned %.1f%%\n",
            n_vus, 100*n_pp3/n_vus, 100*n_unc/n_vus, 100*n_bp4/n_vus,
            100*(n_pp3+n_bp4)/n_vus))

# benign tiers left, pathogenic right; darker = stronger
TIER_LV  <- c("BP4_Strong","BP4_Moderate","BP4_Supporting","Uncertain",
              "PP3_Supporting","PP3_Moderate","PP3_Strong")
TIER_COL <- c(BP4_Strong="#3E5FA6", BP4_Moderate="#7294D4", BP4_Supporting="#AEC4E6",
              Uncertain="#D9DCE0", PP3_Supporting="#F2C0D3", PP3_Moderate="#E47DA3",
              PP3_Strong="#C0507A")
TIER_LAB <- c(BP4_Strong="BP4_Strong (−4)", BP4_Moderate="BP4_Moderate (−2/−3)",
              BP4_Supporting="BP4_Supporting (−1)", Uncertain="Indeterminate",
              PP3_Supporting="PP3_Supporting (+1)", PP3_Moderate="PP3_Moderate (+2/+3)",
              PP3_Strong="PP3_Strong (+4)")
tt <- vus[, .(n = .N), by = .(tier = acmg_call)]
tt <- tt[match(TIER_LV, tier)][!is.na(n)]
tt[, tier := factor(tier, levels = TIER_LV)]; setorder(tt, tier)
tt[, pct := 100 * n / sum(n)][, xend := cumsum(pct)][, xstart := xend - pct]
.uc <- tt[tier == "Uncertain", (xstart + xend) / 2]      # centre Indeterminate on 0
tt[, `:=`(xstart = xstart - .uc, xend = xend - .uc)][, xmid := (xstart + xend) / 2]
YB <- 0.55
xL  <- tt[tier=="BP4_Strong", xstart];    xLm <- tt[tier=="BP4_Supporting", xend]
xR  <- tt[tier=="PP3_Strong", xend];      xRm <- tt[tier=="PP3_Supporting", xstart]

pVus <- ggplot(tt) +
  annotate("segment", x = xL+0.4, xend = xLm-0.4, y = YB+0.20, yend = YB+0.20,
           colour = COL_BENIGN, linewidth = 0.4) +
  annotate("segment", x = xRm+0.4, xend = xR-0.4, y = YB+0.20, yend = YB+0.20,
           colour = COL_PATH, linewidth = 0.4) +
  annotate("text", x = (xL+xLm)/2, y = YB+0.46, colour = COL_BENIGN, size = 2.6,
           family = "Inter", fontface = "bold", lineheight = 0.85,
           label = sprintf("BP4-assigned\n%.1f%%", 100*n_bp4/n_vus)) +
  annotate("text", x = (xRm+xR)/2, y = YB+0.46, colour = COL_PATH, size = 2.6,
           family = "Inter", fontface = "bold", lineheight = 0.85,
           label = sprintf("PP3-assigned\n%.1f%%", 100*n_pp3/n_vus)) +
  geom_rect(aes(xmin = xstart, xmax = xend, ymin = -YB, ymax = YB, fill = tier),
            colour = "white", linewidth = 0.5) +
  geom_text(data = tt[pct >= 4],
            aes(x = xmid, y = 0, label = sprintf("%.1f", pct)),
            size = 2.4, colour = COL_TEXT, family = "Inter") +
  geom_segment(data = tt[pct < 4],
               aes(x = xmid, xend = xmid, y = -YB, yend = -YB - 0.13),
               colour = COL_TEXT, linewidth = 0.25) +
  geom_text(data = tt[pct < 4],
            aes(x = xmid, y = -YB - 0.24, label = sprintf("%.1f", pct)),
            size = 2.1, colour = COL_TEXT, family = "Inter") +
  scale_fill_manual(values = TIER_COL, breaks = TIER_LV, labels = TIER_LAB,
                    name = NULL, guide = guide_legend(nrow = 2, byrow = TRUE)) +
  coord_cartesian(ylim = c(-YB-0.35, YB+0.85), clip = "off") +
  labs(x = NULL, y = NULL,
       title = sprintf("ClinVar VUS in-frame indels (n = %s)",
                       formatC(n_vus, big.mark = ",", format = "d"))) +
  theme_void(base_family = "Inter") +
  theme(plot.title = element_text(size = 9.6, face = "bold", hjust = 0.5,
                                  colour = COL_TEXT, margin = margin(b = 8)),
        legend.position = "bottom",
        legend.text = element_text(size = 6.6, colour = COL_TEXT),
        legend.key.size = unit(3, "mm"), legend.margin = margin(t = 6),
        plot.margin = margin(8, 14, 6, 14))

# Panel B: INDELVAR vs DMS readout, three disease indel scans
# DMS lookup: precomputed DB score, deletion assays, TP53 leakage-free (n=323)
pgv <- fread(file.path(PROJ, "dataset/table_fig_data/dms_db_scored.tsv"))
setnames(pgv, "indelvar_db", "indelvar_retrain")
sp  <- fread(file.path(PROJ, "dataset/table_fig_data/dms_db_spearman.tsv"))
E_ASSAYS <- c("P53_HUMAN_Kotler_2018_indels", "A4_HUMAN_Seuma_2022_indels",
              "S22A1_HUMAN_Yee_2023_abundance_indels")
E_GENE <- c(P53_HUMAN_Kotler_2018_indels = "TP53",
            A4_HUMAN_Seuma_2022_indels = "APP",
            S22A1_HUMAN_Yee_2023_abundance_indels = "SLC22A1")
rho_map <- setNames(sp[type == "del"]$rho, sp[type == "del"]$gene)   # keyed by gene
es <- pgv[DMS_id %in% E_ASSAYS & type == "del" & leakage_overlap == FALSE &
          is.finite(indelvar_retrain) & is.finite(DMS_score)]
elab <- es[, .(n = .N), by = DMS_id]
elab[, flab := E_GENE[DMS_id]]
elab[, rlab := sprintf("***ρ* = %.2f**", rho_map[flab])]
es <- merge(es, elab[, .(DMS_id, flab)], by = "DMS_id")
es[,   flab := factor(flab, levels = E_GENE[E_ASSAYS])]
elab[, flab := factor(flab, levels = E_GENE[E_ASSAYS])]
cat(sprintf("[6B] disease-assay functional: %s\n", paste(elab$flab, collapse = " | ")))
pDms <- ggplot(es, aes(x = indelvar_retrain, y = DMS_score)) +
  geom_point(alpha = 0.50, size = 0.95, colour = COL_BENIGN, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              colour = COL_PATH, linewidth = 1.0) +
  geom_richtext(data = elab, aes(x = -Inf, y = Inf, label = rlab),
            hjust = -0.1, vjust = 1.7, size = 2.6, colour = COL_TEXT,
            family = "Inter", fill = NA, label.color = NA,
            label.padding = unit(0, "pt"), inherit.aes = FALSE) +
  facet_wrap(~ flab, nrow = 1, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = "INDELVAR score", y = "DMS functional score") +
  theme_indelvar() +
  theme(strip.text = element_text(size = 8.6, face = "bold"),
        axis.text = element_text(size = 7), axis.title = element_text(size = 8.5))

# Panel C: deployment-scale score distribution, ACMG bands + density
# deletions don't reach PP3_Strong (+4); insertions do
set.seed(42)
allcut <- fread(file.path(PROJ,
  "model/evidence_cutoffs.tsv"))
PC_TYPES <- list(
  deletion  = list(tag = "deletions",           lab = "Deletions",  sub = "del"),
  insertion = list(tag = "insertions",          lab = "Insertions", sub = "ins"))
PC_FLEV <- c("Deletions", "Insertions")
band_names <- c("BP4_Strong","BP4_Moderate3","BP4_Moderate","BP4_Supporting",
                "No_Evidence","PP3_Supporting","PP3_Moderate","PP3_Moderate3","PP3_Strong")
ev_order   <- c("BP4_Strong","BP4_Moderate3","BP4_Moderate","BP4_Supporting",
                "PP3_Supporting","PP3_Moderate","PP3_Moderate3","PP3_Strong")
dbbg_list <- list(); band_list <- list(); vline_list <- list()
for (ty in names(PC_TYPES)) {
  tg <- PC_TYPES[[ty]]; lab <- tg$lab
  cc <- setNames(allcut[subset == tg$sub]$cutoff, allcut[subset == tg$sub]$evidence)
  # needs Zenodo DB parquet, else ACMG bands only
  pqf <- file.path(PROJ, sprintf("data/precomputed_%s.parquet", tg$tag))
  if (file.exists(pqf)) {
    sc <- as.numeric(read_parquet(pqf, col_select = "indelvar_score")$indelvar_score)
    sc <- sc[is.finite(sc)]
    fpp3 <- mean(sc >= cc["PP3_Supporting"]); fbp4 <- mean(sc <= cc["BP4_Supporting"])
    cat(sprintf("[6C] deployment %s %s: PP3>=Sup %.1f%%, BP4<=Sup %.1f%%, no-evidence %.1f%%\n",
                lab, format(length(sc), big.mark = ","),
                100*fpp3, 100*fbp4, 100*(1-fpp3-fbp4)))
    dbbg_list[[ty]] <- data.table(score = sc[sample(length(sc), min(3e5, length(sc)))], facet = lab)
    rm(sc); invisible(gc(verbose = FALSE))
  } else {
    cat(sprintf("[6C] %s: precomputed DB absent (download the Zenodo release to data/precomputed_%s.parquet); deployment density skipped, band landscape shown\n", lab, tg$tag))
  }
  # ACMG bands; drop unreachable PP3_Strong (deletions)
  if (is.na(cc["PP3_Strong"])) {
    edges <- c(0, cc[ev_order[1:7]], 1); bn <- band_names[1:8]
  } else {
    edges <- c(0, cc[ev_order], 1); bn <- band_names
  }
  band_list[[ty]] <- data.table(xmin = head(edges, -1), xmax = tail(edges, -1),
    stratum = factor(bn, levels = PEJAVER_ORDER), facet = lab)
  # BP4/PP3 Supporting = B/LB, P/LP boundaries
  vx <- c(cc["BP4_Supporting"], cc["PP3_Supporting"])
  vd <- data.table(x = as.numeric(vx), cls = c("BP4 cutoff", "PP3 cutoff"),
                   lcol = c("#3E5FA6", "#C0507A"), hj = c(1.15, -0.15),
                   facet = lab)
  vline_list[[ty]] <- vd[!is.na(x)]
}
dbbg_O  <- if (length(dbbg_list)) rbindlist(dbbg_list)[, facet := factor(facet, levels = PC_FLEV)] else data.table(score = numeric(0), facet = factor(character(0), levels = PC_FLEV))
band_dt <- rbindlist(band_list)[,  facet := factor(facet, levels = PC_FLEV)]
vline_dt<- rbindlist(vline_list)[, facet := factor(facet, levels = PC_FLEV)]

# observed test set, per type (HGMD-licensed; optional)
TSPATH <- file.path(PROJ, "dataset/test_data/test_set_scored.tsv")
if (file.exists(TSPATH)) {
  obs_O <- fread(TSPATH)
  obs_O[, facet := factor(ifelse(indel_type == "deletion", "Deletions", "Insertions"),
                          levels = PC_FLEV)]
  obs_O2 <- rbind(
    obs_O[label == "Benign",     .(score = indelvar, facet,
                                   src = "Observed benign (gnomAD-B)")],
    obs_O[label == "Pathogenic", .(score = indelvar, facet,
                                   src = "Observed pathogenic (HGMD-P)")])
} else {
  cat("[6C] test_set_scored.tsv absent (HGMD-licensed; not shipped): observed-variant overlay skipped\n")
  obs_O2 <- data.table(score = numeric(0), facet = factor(character(0), levels = PC_FLEV),
                       src = character(0))
}
SRC_COL <- c("Observed benign (gnomAD-B)" = COL_BENIGN,
             "Observed pathogenic (HGMD-P)" = COL_PATH)
OBS_BRK <- c("Observed benign (gnomAD-B)", "Observed pathogenic (HGMD-P)")
OBS_LAB <- c("Observed B/LB", "Observed P/LP")

pDist <- ggplot() +
  # all-possible deployment density (filled grey)
  geom_density(data = dbbg_O, aes(score, after_stat(density), fill = "ref"),
               colour = "#7F8A95", alpha = 0.45, linewidth = 0.4) +
  scale_fill_manual(name = NULL, values = c(ref = "#A9B2BD"),
                    labels = "All possible indels (proteome-wide)",
                    guide = guide_legend(order = 1)) +
  # observed classes as density lines
  geom_density(data = obs_O2, aes(score, after_stat(density), colour = src),
               fill = NA, linewidth = 0.9, key_glyph = "path") +
  scale_colour_manual(name = NULL, values = SRC_COL, breaks = OBS_BRK,
                      labels = OBS_LAB, guide = guide_legend(order = 2)) +
  geom_vline(data = vline_dt, aes(xintercept = x), linetype = "dashed",
             colour = COL_TEXT, linewidth = 0.25, alpha = 0.55) +
  geom_text(data = vline_dt, aes(x = x, label = cls), y = Inf, vjust = 1.4,
            hjust = vline_dt$hj, size = 2.5, fontface = "bold",
            colour = vline_dt$lcol, family = "Inter") +
  facet_wrap(~ facet, nrow = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                     expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = "INDELVAR score", y = "Probability density") +
  theme_indelvar() +
  theme(strip.text = element_text(size = 8.6, face = "bold"),
        axis.text = element_text(size = 7), axis.title = element_text(size = 8.5),
        panel.spacing.x = unit(5, "mm"),
        legend.position = "bottom", legend.box = "horizontal",
        legend.title = element_blank(), legend.key = element_blank(),
        legend.text = element_text(size = 6.6, colour = COL_TEXT),
        legend.key.size = unit(3.2, "mm"),
        legend.margin = margin(t = 1, b = 0),
        legend.spacing.y = unit(0.5, "mm"))

# assemble: top row A|B, full-width C; free(pVus) for axis-free A
pVus  <- pVus  + labs(tag = "A")   # VUS computational evidence
pDms  <- pDms  + labs(tag = "B")   # DMS functional validation
pDist <- pDist + labs(tag = "C")   # proteome-scale deployment

fig6 <- free(pVus) + pDms + pDist +
  plot_layout(design = "AB\nCC", widths = c(1, 1.35), heights = c(1.05, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(8, 8, 8, 8))) &
  theme(plot.tag = element_text(face = "bold", size = 14, family = "Inter"),
        plot.tag.position = c(0, 1),
        plot.margin = margin(14, 16, 10, 16))

FIG_NAME <- "fig06_application"
pdf_path <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
# probe cairo_pdf (needs XQuartz); skip PDF if unavailable
cairo_ok <- tryCatch({ .tf <- tempfile(fileext = ".pdf"); cairo_pdf(.tf); dev.off()
                       ok <- file.exists(.tf) && file.info(.tf)$size > 0; unlink(.tf); ok },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (cairo_ok) ggsave(pdf_path, fig6, width = 290, height = 200, units = "mm",
                     device = cairo_pdf, dpi = 600) else
  message("[pdf] SKIPPED ", basename(pdf_path),
          " - cairo_pdf unavailable (install XQuartz for vector PDF); 600 dpi PNG is the deliverable")
ggsave(png_path, fig6, width = 290, height = 200, units = "mm",
       device = ragg::agg_png, dpi = 600)
cat(sprintf("[done] wrote %s\n[done] wrote %s\n", pdf_path, png_path))
