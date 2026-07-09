# Figure 5: local-posterior calibration of INDELVAR: PP3 pathogenic, BP4 benign.
suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales); library(patchwork)
})
PROJ <- Sys.getenv("INDELVAR_ROOT", "."); setwd(PROJ)
set.seed(42)
E <- 0.05; NBOOT <- 1000; PC <- 0.5; GRID <- seq(0, 1, by = 0.005)

# PP3/BP4 threshold colors: strong red/blue
COL_INDELVAR <- "#FD6467"; COL_PATH <- "#D62728"; COL_BENIGN <- "#2166AC"
COL_TEXT <- "#1F2933"
theme_indelvar <- function(base_size = 9) {
  theme_minimal(base_size = base_size, base_family = "Inter") +
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold", colour = COL_TEXT),
      axis.title       = element_text(size = base_size, face = "plain", colour = COL_TEXT),
      axis.text        = element_text(size = base_size - 1, colour = COL_TEXT),
      strip.text       = element_text(size = base_size, face = "bold", colour = COL_TEXT),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#E5E7EB", linewidth = 0.3),
      panel.spacing    = unit(1.1, "lines"),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(colour = "#1F2933", fill = NA, linewidth = 0.4),
      panel.background = element_rect(fill = "white", colour = NA))
}
r3 <- function(z) sprintf("%.3f", floor(z * 1000 + 0.5 + 1e-6) / 1000)

# adaptive-window local posterior: pos/(pos + w*neg); one-sided bound = bootstrap 5th pct
minpoints <- 100; gft <- 0.03; increment <- 0.001; HGRID <- seq(0, 1, by = increment)
GRID <- seq(0, 1, by = 0.005); NBOOT <- 500
source(file.path(PROJ, "code", "calib_params.R"))   # alpha/OP
PARAM <- CALIB_PARAM

oof  <- fread("dataset/train_data/train_oof_predictions.tsv")
setnames(oof, grep("indelvar|score|prob", names(oof), value = TRUE)[1], "sc")
full <- fread("dataset/train_data/train_set.tsv")               # variant_id + indel_type
oof  <- merge(oof, full[, .(variant_id, indel_type)], by = "variant_id", all.x = TRUE)
abd  <- fread("dataset/calibration/gnomad_benign_scored.tsv")        # population-benign reference
abd[, gt := ifelse(indel_type == "deletion", "deletion", "insertion")]

local_post <- function(x, y, g, thrs, w, maxt, mint) {
  Ng <- length(g); xp <- sort(x[y == 1]); xn <- sort(x[y == 0]); gs <- sort(g)
  cnt <- function(v, lo, hi) findInterval(hi, v) - findInterval(lo - 1e-12, v)
  vapply(thrs, function(t) {
    ok <- function(h) { lo <- t - h; hi <- t + h
      cc <- if (hi > maxt) (maxt - lo) / (hi - lo) else if (lo < mint) (hi - mint) / (hi - lo) else 1
      if (cc <= 0) cc <- 1e-9
      (cnt(xp, lo, hi) + cnt(xn, lo, hi) >= cc * minpoints) && (cnt(gs, lo, hi) >= gft * cc * Ng) }
    lo <- 1L; hi <- length(HGRID); if (!ok(HGRID[hi])) return(NA_real_)
    while (lo < hi) { m <- (lo + hi) %/% 2L; if (ok(HGRID[m])) hi <- m else lo <- m + 1L }
    h <- HGRID[lo]; pos <- cnt(xp, t - h, t + h); neg <- cnt(xn, t - h, t + h)
    pos / (pos + w * neg) }, numeric(1))
}
mk <- function(tt) {
  alpha <- PARAM[[tt]]$alpha
  d <- oof[!is.na(sc) & !is.na(y) & indel_type == tt]; x <- d$sc; yy <- d$y
  g <- abd[gt == tt & !is.na(indelvar_score), indelvar_score]
  w <- (1 - alpha) * sum(yy == 1) / (sum(yy == 0) * alpha)
  maxt <- max(c(x, g)); mint <- min(c(x, g))
  pt <- local_post(x, yy, g, GRID, w, maxt, mint)
  BM <- matrix(NA_real_, NBOOT, length(GRID))
  for (b in 1:NBOOT) {
    qx <- sample.int(length(x), replace = TRUE); qg <- sample.int(length(g), replace = TRUE)
    BM[b, ] <- local_post(x[qx], yy[qx], g[qg], GRID, w, maxt, mint)
  }
  data.table(type = tt, score = GRID, pst = pt,
             lo = apply(BM, 2, quantile, .05, na.rm = TRUE),
             hi = apply(BM, 2, quantile, .95, na.rm = TRUE))
}
lrc <- rbindlist(list(mk("deletion"), mk("insertion")))
lrc[, xpos := score]
flev <- c("Deletions", "Insertions")
lrc[, facet := factor(ifelse(type == "deletion", flev[1], flev[2]), levels = flev)]

# Table 1 cutoffs (del no +4); LR target OP^(pts/8)
ct <- fread("model/evidence_cutoffs.tsv")   # subset/evidence/points/cutoff
ct[, type := ifelse(subset == "del", "deletion", "insertion")]
ct[, facet := factor(ifelse(type == "deletion", flev[1], flev[2]), levels = flev)]
ct[, side := ifelse(grepl("PP3", evidence), "PP3", "BP4")]
ct[, col := ifelse(side == "BP4", COL_BENIGN, COL_PATH)]
LAB_MAP <- c(
  PP3_Supporting = "+1", PP3_Moderate = "+2",
  PP3_Moderate3  = "+3", PP3_Strong   = "+4",
  BP4_Supporting = "−1", BP4_Moderate = "−2",
  BP4_Moderate3  = "−3", BP4_Strong   = "−4")
ct[, lab := LAB_MAP[evidence]]
ct[, alpha := vapply(seq_len(.N), function(i) PARAM[[type[i]]]$alpha, numeric(1))]
ct[, target := vapply(seq_len(.N), function(i) PARAM[[type[i]]]$OP ^ (points[i] / 8), numeric(1))]
ct[, ypost := target * alpha / (target * alpha + (1 - alpha))]   # LR target -> posterior threshold
LTY_MAP <- c(PP3_Supporting = "dotted", PP3_Moderate = "dashed",
             PP3_Moderate3 = "dotdash", PP3_Strong = "solid",
             BP4_Supporting = "dotted", BP4_Moderate = "dashed",
             BP4_Moderate3 = "dotdash", BP4_Strong = "solid")
LR_T <- ct[, .(facet, side, col, evidence, y = ypost,
               lab = LAB_MAP[evidence], lty = LTY_MAP[evidence])]
# plot height: PP3 = ypost, BP4 = 1 - ypost
TH <- copy(LR_T); TH[, yplot := ifelse(side == "PP3", y, 1 - y)]
# curve data: PP3 = P(path), BP4 = 1 - P(path)
pp3d <- lrc[, .(facet, xpos, val = pst,     bnd = lo)]
bp4d <- lrc[, .(facet, xpos, val = 1 - pst, bnd = 1 - hi)]

# black = point estimate, grey = 95% bound; line style = tier strength
# ── Panel A: PP3, pathogenic ──
pPP3 <- ggplot(pp3d, aes(xpos)) +
  geom_hline(data = TH[side == "PP3"],
             aes(yintercept = yplot, colour = col, linetype = lty), linewidth = .34) +
  scale_linetype_identity() +
  geom_vline(data = ct[side == "PP3" & !is.na(cutoff)],
             aes(xintercept = cutoff, colour = col),
             linetype = "dotted", linewidth = .3, alpha = .6) +
  geom_line(aes(y = bnd), colour = "grey62", linewidth = .5) +
  geom_line(aes(y = val), colour = "black",  linewidth = .65) +
  scale_colour_identity() +
  geom_text(data = ct[side == "PP3" & !is.na(cutoff)],
            aes(x = cutoff, y = Inf, label = r3(cutoff), colour = col),
            angle = 90, hjust = 1.05, vjust = -.3, size = 2.2, fontface = "bold",
            inherit.aes = FALSE, family = "Inter") +
  facet_wrap(~facet) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1),
                     expand = expansion(mult = .015)) +
  scale_y_continuous(breaks = seq(0, 1, .25), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 1), clip = "off") +
  labs(x = NULL, y = "Posterior P(pathogenic)") +
  theme_indelvar(10) +
  theme(axis.title.y = element_text(size = 9, colour = COL_TEXT, margin = margin(r = 1)),
        strip.text = element_text(size = 11, face = "bold", colour = COL_TEXT))

# ── Panel B: BP4, benign ──
BP4_LO <- 0.95
pBP4 <- ggplot(bp4d, aes(xpos)) +
  geom_hline(data = TH[side == "BP4"],
             aes(yintercept = yplot, colour = col, linetype = lty), linewidth = .34) +
  scale_linetype_identity() +
  geom_vline(data = ct[side == "BP4" & !is.na(cutoff)],
             aes(xintercept = cutoff, colour = col),
             linetype = "dotted", linewidth = .3, alpha = .6) +
  geom_line(data = bp4d[bnd >= BP4_LO], aes(y = bnd), colour = "grey62", linewidth = .5) +
  geom_line(data = bp4d[val >= BP4_LO], aes(y = val), colour = "black",  linewidth = .65) +
  scale_colour_identity() +
  geom_text(data = ct[side == "BP4" & !is.na(cutoff)],
            aes(x = cutoff, y = Inf, label = r3(cutoff), colour = col),
            angle = 90, hjust = 1.05, vjust = -.3, size = 2.2, fontface = "bold",
            inherit.aes = FALSE, family = "Inter") +
  facet_wrap(~facet) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1),
                     expand = expansion(mult = .015)) +
  scale_y_continuous(breaks = seq(BP4_LO, 1, .01), expand = c(0, 0)) +
  coord_cartesian(ylim = c(BP4_LO, 1), clip = "off") +
  labs(x = "INDELVAR score", y = "Posterior P(benign)") +
  theme_indelvar(10) +
  theme(axis.title.y = element_text(size = 9, colour = COL_TEXT, margin = margin(r = 1)),
        strip.text = element_blank())

# ── shared bottom legend: point/bound curves + tier line styles ──
LG <- function(y, x0, xt, col, lty, lab)
  data.table(y = y, x0 = x0, x1 = x0 + 0.085, xt = xt, col = col, lty = lty, lab = lab)
LX <- 0.20; RX <- 0.50                             # left/right column starts
leg <- rbindlist(list(
  LG(5, LX, LX + 0.10, "black",    "solid",   "Point estimate"),
  LG(5, RX, RX + 0.10, "grey62",   "solid",   "One-sided confidence bound"),
  LG(4, LX, LX + 0.10, COL_PATH,   "dotted",  "+1 (PP3_Supporting)"),
  LG(4, RX, RX + 0.10, COL_BENIGN, "dotted",  "−1 (BP4_Supporting)"),
  LG(3, LX, LX + 0.10, COL_PATH,   "dashed",  "+2 (PP3_Moderate)"),
  LG(3, RX, RX + 0.10, COL_BENIGN, "dashed",  "−2 (BP4_Moderate)"),
  LG(2, LX, LX + 0.10, COL_PATH,   "dotdash", "+3 (PP3_Moderate)"),
  LG(2, RX, RX + 0.10, COL_BENIGN, "dotdash", "−3 (BP4_Moderate)"),
  LG(1, LX, LX + 0.10, COL_PATH,   "solid",   "+4 (PP3_Strong)"),
  LG(1, RX, RX + 0.10, COL_BENIGN, "solid",   "−4 (BP4_Strong)")))
pLeg <- ggplot(leg) +
  geom_segment(aes(x = x0, xend = x1, y = y, yend = y, colour = col, linetype = lty),
               linewidth = 0.75) +
  geom_text(aes(x = xt, y = y, label = lab), hjust = 0, size = 2.8,
            family = "Inter", colour = COL_TEXT) +
  scale_colour_identity() + scale_linetype_identity() +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 5.5), expand = c(0, 0)) +
  theme_void() +
  theme(plot.tag = element_blank(), plot.margin = margin(2, 6, 2, 6))

fig <- pPP3 / pBP4 / pLeg + plot_layout(heights = c(1, 1, 0.34)) +
  plot_annotation(tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 13)))
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
ggsave("figures/fig05_lr_curve.png", fig, width = 9.5, height = 7.0, dpi = 600)
ggsave("figures/fig05_lr_curve.pdf", fig, width = 9.5, height = 7.0)
cat("wrote figures/fig05_lr_curve.{png,pdf}",
    "(A PP3 pathogenic | B BP4 benign posterior + shared tier legend)\n")
