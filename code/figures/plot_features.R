# Figure 3: feature characterization: fold-stability burden, feature violins, deleted-residue composition and biophysics.

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

# feature groups: S structure, C conservation, V sequence/indel, G constraint, P protein
FEATURE_GROUP <- c(
  pli = "G", loeuf = "G", mis_z = "G", cds_length = "G",
  af_plddt_mean = "S", af_helix_frac = "S", af_sheet_frac = "S",
  af_rsa_mean = "S", low_plddt_region = "S",
  wcn_mean = "S", contact_density = "S", bridging_contacts = "S",
  ss_element_break = "S", ss_element_frac_affected = "S",
  helix_register_preserved = "S", hbond_density = "S",
  hydrophobic_exposure_risk = "S", local_ss_class = "S",
  in_domain = "P", n_domains = "P", in_repeat = "P",
  ptm_nearby = "P", ptm_count_near = "P", disulfide_nearby = "P",
  in_active_site = "P", indel_size_aa = "V", indel_type = "V",
  indel_size_relative = "V",
  phylop100_mean = "C", gerp_mean = "C",
  local_entropy = "V", dist_to_splice = "V",
  subcell_Cytoplasm = "P", subcell_Membrane = "P",
  subcell_Mitochondrion = "P", subcell_Nucleus = "P",
  subcell_Other = "P", subcell_Secreted = "P", subcell_Unknown = "P"
)

# publication feature labels
FEATURE_NAME <- c(
  phylop100_mean         = "phyloP",
  af_plddt_mean          = "pLDDT",
  af_rsa_mean            = "RSA",
  contact_density        = "Contact density",
  wcn_mean               = "WCN",
  gerp_mean              = "GERP",
  local_entropy          = "Local entropy",
  dist_to_splice         = "Dist. to splice",
  mis_z                  = "Missense Z",
  cds_length             = "CDS length",
  low_plddt_region       = "Low-pLDDT region",
  indel_size_relative    = "Indel size (relative)",
  loeuf                  = "LOEUF",
  pli                    = "pLI",
  hbond_density          = "H-bond density",
  indel_size_aa          = "Indel size (aa)",
  ss_element_frac_affected   = "SS element affected",
  local_ss_class      = "Local SS class",
  subcell_Secreted       = "Secreted",
  subcell_Membrane       = "Membrane",
  hydrophobic_exposure_risk = "Hydrophobic exposure",
  n_domains              = "Domain count",
  af_helix_frac          = "Helix fraction",
  af_sheet_frac          = "Sheet fraction",
  ss_element_break       = "SS element break",
  helix_register_preserved = "Helix register",
  bridging_contacts      = "Bridging contacts",
  in_domain              = "In domain",
  in_repeat              = "In repeat",
  ptm_nearby             = "PTM nearby",
  ptm_count_near         = "PTM count",
  disulfide_nearby       = "Disulfide nearby",
  in_active_site         = "Active site",
  indel_type             = "Indel type",
  subcell_Cytoplasm      = "Cytoplasm",
  subcell_Mitochondrion  = "Mitochondrion",
  subcell_Nucleus        = "Nucleus",
  subcell_Other          = "Subcell: other",
  subcell_Unknown        = "Subcell: unknown"
)
display_name <- function(f) ifelse(is.na(FEATURE_NAME[f]), f, FEATURE_NAME[f])

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

dat <- as.data.table(fread(file.path(PROJ,
  "dataset/train_data/train_set.tsv")))            # 39 features + label per variant
dat[, label := factor(label, levels = c("Benign", "Pathogenic"))]
cat(sprintf("[load] cohort total=%d  P=%d  B=%d\n",
            nrow(dat),
            sum(dat$label == "Pathogenic"),
            sum(dat$label == "Benign")))


# Panel B: mechanisms disrupted (0-4), per-class fraction
casc_dt <- copy(dat)
q75 <- function(x) quantile(as.numeric(x), 0.75, na.rm = TRUE)
casc_dt[, `:=`(
  s1 = as.integer(as.numeric(contact_density)           >= q75(contact_density)),
  s2 = as.integer(as.numeric(ss_element_break)          >= 1),
  s3 = as.integer(as.numeric(hbond_density)             >= q75(hbond_density)),
  s4 = as.integer(as.numeric(hydrophobic_exposure_risk) >= 1))]
for (s in c("s1", "s2", "s3", "s4")) casc_dt[is.na(get(s)), (s) := 0L]
casc_dt[, ncount := s1 + s2 + s3 + s4]
casc_dt[, cls := factor(label, levels = c("Benign", "Pathogenic"))]
cnt <- casc_dt[, .N, by = .(cls, ncount)]
cnt[, pct := N / sum(N), by = cls]           # within-class fraction
nB <- casc_dt[cls == "Benign", .N]
nP <- casc_dt[cls == "Pathogenic", .N]
mB <- casc_dt[cls == "Benign", mean(ncount)]
mP <- casc_dt[cls == "Pathogenic", mean(ncount)]
cat(sprintf("[B] mechanisms disrupted: mean B=%.2f P=%.2f | B 0-hit=%.0f%% P>=1=%.0f%% | nB=%d nP=%d\n",
            mB, mP, 100 * cnt[cls == "Benign" & ncount == 0]$pct,
            100 * casc_dt[cls == "Pathogenic", mean(ncount >= 1)], nB, nP))

pA <- ggplot(cnt, aes(x = factor(ncount), y = pct, fill = cls)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64,
           colour = "white", linewidth = 0.3, alpha = 0.92) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * pct)),
            position = position_dodge(width = 0.72), vjust = -0.45,
            size = 2.3, colour = COL_TEXT, family = "Inter") +
  scale_fill_manual(values = c(Benign = COL_BENIGN, Pathogenic = COL_PATH),
                    labels = c(Benign     = sprintf("B/LB (n=%s)", format(nB, big.mark = ",")),
                               Pathogenic = sprintf("P/LP (n=%s)", format(nP, big.mark = ","))),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 0.95), breaks = seq(0, 0.8, 0.2),
                     labels = scales::percent_format(accuracy = 1, suffix = ""),
                     expand = c(0, 0)) +
  labs(x = "Number of mechanisms disrupted",
       y = "Fraction within class (%)") +
  theme_indelvar() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = c(0.99, 0.96), legend.justification = c(1, 1),
        legend.background = element_rect(fill = "white", colour = NA),
        legend.key.size = unit(0.34, "cm"),
        legend.text = element_text(size = 7.5))

# Cohen's d for every continuous-looking feature
cohen_d <- function(x, g) {
  if (is.factor(x))    x <- as.integer(x) - 1L
  if (is.character(x)) x <- as.integer(as.factor(x)) - 1L
  if (is.logical(x))   x <- as.integer(x)
  x <- as.numeric(x)
  x1 <- x[g == "Pathogenic"]; x0 <- x[g == "Benign"]
  x1 <- x1[is.finite(x1)];     x0 <- x0[is.finite(x0)]
  if (length(x1) < 5 || length(x0) < 5) return(NA_real_)
  v1 <- var(x1); v0 <- var(x0)
  if (!is.finite(v1) || !is.finite(v0) || v1 + v0 == 0) return(NA_real_)
  sp <- sqrt(((length(x1) - 1) * v1 + (length(x0) - 1) * v0) /
               (length(x1) + length(x0) - 2))
  (mean(x1) - mean(x0)) / sp
}

feat_vec <- names(FEATURE_GROUP)
feat_vec <- intersect(feat_vec, names(dat))

# indel_type -> numeric for Cohen's d
if ("indel_type" %in% names(dat) && is.character(dat$indel_type)) {
  dat[, indel_type := as.integer(indel_type == "insertion")]
}

ds <- vapply(feat_vec,
             function(f) cohen_d(dat[[f]], dat$label),
             numeric(1))
cohen_tbl <- data.table(feature = feat_vec, d = ds)
cohen_tbl[, group := FEATURE_GROUP[feature]]
cohen_tbl <- cohen_tbl[is.finite(d)]
cat("[features] Cohen's d ranking (top 15):\n")
print(cohen_tbl[order(-abs(d))][1:15])

# top feature per group -> violins below

# Panel B: violins, six leading continuous features
panelB_feats <- c("af_plddt_mean", "af_rsa_mean", "contact_density",
                  "hbond_density", "phylop100_mean", "gerp_mean")
top_per_group <- cohen_tbl[match(panelB_feats, feature)][!is.na(feature)]
cat("[B] top feature per group:\n"); print(top_per_group)

mk_violin <- function(feature_name, group_letter, d_val) {
  v <- dat[[feature_name]]
  if (is.character(v)) v <- as.numeric(as.factor(v))
  vd <- data.table(label = dat$label, val = as.numeric(v))[is.finite(val)]
  ggplot(vd, aes(x = label, y = val, fill = label)) +
    geom_violin(colour = NA, alpha = 0.55, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.15, outlier.shape = NA,
                 colour = COL_TEXT, fill = "white", linewidth = 0.3) +
    geom_richtext(data = data.frame(x = -Inf, y = Inf,
                    label = sprintf("***d* = %+.2f**", d_val)),
                  aes(x = x, y = y, label = label), vjust = 1.5, hjust = -0.1,
                  size = 3.3, colour = COL_TEXT, family = "Inter", fill = NA,
                  label.color = NA, label.padding = unit(0, "pt"),
                  inherit.aes = FALSE) +
    scale_fill_manual(values = c(Benign = COL_BENIGN, Pathogenic = COL_PATH),
                      guide = "none") +
    scale_x_discrete(labels = c(Benign = "B/LB", Pathogenic = "P/LP")) +
    labs(x = NULL, y = NULL, title = display_name(feature_name)) +
    theme_indelvar() +
    theme(axis.text.x  = element_text(size = 7),
          plot.title = element_text(size = 9.6, face = "bold", hjust = 0.5,
                                    colour = COL_TEXT,
                                    margin = margin(0, 0, 3, 0)))
}

violin_list <- mapply(mk_violin,
                       top_per_group$feature,
                       as.character(top_per_group$group),
                       top_per_group$d,
                       SIMPLIFY = FALSE)
pB <- wrap_plots(violin_list, nrow = 1, ncol = 6)

# Panel C: deleted-residue AA composition
comp <- fread(file.path(PROJ,
  "dataset/table_fig_data/deleted_residue_composition.tsv"))
comp[, aa := factor(aa, levels = comp[order(log2FC)]$aa)]
pComp <- ggplot(comp, aes(x = log2FC, y = aa, colour = enrich)) +
  geom_vline(xintercept = 0, colour = COL_NEUT, linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2FC, yend = aa), linewidth = 0.7) +
  geom_point(size = 2.0) +
  scale_colour_manual(values = c(Pathogenic = COL_PATH, Benign = COL_BENIGN),
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(x = expression(log[2] * "FC (P/LP versus B/LB)"),
       y = "Deleted residue") +
  theme_indelvar() +
  theme(axis.text.y = element_text(size = 6.4),
        panel.grid.major.y = element_blank())

# Panel D: bulkiness (Zimmerman), disorder propensity (TOP-IDP)
prop <- fread(file.path(PROJ,
  "dataset/table_fig_data/deleted_residue_properties.tsv"))
prop[, label := factor(label, levels = c("Benign", "Pathogenic"))]
pv_b <- wilcox.test(bulk ~ label, prop)$p.value
pv_t <- wilcox.test(topidp ~ label, prop)$p.value
cat(sprintf("[F] deleted-residue bulkiness p=%s ; disorder-propensity p=%s\n",
            format.pval(pv_b), format.pval(pv_t)))
plong <- melt(prop[, .(label, Bulkiness = bulk, `Disorder propensity` = topidp)],
              id.vars = "label", variable.name = "scale", value.name = "val")
plab <- data.table(
  scale = factor(c("Bulkiness", "Disorder propensity"),
                 levels = levels(plong$scale)),
  lab = c("***P* < 0.001**", "***P* < 0.001**"))
pBio <- ggplot(plong, aes(x = label, y = val, fill = label)) +
  geom_violin(colour = NA, alpha = 0.55, scale = "width") +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_richtext(data = plab, aes(x = -Inf, y = Inf, label = lab), vjust = 1.5,
            hjust = -0.1, size = 3.3, colour = COL_TEXT, family = "Inter",
            fill = NA, label.color = NA, label.padding = unit(0, "pt"),
            inherit.aes = FALSE) +
  facet_wrap(~ scale, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(Benign = COL_BENIGN, Pathogenic = COL_PATH),
                    guide = "none") +
  scale_x_discrete(labels = c(Benign = "B/LB", Pathogenic = "P/LP")) +
  labs(x = NULL, y = "Mean over deleted residues") +
  theme_indelvar() +
  theme(axis.text.x = element_text(size = 7),
        strip.text = element_text(size = 8.6, face = "bold"))

# assemble: A violins, B fold destabilization, C-D deleted-residue signature
pA <- pA + labs(tag = "B", title = "Fold destabilization")
pB <- pB + plot_annotation(tag_levels = NULL) +
  plot_layout(guides = "collect")
pB_wrap <- wrap_elements(full = pB) + labs(tag = "A") +
  theme(plot.tag = element_text(face = "bold", size = 14,
                                family = "Inter"),
        plot.tag.position = c(0, 1),
        plot.margin = margin(14, 14, 10, 14))
pComp <- pComp + labs(tag = "C", title = "Deleted residue composition")
pBio  <- pBio  + labs(tag = "D", title = "Deleted residue biophysics")

row_cd <- (pComp + pBio) + plot_layout(widths = c(1, 1.05))

fig3 <- pB_wrap / pA / row_cd +
  plot_layout(heights = c(1.0, 0.72, 0.95)) +
  plot_annotation(
    theme = theme(plot.margin = margin(8, 8, 8, 8))
  ) &
  theme(
    plot.tag = element_text(face = "bold", size = 14, family = "Inter"),
    plot.tag.position = c(0, 1),
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0.5,
                              colour = COL_TEXT, margin = margin(0, 0, 4, 0)),
    plot.margin = margin(13, 12, 9, 12)
  )

FIG_NAME  <- "fig03_features"
pdf_path  <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
ggsave(pdf_path,  fig3, width = 297, height = 250, units = "mm",
       device = "pdf",  dpi = 600)
ggsave(png_path, fig3, width = 297, height = 250, units = "mm",
       device = "png", dpi = 600)
cat(sprintf("[done] wrote %s\n[done] wrote %s\n", pdf_path, png_path))
