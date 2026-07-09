# Figure S2: feature importance ranking: Gini and permutation-Shapley attribution.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
COL_TEXT <- "#1F2933"

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

# feature groups: S structure, C conservation, V sequence context, G gene constraint, P protein annotation
grp <- list(
  S = c("af_plddt_mean","af_helix_frac","af_sheet_frac","af_rsa_mean","low_plddt_region",
        "wcn_mean","contact_density","bridging_contacts","ss_element_break",
        "ss_element_frac_affected","helix_register_preserved","hbond_density",
        "hydrophobic_exposure_risk","local_ss_class"),
  C = c("phylop100_mean","gerp_mean"),
  V = c("local_entropy","dist_to_splice","indel_type","indel_size_aa","indel_size_relative"),
  G = c("pli","loeuf","mis_z","cds_length"),
  P = c("in_domain","n_domains","in_repeat","ptm_nearby","ptm_count_near",
        "disulfide_nearby","in_active_site",
        "subcell_Cytoplasm","subcell_Membrane","subcell_Mitochondrion",
        "subcell_Nucleus","subcell_Secreted","subcell_Other","subcell_Unknown"))
feat2grp <- unlist(lapply(names(grp), function(g) setNames(rep(g, length(grp[[g]])), grp[[g]])))

# Darjeeling2 palette
GRP_LAB <- c(S="Protein structure", C="Conservation", V="Sequence context",
             G="Gene constraint", P="Protein annotation")
GRP_COL <- c(S="#046C9A", C="#D69C4E", V="#ABDDDE", G="#4D4D4D", P="#ECCBAE")

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

# ranked bar panel from a (feature, value) table
make_panel <- function(dt, valcol, xlab, ttl, show_legend) {
  d <- copy(dt)
  d[, group := feat2grp[feature]][is.na(group), group := "P"]
  d[, group_lab := factor(GRP_LAB[group], levels = GRP_LAB)]
  d[, pct := 100 * get(valcol) / sum(get(valcol))]
  d[, flab := FEATURE_NAME[feature]]
  d[is.na(flab), flab := gsub("_", " ", feature)]
  d <- d[order(pct)]
  d[, flab := factor(flab, levels = flab)]
  p <- ggplot(d, aes(pct, flab, fill = group_lab)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = sprintf("%.1f", pct)),
              hjust = -0.15, size = 2.2, colour = COL_TEXT) +
    scale_fill_manual(values = setNames(GRP_COL[match(levels(d$group_lab),
                      GRP_LAB)], levels(d$group_lab)), name = "Group") +
    scale_x_continuous(expand = expansion(c(0, 0.14)), name = xlab) +
    labs(y = NULL, title = ttl) +
    theme_indelvar(base_size = 9) +
    theme(axis.text.y = element_text(size = 7, colour = COL_TEXT),
          axis.text.x = element_text(size = 8, colour = COL_TEXT),
          axis.title.x = element_text(size = 9, face = "plain", colour = COL_TEXT),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = "#E5E7EB", linewidth = 0.3))
  if (show_legend) {
    p <- p + theme(legend.position = c(0.74, 0.28),
                   legend.text = element_text(size = 7.5),
                   legend.title = element_text(size = 8, face = "bold"),
                   legend.key.size = unit(3.5, "mm"))
  } else p <- p + theme(legend.position = "none")
  p
}

gini <- fread(file.path(ROOT, "dataset/table_fig_data/feature_importance_retrain.tsv"))
shap <- fread(file.path(ROOT, "dataset/table_fig_data/shap_global_feature.tsv"))

pA <- make_panel(shap, "mean_abs_shap", "Relative SHAP attribution (%)",
                 "Permutation Shapley", show_legend = TRUE)
pB <- make_panel(gini, "gini", "Relative Gini importance (%)",
                 "Mean decrease Gini", show_legend = FALSE)

# SHAP-vs-Gini rank concordance (console log only)
rho_gs <- cor(gini[order(feature)]$gini, shap[order(feature)]$mean_abs_shap, method = "spearman")

figS2 <- (pA | pB) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 12, face = "bold", colour = COL_TEXT))

FIG_NAME <- "figS2_importance"
ggsave(file.path(OUT, paste0(FIG_NAME, ".pdf")), figS2, width = 300, height = 240,
       units = "mm", device = "pdf", dpi = 600)
ggsave(file.path(OUT, paste0(FIG_NAME, ".png")), figS2, width = 300, height = 240,
       units = "mm", device = "png", dpi = 600)
cat(sprintf("[done] wrote figS2_importance (SHAP-Gini Spearman rho=%.2f)\n", rho_gs))
