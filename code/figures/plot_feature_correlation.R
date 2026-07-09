# Figure S3: Pearson correlation matrix of INDELVAR features on the training set.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

ROOT <- Sys.getenv("INDELVAR_ROOT", ".")
OUT  <- file.path(ROOT, "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

COL_INDELVAR <- "#FD6467"
COL_PATH   <- "#E47DA3"
COL_BENIGN <- "#7294D4"
COL_TEXT   <- "#1F2933"

# feature groups (Darjeeling2 palette)
GROUP_COLOURS <- c(S = "#046C9A", C = "#B0791E", V = "#3A8E98",
                   G = "#4D4D4D", P = "#AD7C4A")

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

FEATURE_GROUP <- c(
  # G: gene-level constraint (4)
  pli="G", loeuf="G", mis_z="G", cds_length="G",
  # S: protein structure, raw + derived (14)
  af_plddt_mean="S", af_helix_frac="S", af_sheet_frac="S",
  af_rsa_mean="S", low_plddt_region="S",
  wcn_mean="S", contact_density="S", bridging_contacts="S",
  ss_element_break="S", ss_element_frac_affected="S",
  helix_register_preserved="S", hbond_density="S",
  hydrophobic_exposure_risk="S", local_ss_class="S",
  # P: protein annotation: domains/PTMs/sites + subcellular localization (14)
  in_domain="P", n_domains="P", in_repeat="P",
  ptm_nearby="P", ptm_count_near="P", disulfide_nearby="P",
  in_active_site="P",
  subcell_Cytoplasm="P", subcell_Membrane="P",
  subcell_Mitochondrion="P", subcell_Nucleus="P",
  subcell_Other="P", subcell_Secreted="P", subcell_Unknown="P",
  # C: conservation (2)
  phylop100_mean="C", gerp_mean="C",
  # V: sequence context & indel descriptors (5)
  local_entropy="V", dist_to_splice="V",
  indel_size_aa="V", indel_type="V", indel_size_relative="V"
)

FEATURE_NAME <- c(
  pli                      = "pLI",
  loeuf                    = "LOEUF",
  mis_z                    = "Missense Z",
  cds_length               = "CDS length",
  af_plddt_mean            = "pLDDT",
  af_helix_frac            = "Helix fraction",
  af_sheet_frac            = "Sheet fraction",
  af_rsa_mean              = "RSA",
  wcn_mean                 = "WCN",
  contact_density          = "Contact density",
  bridging_contacts        = "Bridging contacts",
  ss_element_break         = "SS element break",
  ss_element_frac_affected     = "SS element affected",
  helix_register_preserved = "Helix register",
  hbond_density            = "H-bond density",
  hydrophobic_exposure_risk= "Hydrophobic exposure",
  local_ss_class        = "Local SS class",
  low_plddt_region         = "Low-pLDDT region",
  in_domain                = "In domain",
  n_domains                = "Domain count",
  in_repeat                = "In repeat",
  ptm_nearby               = "PTM nearby",
  ptm_count_near           = "PTM count",
  disulfide_nearby         = "Disulfide nearby",
  in_active_site           = "Active site",
  indel_size_aa            = "Indel size (aa)",
  indel_type               = "Indel type",
  indel_size_relative      = "Indel size (relative)",
  phylop100_mean           = "phyloP",
  gerp_mean                = "GERP",
  local_entropy            = "Local entropy",
  dist_to_splice           = "Dist. to splice",
  subcell_Cytoplasm        = "Cytoplasm",
  subcell_Membrane         = "Membrane",
  subcell_Mitochondrion    = "Mitochondrion",
  subcell_Nucleus          = "Nucleus",
  subcell_Other            = "Subcell: other",
  subcell_Secreted         = "Secreted",
  subcell_Unknown          = "Subcell: unknown"
)

# Load training set
train <- fread(file.path(ROOT, "dataset/train_data/train_set.tsv"))
cd    <- readRDS(file.path(ROOT, "model/col_definitions.rds"))
feats <- cd$feature_cols
cat(sprintf("[load] train n = %d; features = %d\n", nrow(train), length(feats)))

# factors / characters -> numeric
X <- as.data.frame(train)[, feats, drop = FALSE]
for (col in names(X)) {
  if (is.factor(X[[col]]))      X[[col]] <- as.numeric(X[[col]]) - 1
  if (is.character(X[[col]]))   X[[col]] <- as.numeric(as.factor(X[[col]])) - 1
  if (is.logical(X[[col]]))     X[[col]] <- as.numeric(X[[col]])
}

M <- cor(X, use = "pairwise.complete.obs", method = "pearson")

# Order features by group, then by name
group_order <- c("S", "C", "V", "G", "P")
feat_dt <- data.table(feature = feats,
                      group   = FEATURE_GROUP[feats],
                      label   = FEATURE_NAME[feats])
feat_dt[, group := factor(group, levels = group_order)]
feat_dt <- feat_dt[order(group, feature)]
ordered_feats  <- feat_dt$feature
ordered_labels <- feat_dt$label
ordered_groups <- feat_dt$group

M <- M[ordered_feats, ordered_feats]
rownames(M) <- ordered_labels
colnames(M) <- ordered_labels

# lower-triangle long format
n <- length(ordered_labels)
df <- as.data.table(reshape2::melt(M, varnames = c("row", "col"),
                                   value.name = "r"))
df[, row := factor(row, levels = ordered_labels)]
df[, col := factor(col, levels = ordered_labels)]
df[, row_i := as.integer(row)]
df[, col_i := as.integer(col)]
df <- df[col_i <= row_i]

# axis-label colours by feature group; y reversed for diagonal
x_lab_cols <- GROUP_COLOURS[as.character(ordered_groups)]
y_lab_cols <- GROUP_COLOURS[as.character(rev(ordered_groups))]
grp_key <- data.table(grp = factor(group_order, levels = group_order),
                      col = factor(ordered_labels[1], levels = ordered_labels),
                      row = factor(ordered_labels[1], levels = ordered_labels))

hm <- ggplot(df, aes(x = col, y = row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  geom_point(data = grp_key, aes(x = col, y = row, colour = grp),
             size = 0, stroke = 0, inherit.aes = FALSE) +
  scale_colour_manual(values = GROUP_COLOURS, name = "Group",
                      breaks = group_order) +
  scale_fill_gradient2(low      = "#8B5FA3",
                       mid      = "#F4F5F7",
                       high     = "#446455",
                       midpoint = 0,
                       limits   = c(-1, 1),
                       breaks   = c(-1, -0.5, 0, 0.5, 1),
                       name     = "Pearson r") +
  scale_y_discrete(limits = rev(ordered_labels)) +
  scale_x_discrete(position = "bottom") +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  guides(fill = guide_colourbar(order = 1),
         colour = guide_legend(order = 2,
                               override.aes = list(size = 3, stroke = 0.7))) +
  theme_indelvar(base_size = 8) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1,
                                      size = 6.5, colour = x_lab_cols),
    axis.text.y        = element_text(size = 6.5, colour = y_lab_cols),
    panel.grid         = element_blank(),
    legend.position    = "right",
    legend.key.height  = unit(10, "mm"),
    legend.key.width   = unit(2.5, "mm"),
    legend.title       = element_text(size = 7.5, face = "bold",
                                      colour = COL_TEXT),
    legend.text        = element_text(size = 6.8, colour = COL_TEXT),
    plot.background    = element_rect(fill = "white", colour = NA),
    panel.background   = element_rect(fill = "white", colour = NA),
    plot.margin        = margin(6, 6, 6, 6)
  )

pS3 <- hm

FIG_NAME  <- "figS3_feature_correlation"
pdf_path  <- file.path(OUT, paste0(FIG_NAME, ".pdf"))
png_path <- file.path(OUT, paste0(FIG_NAME, ".png"))
ggsave(png_path, pS3, width = 220, height = 200, units = "mm",
       device = ragg::agg_png, dpi = 600)
# probe cairo_pdf (needs XQuartz); skip PDF if unavailable
cairo_ok <- tryCatch({ .tf <- tempfile(fileext = ".pdf"); cairo_pdf(.tf); dev.off()
                       ok <- file.exists(.tf) && file.info(.tf)$size > 0; unlink(.tf); ok },
                     error = function(e) FALSE, warning = function(w) FALSE)
if (cairo_ok) ggsave(pdf_path, pS3, width = 220, height = 200, units = "mm",
                     device = cairo_pdf, dpi = 600) else
  message("[pdf] SKIPPED ", basename(pdf_path),
          " - cairo_pdf unavailable (install XQuartz for vector PDF); 600 dpi PNG is the deliverable")
cat(sprintf("[done] wrote %s\n[done] wrote %s\n", pdf_path, png_path))
