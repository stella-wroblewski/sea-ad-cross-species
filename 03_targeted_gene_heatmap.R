#!/usr/bin/env Rscript
# =============================================================================
# 03_targeted_gene_heatmap.R
# -----------------------------------------------------------------------------
# Side-by-side logFC heatmap for ~30 curated genes of interest (BBB junctions,
# endothelial identity, transporters, pericyte markers, ISGs). Mouse vs human
# in matched columns.
#
# This is the most directly interpretable analysis. Unlike the genome-wide
# concordance, this asks: do specific, biologically meaningful genes show the
# expected directional pattern?
#
# Author: Stella Wroblewski
# =============================================================================

if (!exists("load_sea_ad")) source("scripts/01_load_sea_ad.R")

suppressPackageStartupMessages({
  library(tibble)
})

# =============================================================================
# CURATED GENE LIST
# =============================================================================

GENES <- tribble(
  ~mouse,    ~human,    ~category,           ~human_cell,
  # BBB junctions
  "Cldn5",   "CLDN5",   "BBB junction",      "Endothelial",
  "Ocln",    "OCLN",    "BBB junction",      "Endothelial",
  "Tjp1",    "TJP1",    "BBB junction",      "Endothelial",
  "Cdh5",    "CDH5",    "BBB junction",      "Endothelial",
  # Endothelial identity
  "Pecam1",  "PECAM1",  "Endothelial ID",    "Endothelial",
  "Flt1",    "FLT1",    "Endothelial ID",    "Endothelial",
  "Kdr",     "KDR",     "Endothelial ID",    "Endothelial",
  "Vwf",     "VWF",     "Endothelial ID",    "Endothelial",
  # Transporters
  "Slc2a1",  "SLC2A1",  "Transport",         "Endothelial",
  "Slc7a5",  "SLC7A5",  "Transport",         "Endothelial",
  "Abcb1",   "ABCB1",   "Transport",         "Endothelial",
  # Pericyte
  "Pdgfrb",  "PDGFRB",  "Pericyte",          "Pericyte",
  "Rgs5",    "RGS5",    "Pericyte",          "Pericyte",
  "Notch3",  "NOTCH3",  "Pericyte",          "Pericyte",
  # ISGs
  "Stat1",   "STAT1",   "ISG",               "Endothelial",
  "Ifit1",   "IFIT1",   "ISG",               "Endothelial",
  "Ifit3",   "IFIT3",   "ISG",               "Endothelial",
  "Isg15",   "ISG15",   "ISG",               "Endothelial",
  "Ifitm3",  "IFITM3",  "ISG",               "Endothelial",
  # Microglial homeostatic
  "P2ry12",  "P2RY12",  "Homeostatic uG",    "Microglia",
  "Tmem119", "TMEM119", "Homeostatic uG",    "Microglia",
  "Cx3cr1",  "CX3CR1",  "Homeostatic uG",    "Microglia",
  # DAM
  "Tyrobp",  "TYROBP",  "DAM",               "Microglia",
  "Apoe",    "APOE",    "DAM",               "Microglia",
  "Trem2",   "TREM2",   "DAM",               "Microglia",
  "Cst7",    "CST7",    "DAM",               "Microglia",
  "Lpl",     "LPL",     "DAM",               "Microglia"
)

# =============================================================================
# LOAD
# =============================================================================

sea_ad <- load_sea_ad(file.path(PATHS$sea_ad_dir, "All"), "All")

# All mouse cell classes that might be needed
mouse_de_all <- bind_rows(
  tryCatch(load_mouse_de(cell_class = "Vascular"), error = function(e) NULL),
  tryCatch(load_mouse_de(cell_class = "Immune"),   error = function(e) NULL)
)

# =============================================================================
# BUILD MATRIX
# =============================================================================

# For each gene, look up mouse logFC and human logFC. Determine which mouse
# cell class to read based on whether it's a microglial or vascular gene.
GENES$mouse_class <- ifelse(GENES$human_cell == "Microglia", "Immune", "Vascular")

heatmap_df <- GENES %>%
  rowwise() %>%
  mutate(
    mouse_logFC = {
      m <- mouse_de_all %>%
        filter(.data$cell_class == .env$mouse_class) %>%
        filter(.data$gene == .env$mouse)
      if (nrow(m) == 0) NA_real_ else m$logFC[1]
    },
    mouse_FDR = {
      m <- mouse_de_all %>%
        filter(.data$cell_class == .env$mouse_class) %>%
        filter(.data$gene == .env$mouse)
      if (nrow(m) == 0) NA_real_ else m$FDR[1]
    },
    human_logFC = {
      h <- sea_ad$by_broad[[.env$human_cell]]
      if (is.null(h)) return(NA_real_)
      r <- h %>% filter(.data$gene == .env$human)
      if (nrow(r) == 0) NA_real_ else r$logFC[1]
    },
    human_FDR = {
      h <- sea_ad$by_broad[[.env$human_cell]]
      if (is.null(h)) return(NA_real_)
      r <- h %>% filter(.data$gene == .env$human)
      if (nrow(r) == 0) NA_real_ else r$FDR[1]
    }
  ) %>%
  ungroup() %>%
  as.data.frame()

# Summary stats
n_with_data <- sum(!is.na(heatmap_df$mouse_logFC) & !is.na(heatmap_df$human_logFC))
concordant  <- sign(heatmap_df$mouse_logFC) == sign(heatmap_df$human_logFC)
pct_conc    <- mean(concordant, na.rm = TRUE) * 100
cat(sprintf("\n%d/%d genes have data in both species\n",
            n_with_data, nrow(heatmap_df)))
cat(sprintf("%.1f%% of those are directionally concordant\n", pct_conc))

save_csv("03_targeted_genes", heatmap_df)

# =============================================================================
# PLOT
# =============================================================================

# Long format for ggplot
plot_df <- heatmap_df %>%
  select(.data$mouse, .data$category, .data$mouse_logFC,
         .data$human_logFC, .data$mouse_FDR, .data$human_FDR) %>%
  pivot_longer(cols = c(.data$mouse_logFC, .data$human_logFC),
               names_to = "species", values_to = "logFC") %>%
  mutate(
    species = ifelse(.data$species == "mouse_logFC", "Mouse CMV", "Human AD"),
    FDR     = ifelse(.data$species == "Mouse CMV", .data$mouse_FDR,
                                                   .data$human_FDR),
    sig     = ifelse(is.na(.data$FDR), "",
                     ifelse(.data$FDR < 0.05, "*",
                            ifelse(.data$FDR < 0.25, "+", "")))
  )

# Order genes by category
plot_df$mouse    <- factor(plot_df$mouse,    levels = rev(GENES$mouse))
plot_df$species  <- factor(plot_df$species,  levels = c("Mouse CMV", "Human AD"))
plot_df$category <- factor(plot_df$category, levels = unique(GENES$category))

fig <- ggplot(plot_df, aes(x = .data$species, y = .data$mouse,
                            fill = .data$logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = .data$sig), size = 5, fontface = "bold") +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white",
                       high = "#C93312", midpoint = 0,
                       limits = c(-2, 2), oob = scales::squish,
                       name = "logFC", na.value = "gray90") +
  facet_grid(rows = vars(.data$category), scales = "free_y", space = "free_y",
             switch = "y") +
  theme_pub(12) +
  theme(
    strip.text.y.left  = element_text(angle = 0, hjust = 1, face = "bold"),
    strip.placement    = "outside",
    axis.text.y        = element_text(face = "italic")
  ) +
  labs(title = "Cross-species logFC heatmap for targeted genes",
       subtitle = sprintf("%.0f%% directional concordance | * FDR<0.05 + FDR<0.25",
                          pct_conc),
       x = "", y = "")

save_fig("03_targeted_heatmap", fig, w = 7, h = 11)

cat("\nDone.\n")
