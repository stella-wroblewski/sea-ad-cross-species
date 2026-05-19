#!/usr/bin/env Rscript
# =============================================================================
# 05_dam_microglia_comparison.R
# -----------------------------------------------------------------------------
# Disease-Associated Microglia (DAM) gene-by-gene comparison.
#
# Background: in mouse models of AD, microglia transition from a homeostatic
# state through DAM Stage 1 (Tyrobp+ Apoe+) into DAM Stage 2 (Trem2+ Cst7+ Lpl+
# Cd9+). Keren-Shaul et al. 2017 (Cell) characterized this. In human AD,
# similar but not identical microglial states are observed and reported in
# SEA-AD as elevated coefficients on DAM-associated genes across CPS.
#
# This script compares the mouse logFC to the SEA-AD microglial logFC for each
# of these canonical genes, in three modules:
#   - Homeostatic (should DOWN in DAM)
#   - DAM Stage 1
#   - DAM Stage 2
#
# Output: faceted bar plot + summary table of which genes are concordant.
#
# Author: Stella Wroblewski
# =============================================================================

if (!exists("load_sea_ad")) source("scripts/01_load_sea_ad.R")

suppressPackageStartupMessages({
  library(tibble)
})

# =============================================================================
# DAM GENES
# =============================================================================

DAM_GENES <- tribble(
  ~mouse,    ~human,    ~module,
  # Homeostatic — expected DOWN in DAM
  "P2ry12",  "P2RY12",  "Homeostatic",
  "Tmem119", "TMEM119", "Homeostatic",
  "Cx3cr1",  "CX3CR1",  "Homeostatic",
  "Hexb",    "HEXB",    "Homeostatic",
  # DAM Stage 1
  "Tyrobp",  "TYROBP",  "DAM Stage 1",
  "Apoe",    "APOE",    "DAM Stage 1",
  "B2m",     "B2M",     "DAM Stage 1",
  "Fth1",    "FTH1",    "DAM Stage 1",
  # DAM Stage 2
  "Trem2",   "TREM2",   "DAM Stage 2",
  "Cst7",    "CST7",    "DAM Stage 2",
  "Lpl",     "LPL",     "DAM Stage 2",
  "Cd9",     "CD9",     "DAM Stage 2"
)

# =============================================================================
# LOAD
# =============================================================================

sea_ad <- load_sea_ad(file.path(PATHS$sea_ad_dir, "All"), "All")
micro_h <- sea_ad$by_broad[["Microglia"]]
mouse_imm <- load_mouse_de(cell_class = "Immune")

# =============================================================================
# JOIN
# =============================================================================

dam <- DAM_GENES %>%
  left_join(
    mouse_imm %>%
      select(mouse = .data$gene, mouse_logFC = .data$logFC,
             mouse_FDR = .data$FDR) %>%
      distinct(.data$mouse, .keep_all = TRUE),
    by = "mouse"
  ) %>%
  left_join(
    micro_h %>%
      select(human = .data$gene, human_logFC = .data$logFC,
             human_FDR = .data$FDR) %>%
      distinct(.data$human, .keep_all = TRUE),
    by = "human"
  ) %>%
  mutate(
    has_data    = !is.na(.data$mouse_logFC) & !is.na(.data$human_logFC),
    concordant  = sign(.data$mouse_logFC) == sign(.data$human_logFC),
    in_dam_dir  = case_when(
      .data$module == "Homeostatic"     ~ .data$mouse_logFC < 0 & .data$human_logFC < 0,
      .data$module %in% c("DAM Stage 1", "DAM Stage 2") ~
                                          .data$mouse_logFC > 0 & .data$human_logFC > 0,
      TRUE                              ~ NA
    )
  ) %>%
  as.data.frame()

save_csv("05_DAM_genes", dam)

# Summary
cat("\n  DAM concordance summary:\n")
with_data <- dam[dam$has_data, ]
cat(sprintf("    Genes with data in both species: %d / %d\n",
            nrow(with_data), nrow(dam)))
cat(sprintf("    Sign-concordant: %d / %d (%.0f%%)\n",
            sum(with_data$concordant), nrow(with_data),
            mean(with_data$concordant) * 100))
cat(sprintf("    In expected DAM direction in both: %d / %d (%.0f%%)\n",
            sum(with_data$in_dam_dir), nrow(with_data),
            mean(with_data$in_dam_dir) * 100))
cat("\n  By module:\n")
dam %>%
  group_by(.data$module) %>%
  summarise(
    n         = sum(.data$has_data),
    n_conc    = sum(.data$concordant, na.rm = TRUE),
    pct_conc  = round(mean(.data$concordant, na.rm = TRUE) * 100, 0),
    n_in_dam  = sum(.data$in_dam_dir, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  print(row.names = FALSE)

# =============================================================================
# PLOT
# =============================================================================

plot_df <- dam %>%
  filter(.data$has_data) %>%
  select(.data$mouse, .data$module, .data$mouse_logFC, .data$human_logFC,
         .data$mouse_FDR, .data$human_FDR) %>%
  pivot_longer(
    cols = c(.data$mouse_logFC, .data$human_logFC),
    names_to = "species", values_to = "logFC"
  ) %>%
  mutate(
    species = ifelse(.data$species == "mouse_logFC", "Mouse CMV", "Human AD"),
    FDR     = ifelse(.data$species == "Mouse CMV",  .data$mouse_FDR,
                                                     .data$human_FDR),
    sig_marker = ifelse(.data$FDR < 0.05, "*", "")
  )

plot_df$mouse   <- factor(plot_df$mouse,   levels = rev(DAM_GENES$mouse))
plot_df$species <- factor(plot_df$species, levels = c("Mouse CMV", "Human AD"))
plot_df$module  <- factor(plot_df$module,  levels = c("Homeostatic",
                                                       "DAM Stage 1",
                                                       "DAM Stage 2"))

fig <- ggplot(plot_df, aes(x = .data$mouse, y = .data$logFC,
                            fill = .data$species)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  geom_text(aes(label = .data$sig_marker, y = .data$logFC * 1.05),
            position = position_dodge(0.7), size = 4, fontface = "bold") +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = SPECIES_COLORS) +
  facet_grid(rows = vars(.data$module), scales = "free_y", space = "free_y") +
  coord_flip() +
  theme_pub(13) +
  theme(strip.text.y = element_text(angle = 0, face = "bold")) +
  labs(title = "Microglial DAM program: mouse CMV vs human AD",
       subtitle = "* indicates FDR < 0.05 within each species",
       x = "", y = "logFC")

save_fig("05_DAM_comparison", fig, w = 9, h = 8)

cat("\nDone.\n")
