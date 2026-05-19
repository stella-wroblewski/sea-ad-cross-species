#!/usr/bin/env Rscript
# =============================================================================
# 04_cross_species_gsea.R
# -----------------------------------------------------------------------------
# Hallmark pathway GSEA on SEA-AD per-broad-cell-type DE results, with
# side-by-side comparison to mouse GSEA on matched cell types.
#
# Why this analysis matters: many cross-species concordance effects show up
# more cleanly at the pathway level than at the gene level, because pathway
# enrichment integrates information across hundreds of genes. The OXPHOS
# depletion finding (conserved across mouse and human glia) emerged from
# this analysis, not from the gene-level overlap.
#
# Two outputs:
#   1. NES heatmap (10 Hallmark pathways x 4 cell types x 2 species)
#   2. CSV with full NES + FDR per pathway / cell type / species
#
# Author: Stella Wroblewski
# =============================================================================

if (!exists("load_sea_ad")) source("scripts/01_load_sea_ad.R")

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
})

# =============================================================================
# CONFIG
# =============================================================================

# Pathways to plot in the heatmap (subset of all Hallmark). Edit to taste.
PATHWAYS_TO_PLOT <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_COAGULATION",
  "HALLMARK_ANGIOGENESIS"
)

# Matched mouse class -> human broad type
CELL_PAIRS <- list(
  list(mouse = "Vascular",   human = "Endothelial"),
  list(mouse = "Immune",     human = "Microglia"),
  list(mouse = "Astro-Epen", human = "Astrocyte"),
  list(mouse = "OPC-Oligo",  human = "Oligodendrocyte")
)

# =============================================================================
# LOAD GENE SETS
# =============================================================================

cat("Loading MSigDB Hallmark gene sets...\n")
hallmark_mm <- msigdbr(species = "Mus musculus", collection = "H") %>%
  select(.data$gs_name, gene_symbol = .data$gene_symbol) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

hallmark_hs <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(.data$gs_name, gene_symbol = .data$gene_symbol) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

# =============================================================================
# GSEA HELPER
# =============================================================================

#' Run fgsea against Hallmark for one DE table.
#' Ranking metric = sign(logFC) * -log10(pvalue), the standard choice for
#' fgsea on DE results.
run_gsea <- function(de_df, pathways, label, gene_col = "gene") {
  ranks <- with(de_df, {
    r <- sign(logFC) * -log10(pvalue + 1e-300)
    setNames(r, get(gene_col))
  })
  ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
  ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)

  res <- fgsea(pathways = pathways, stats = ranks,
               minSize = 15, maxSize = 500, nPermSimple = 10000)
  res <- as.data.frame(res)
  # fgsea returns list columns for leadingEdge; flatten for CSV write
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
  res$cell_type   <- label
  res
}

# =============================================================================
# RUN GSEA -- HUMAN SIDE
# =============================================================================

cat("\n", strrep("=", 70), "\n  HUMAN AD GSEA (SEA-AD All Donors)\n",
    strrep("=", 70), "\n", sep = "")

sea_ad <- load_sea_ad(file.path(PATHS$sea_ad_dir, "All"), "All")

human_gsea <- bind_rows(lapply(CELL_PAIRS, function(pair) {
  bt <- pair$human
  d  <- sea_ad$by_broad[[bt]]
  if (is.null(d)) return(NULL)
  cat("  ", bt, " (n=", nrow(d), " genes)\n", sep = "")
  run_gsea(d, hallmark_hs, label = bt, gene_col = "gene")
}))
human_gsea$species <- "Human AD"

# =============================================================================
# RUN GSEA -- MOUSE SIDE
# =============================================================================

cat("\n", strrep("=", 70), "\n  MOUSE CMV GSEA\n",
    strrep("=", 70), "\n", sep = "")

mouse_gsea <- bind_rows(lapply(CELL_PAIRS, function(pair) {
  d <- tryCatch(load_mouse_de(cell_class = pair$mouse),
                error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  cat("  ", pair$mouse, " (n=", nrow(d), " genes)\n", sep = "")
  run_gsea(d, hallmark_mm, label = pair$mouse, gene_col = "gene")
}))
mouse_gsea$species <- "Mouse CMV"

# Translate mouse cell_type labels to the human "broad" naming so the heatmap
# pairs up cleanly
mouse_to_human_celltype <- setNames(
  sapply(CELL_PAIRS, `[[`, "human"),
  sapply(CELL_PAIRS, `[[`, "mouse")
)
mouse_gsea$cell_type <- mouse_to_human_celltype[mouse_gsea$cell_type]

# =============================================================================
# COMBINED OUTPUT
# =============================================================================

combined <- bind_rows(mouse_gsea, human_gsea)
save_csv("04_cross_species_gsea_full", combined)

# Save a NES-only summary, easier to scan
nes_wide <- combined %>%
  filter(.data$pathway %in% PATHWAYS_TO_PLOT) %>%
  select(.data$pathway, .data$cell_type, .data$species, .data$NES, .data$padj) %>%
  pivot_wider(names_from = c(.data$species, .data$cell_type),
              values_from = c(.data$NES, .data$padj))
save_csv("04_cross_species_gsea_nes_wide", nes_wide)

# =============================================================================
# HEATMAP
# =============================================================================

plot_df <- combined %>%
  filter(.data$pathway %in% PATHWAYS_TO_PLOT) %>%
  mutate(
    pathway_short = gsub("HALLMARK_", "", .data$pathway),
    pathway_short = gsub("_", " ", .data$pathway_short),
    label = sprintf("%.2f%s", .data$NES,
                    ifelse(.data$padj < 0.05, "*",
                           ifelse(.data$padj < 0.25, "+", "")))
  )

plot_df$pathway_short <- factor(plot_df$pathway_short,
  levels = rev(gsub("_", " ", gsub("HALLMARK_", "", PATHWAYS_TO_PLOT))))
plot_df$species   <- factor(plot_df$species, levels = c("Mouse CMV", "Human AD"))
plot_df$cell_type <- factor(plot_df$cell_type,
  levels = sapply(CELL_PAIRS, `[[`, "human"))

fig <- ggplot(plot_df, aes(x = .data$species, y = .data$pathway_short,
                            fill = .data$NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = .data$label), size = 3.4, fontface = "bold") +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white",
                       high = "#C93312", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish,
                       name = "NES") +
  facet_wrap(~ cell_type, nrow = 1) +
  theme_pub(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cross-species Hallmark GSEA",
       subtitle = "Mouse CMV (1WPI) vs Human AD (SEA-AD All Donors, CPS coefficient)",
       x = "", y = "")

save_fig("04_cross_species_gsea_heatmap", fig, w = 13, h = 8)

# =============================================================================
# CONSOLE SUMMARY: SPECIFIC FINDINGS
# =============================================================================

cat("\n", strrep("=", 70), "\n  KEY FINDINGS\n",
    strrep("=", 70), "\n", sep = "")

for (path in c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
               "HALLMARK_INTERFERON_ALPHA_RESPONSE",
               "HALLMARK_TNFA_SIGNALING_VIA_NFKB")) {
  cat("\n  ", path, "\n", sep = "")
  combined %>%
    filter(.data$pathway == path) %>%
    arrange(.data$cell_type, .data$species) %>%
    mutate(NES = round(.data$NES, 2), padj = signif(.data$padj, 2)) %>%
    select(.data$cell_type, .data$species, .data$NES, .data$padj) %>%
    print(row.names = FALSE)
}

cat("\nDone.\n")
