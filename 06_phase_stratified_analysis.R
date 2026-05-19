#!/usr/bin/env Rscript
# =============================================================================
# 06_phase_stratified_analysis.R
# -----------------------------------------------------------------------------
# Compares the mouse CMV signature to SEA-AD Early-phase donors, Late-phase
# donors, and the combined All-donors result side by side.
#
# Why this is the most important script in the pipeline:
# ------------------------------------------------------
# SEA-AD's all-donors Nebula coefficients average across donors spanning the
# entire ADNC spectrum (no-AD through high). The Gabitto et al. paper shows
# that the disease trajectory has two distinct phases:
#
#   Early phase: inflammatory microglia, reactive astrocytes, IFN signatures,
#                vascular stress -- slow pathology accumulation
#   Late phase:  exponential pathology, neuronal loss, IFN signatures REVERSE
#                (anti-inflammatory shift), vascular response inverts
#
# When you average over both phases, opposing signals cancel out. We saw this
# concretely: endothelial cells had 0 genes at FDR<0.05 in all-donors but 20
# genes at FDR<0.05 in early-donors -- not because the early-phase signal is
# weak, but because the late phase actively reverses it.
#
# An acute mouse model in young animals should match early-phase initiating
# events. Using all-donors data dilutes the relevant signal and risks a
# false-negative conclusion. This script demonstrates that.
#
# Outputs:
#   - GSEA NES heatmap across All / Early / Late for matched cell types
#   - Per-gene concordance table
#   - Side-by-side scatter plots
#
# Author: Stella Wroblewski
# =============================================================================

if (!exists("load_sea_ad")) source("scripts/01_load_sea_ad.R")

suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(patchwork)
})

# =============================================================================
# CONFIG
# =============================================================================

PATHWAYS_TO_PLOT <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT"
)

# Matched cell-type pairs
CELL_PAIRS <- list(
  list(mouse = "Vascular",   human = "Endothelial"),
  list(mouse = "Immune",     human = "Microglia"),
  list(mouse = "Astro-Epen", human = "Astrocyte"),
  list(mouse = "OPC-Oligo",  human = "Oligodendrocyte")
)

# =============================================================================
# LOAD ALL THREE STRATIFICATIONS
# =============================================================================

cat("\nLoading SEA-AD All / Early / Late...\n")
sea_all   <- load_sea_ad(file.path(PATHS$sea_ad_dir, "All"),   "All")
sea_early <- load_sea_ad(file.path(PATHS$sea_ad_dir, "Early"), "Early")
sea_late  <- load_sea_ad(file.path(PATHS$sea_ad_dir, "Late"),  "Late")

# =============================================================================
# FDR<0.05 COUNTS PER PHASE
# =============================================================================

cat("\n  Genes at FDR<0.05 by phase x cell type:\n")
cat(sprintf("  %-18s %8s %8s %8s\n", "Cell type", "All", "Early", "Late"))
cat("  ", strrep("-", 50), "\n")
for (bt in c("Endothelial", "Microglia", "Astrocyte",
             "Oligodendrocyte", "Pericyte")) {
  n_a <- if (!is.null(sea_all$by_broad[[bt]]))   sum(sea_all$by_broad[[bt]]$FDR   < 0.05) else NA
  n_e <- if (!is.null(sea_early$by_broad[[bt]])) sum(sea_early$by_broad[[bt]]$FDR < 0.05) else NA
  n_l <- if (!is.null(sea_late$by_broad[[bt]]))  sum(sea_late$by_broad[[bt]]$FDR  < 0.05) else NA
  cat(sprintf("  %-18s %8s %8s %8s\n", bt,
              ifelse(is.na(n_a), "-", n_a),
              ifelse(is.na(n_e), "-", n_e),
              ifelse(is.na(n_l), "-", n_l)))
}

# =============================================================================
# GSEA PER PHASE
# =============================================================================

hallmark_hs <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  select(.data$gs_name, .data$gene_symbol) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

hallmark_mm <- msigdbr(species = "Mus musculus", collection = "H") %>%
  select(.data$gs_name, .data$gene_symbol) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

run_gsea <- function(de_df, pathways, gene_col = "gene") {
  ranks <- with(de_df, {
    r <- sign(logFC) * -log10(pvalue + 1e-300)
    setNames(r, get(gene_col))
  })
  ranks <- ranks[!is.na(ranks) & is.finite(ranks)]
  ranks <- sort(ranks[!duplicated(names(ranks))], decreasing = TRUE)
  as.data.frame(fgsea(pathways = pathways, stats = ranks,
                      minSize = 15, maxSize = 500, nPermSimple = 10000))
}

cat("\nRunning GSEA across phases x cell types...\n")
gsea_results <- list()
for (pair in CELL_PAIRS) {
  bt <- pair$human
  for (phase_name in c("All", "Early", "Late")) {
    phase_data <- switch(phase_name,
                         "All"   = sea_all,
                         "Early" = sea_early,
                         "Late"  = sea_late)
    de <- phase_data$by_broad[[bt]]
    if (is.null(de) || nrow(de) == 0) next
    g <- run_gsea(de, hallmark_hs, gene_col = "gene")
    g$cell_type <- bt
    g$phase     <- phase_name
    g$species   <- "Human AD"
    gsea_results[[paste(phase_name, bt)]] <- g
  }

  # Mouse comparator on the matching cell class
  mouse_de <- tryCatch(load_mouse_de(cell_class = pair$mouse),
                       error = function(e) NULL)
  if (!is.null(mouse_de) && nrow(mouse_de) > 0) {
    g <- run_gsea(mouse_de, hallmark_mm, gene_col = "gene")
    g$cell_type <- bt
    g$phase     <- "Mouse CMV"
    g$species   <- "Mouse CMV"
    gsea_results[[paste("Mouse", bt)]] <- g
  }
}

combined <- bind_rows(gsea_results)
# leadingEdge can be a list column -- drop or collapse for csv
if ("leadingEdge" %in% names(combined)) {
  combined$leadingEdge <- sapply(combined$leadingEdge, paste, collapse = ";")
}
save_csv("06_phase_stratified_gsea_full", combined)

# =============================================================================
# HEATMAP
# =============================================================================

plot_df <- combined %>%
  filter(.data$pathway %in% PATHWAYS_TO_PLOT) %>%
  mutate(
    pathway_short = gsub("HALLMARK_", "", .data$pathway),
    pathway_short = gsub("_", " ", .data$pathway_short),
    label = sprintf("%.2f%s", .data$NES,
                    ifelse(is.na(.data$padj), "",
                           ifelse(.data$padj < 0.05, "*",
                                  ifelse(.data$padj < 0.25, "+", ""))))
  )

plot_df$pathway_short <- factor(plot_df$pathway_short,
  levels = rev(gsub("_", " ", gsub("HALLMARK_", "", PATHWAYS_TO_PLOT))))
plot_df$phase     <- factor(plot_df$phase,
                            levels = c("Mouse CMV", "All", "Early", "Late"))
plot_df$cell_type <- factor(plot_df$cell_type,
                            levels = sapply(CELL_PAIRS, `[[`, "human"))

fig <- ggplot(plot_df, aes(x = .data$phase, y = .data$pathway_short,
                            fill = .data$NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = .data$label), size = 3.2, fontface = "bold") +
  scale_fill_gradient2(low = "#3B6FB6", mid = "white",
                       high = "#C93312", midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish,
                       name = "NES") +
  facet_wrap(~ cell_type, nrow = 1) +
  theme_pub(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Phase-stratified cross-species GSEA",
       subtitle = "Mouse CMV vs Human AD (SEA-AD All / Early / Late donors)",
       x = "", y = "")

save_fig("06_phase_stratified_heatmap", fig, w = 14, h = 7)

# =============================================================================
# PER-GENE PHASE CONCORDANCE -- 17 BBB-RELEVANT GENES
# =============================================================================

cat("\nBuilding per-gene phase concordance table...\n")

BBB_GENES <- tribble(
  ~mouse,    ~human,
  "Cldn5",   "CLDN5",
  "Ocln",    "OCLN",
  "Tjp1",    "TJP1",
  "Cdh5",    "CDH5",
  "Pecam1",  "PECAM1",
  "Flt1",    "FLT1",
  "Kdr",     "KDR",
  "Vwf",     "VWF",
  "Slc2a1",  "SLC2A1",
  "Slc7a5",  "SLC7A5",
  "Abcb1",   "ABCB1",
  "Pdgfrb",  "PDGFRB",
  "Rgs5",    "RGS5",
  "Stat1",   "STAT1",
  "Ifit1",   "IFIT1",
  "Isg15",   "ISG15",
  "Ifitm3",  "IFITM3"
)

mouse_vasc <- load_mouse_de(cell_class = "Vascular")

lookup_human <- function(phase_data, human_gene) {
  d <- phase_data$by_broad[["Endothelial"]]
  if (is.null(d)) return(c(NA_real_, NA_real_))
  r <- d[d$gene == human_gene, ]
  if (nrow(r) == 0) return(c(NA_real_, NA_real_))
  c(r$logFC[1], r$FDR[1])
}

bbb_tbl <- BBB_GENES %>%
  rowwise() %>%
  mutate(
    mouse_logFC = {
      m <- mouse_vasc[mouse_vasc$gene == .data$mouse, ]
      if (nrow(m) == 0) NA_real_ else m$logFC[1]
    },
    All_logFC   = lookup_human(sea_all,   .data$human)[1],
    Early_logFC = lookup_human(sea_early, .data$human)[1],
    Late_logFC  = lookup_human(sea_late,  .data$human)[1]
  ) %>%
  ungroup() %>%
  as.data.frame()

# Concordance per phase
for (col in c("All_logFC", "Early_logFC", "Late_logFC")) {
  has <- !is.na(bbb_tbl$mouse_logFC) & !is.na(bbb_tbl[[col]])
  pct <- mean(sign(bbb_tbl$mouse_logFC[has]) == sign(bbb_tbl[[col]][has])) * 100
  cat(sprintf("  Endothelial BBB concordance (mouse vs %s): %.0f%% (n=%d)\n",
              sub("_logFC", "", col), pct, sum(has)))
}

save_csv("06_phase_BBB_concordance", bbb_tbl)

cat("\nDone.\n")
cat("\nKey output to look at: ", PATHS$fig_dir,
    "/06_phase_stratified_heatmap.png\n", sep = "")
cat("If the early-phase column matches mouse CMV and the late-phase column\n")
cat("reverses, the phase-cancellation interpretation is supported.\n")
