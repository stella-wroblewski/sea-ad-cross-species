#!/usr/bin/env Rscript
# =============================================================================
# 02_hypergeometric_and_concordance.R
# -----------------------------------------------------------------------------
# Cross-species overlap analyses on the All-Donors stratification.
#
# Analysis 1 -- Hypergeometric overlap test:
#   Q: Are mouse DE genes enriched among human AD DE genes more than expected
#      by chance?
#   Tests at FDR < 0.05, 0.10, 0.25 thresholds.
#
# Analysis 2 -- Directional concordance:
#   For genes tested in both species: do logFC signs agree? Spearman
#   correlation on logFC values.
#
# Both analyses are matched-cell-type: mouse vascular vs SEA-AD endothelial,
# mouse immune vs SEA-AD microglia, etc. The matching is hard-coded below.
#
# Author: Stella Wroblewski
# =============================================================================

if (!exists("load_sea_ad")) source("scripts/01_load_sea_ad.R")

# Hard-coded matched pairs: mouse cell_class -> SEA-AD broad cell type
CELL_PAIRS <- list(
  list(mouse = "Vascular",   human = "Endothelial"),
  list(mouse = "Immune",     human = "Microglia"),
  list(mouse = "Astro-Epen", human = "Astrocyte"),
  list(mouse = "OPC-Oligo",  human = "Oligodendrocyte")
)

# =============================================================================
# LOAD
# =============================================================================

sea_ad   <- load_sea_ad(file.path(PATHS$sea_ad_dir, "All"), phase_label = "All")

# =============================================================================
# ANALYSIS 1: HYPERGEOMETRIC OVERLAP
# =============================================================================

hyper_test <- function(mouse_de, human_de, thresholds = c(0.05, 0.10, 0.25)) {
  # Universe = genes tested in both
  universe <- intersect(unique(mouse_de$human), unique(human_de$gene))

  out <- bind_rows(lapply(thresholds, function(thr) {
    ms <- mouse_de %>% filter(.data$FDR < thr) %>%
      pull(.data$human) %>% na.omit() %>% unique() %>% intersect(universe)
    hs <- human_de %>% filter(.data$FDR < thr) %>%
      pull(.data$gene) %>% unique() %>% intersect(universe)
    ov <- length(intersect(ms, hs))
    expected <- length(ms) * length(hs) / max(length(universe), 1)
    p_val <- phyper(ov - 1, length(hs),
                    length(universe) - length(hs),
                    length(ms), lower.tail = FALSE)
    data.frame(
      threshold = sprintf("FDR<%.2f", thr),
      mouse_n   = length(ms),
      human_n   = length(hs),
      overlap   = ov,
      expected  = round(expected, 2),
      fold      = round(ov / max(expected, 0.01), 2),
      p_value   = p_val,
      stringsAsFactors = FALSE
    )
  }))
  list(table = out, universe_size = length(universe))
}

cat("\n", strrep("=", 70), "\n",
    "ANALYSIS 1: HYPERGEOMETRIC OVERLAP\n",
    strrep("=", 70), "\n", sep = "")

hyper_all <- list()
for (pair in CELL_PAIRS) {
  mname <- pair$mouse; hname <- pair$human
  cat("\n  ", mname, " (mouse) <-> ", hname, " (human)\n", sep = "")

  mouse_de <- tryCatch(load_mouse_de(cell_class = mname),
                       error = function(e) { message(conditionMessage(e)); NULL })
  human_de <- sea_ad$by_broad[[hname]]
  if (is.null(mouse_de) || is.null(human_de) || nrow(mouse_de) == 0) {
    cat("  skipping (missing data)\n"); next
  }

  res <- hyper_test(mouse_de, human_de)
  cat("  Universe (genes tested in both): ", res$universe_size, "\n", sep = "")
  print(res$table, row.names = FALSE)
  res$table$pair <- paste0(mname, "_vs_", hname)
  hyper_all[[mname]] <- res$table
}

save_csv("01_hypergeometric", bind_rows(hyper_all))

# =============================================================================
# ANALYSIS 2: DIRECTIONAL CONCORDANCE
# =============================================================================

concordance <- function(mouse_de, human_de, label,
                        sig_filter = "any") {
  # Join on human ortholog
  merged <- mouse_de %>%
    select(.data$human, mouse_logFC = .data$logFC,
                        mouse_FDR   = .data$FDR) %>%
    distinct(.data$human, .keep_all = TRUE) %>%
    inner_join(
      human_de %>%
        select(human = .data$gene,
               human_logFC = .data$logFC,
               human_FDR   = .data$FDR) %>%
        distinct(.data$human, .keep_all = TRUE),
      by = "human"
    ) %>%
    filter(!is.na(.data$mouse_logFC), !is.na(.data$human_logFC))

  if (sig_filter == "mouse_sig") merged <- merged %>% filter(.data$mouse_FDR < 0.05)
  if (sig_filter == "human_sig") merged <- merged %>% filter(.data$human_FDR < 0.05)
  if (sig_filter == "both_nominal")
    merged <- merged %>% filter(.data$mouse_FDR < 0.10, .data$human_FDR < 0.10)

  if (nrow(merged) < 5) {
    cat("  [", label, "] too few genes (n=", nrow(merged), ")\n", sep = "")
    return(NULL)
  }

  rho <- cor.test(merged$mouse_logFC, merged$human_logFC, method = "spearman",
                  exact = FALSE)
  concordant <- sign(merged$mouse_logFC) == sign(merged$human_logFC)
  pct_concordant <- mean(concordant) * 100

  # Binomial test against 50%
  bt <- binom.test(sum(concordant), nrow(merged), p = 0.5)

  cat(sprintf("  [%s] n=%d  rho=%+.3f (p=%.2e)  concordant=%.1f%% (binom p=%.2e)\n",
              label, nrow(merged), rho$estimate, rho$p.value,
              pct_concordant, bt$p.value))

  # Scatter
  fig <- ggplot(merged, aes(x = .data$mouse_logFC, y = .data$human_logFC)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_point(alpha = 0.4, size = 1) +
    geom_smooth(method = "lm", se = FALSE, color = SPECIES_COLORS[["Mouse CMV"]],
                linewidth = 0.7) +
    labs(title = paste0("Cross-species concordance: ", label),
         subtitle = sprintf("rho = %+.3f  |  %.1f%% concordant  |  n = %d",
                            rho$estimate, pct_concordant, nrow(merged)),
         x = "Mouse CMV logFC", y = "Human AD CPS logFC") +
    theme_pub(13)
  save_fig(paste0("02_concordance_", gsub("[^A-Za-z0-9]+", "_", label)), fig,
           w = 7, h = 6)

  data.frame(
    pair             = label,
    n_genes          = nrow(merged),
    rho              = unname(rho$estimate),
    rho_p            = rho$p.value,
    pct_concordant   = pct_concordant,
    binom_p          = bt$p.value,
    stringsAsFactors = FALSE
  )
}

cat("\n", strrep("=", 70), "\n",
    "ANALYSIS 2: DIRECTIONAL CONCORDANCE\n",
    strrep("=", 70), "\n", sep = "")

conc_all <- list()
for (pair in CELL_PAIRS) {
  mname <- pair$mouse; hname <- pair$human
  mouse_de <- tryCatch(load_mouse_de(cell_class = mname),
                       error = function(e) NULL)
  human_de <- sea_ad$by_broad[[hname]]
  if (is.null(mouse_de) || is.null(human_de)) next

  conc_all[[paste0(mname, "_vs_", hname)]] <-
    concordance(mouse_de, human_de, paste0(mname, " vs ", hname))
}

save_csv("02_concordance", bind_rows(conc_all))

cat("\nDone. Tables in ", PATHS$csv_dir, ", figures in ", PATHS$fig_dir, "\n")
