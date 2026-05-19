#!/usr/bin/env Rscript
# =============================================================================
# run_all.R
# -----------------------------------------------------------------------------
# Wrapper that runs all SEA-AD cross-species analyses end to end.
#
# Prerequisites:
#   1. SEA-AD Nebula CSVs downloaded:
#        bash scripts/00_download_sea_ad_data.sh
#   2. Mouse DE table available at PATHS$mouse_de_file (see 01_load_sea_ad.R)
#      with the expected schema. See example_data/mouse_DE_example.csv.
#
# Usage:
#   Rscript run_all.R
#
# Author: Stella Wroblewski
# =============================================================================

# Always run from the repository root so relative paths resolve correctly.
# Confirm we're in the right place by checking for the scripts directory.
if (!dir.exists("scripts")) {
  stop("run_all.R must be run from the repository root. ",
       "Current working dir: ", getwd())
}

t0 <- Sys.time()
cat("\n", strrep("#", 70), "\n",
    "  SEA-AD cross-species analysis pipeline\n",
    "  started ", format(t0), "\n",
    strrep("#", 70), "\n\n", sep = "")

steps <- c(
  "scripts/02_hypergeometric_and_concordance.R",
  "scripts/03_targeted_gene_heatmap.R",
  "scripts/04_cross_species_gsea.R",
  "scripts/05_dam_microglia_comparison.R",
  "scripts/06_phase_stratified_analysis.R"
)

for (step in steps) {
  cat("\n", strrep(">", 70), "\n",
      "  Running ", step, "\n",
      strrep(">", 70), "\n", sep = "")
  source(step, echo = FALSE)
}

t1 <- Sys.time()
cat("\n", strrep("#", 70), "\n",
    "  Pipeline complete\n",
    "  total time: ", round(as.numeric(difftime(t1, t0, units = "mins")), 1),
    " min\n",
    strrep("#", 70), "\n", sep = "")
