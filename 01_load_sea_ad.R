#!/usr/bin/env Rscript
# =============================================================================
# 01_load_sea_ad.R
# -----------------------------------------------------------------------------
# Shared utilities used by the rest of the pipeline. Sourced by every analysis
# script. Provides:
#
#   - Path config (one place to edit input/output directories)
#   - load_sea_ad(): reads all per-supertype Nebula CSVs, computes per-file
#                    FDR, aggregates to broad cell type
#   - load_mouse_de(): reads the user's mouse DE results, attaches human
#                      orthologs
#   - ortho_map / map_to_human(): mouse → human ortholog mapping
#   - publication ggplot2 theme + color palette
#   - save_fig() / save_csv() helpers
#
# This script does not produce any output on its own — it defines functions
# and constants. Run 02_*.R through 06_*.R for analyses.
#
# Author: Stella Wroblewski
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

# =============================================================================
# PATH CONFIG -- edit per local install
# =============================================================================

PATHS <- list(
  # Input: SEA-AD Nebula DE CSVs (output of 00_download_sea_ad_data.sh)
  sea_ad_dir = "./data/SEA-AD",

  # Input: mouse DE table. Single CSV with one row per (gene, cell_class,
  # timepoint), columns: gene, cell_class, timepoint, logFC, pvalue, FDR.
  # See example_data/mouse_DE_example.csv for the expected schema.
  mouse_de_file = "./data/mouse_DE.csv",

  # Mouse timepoint to compare (e.g. "1WPI"); set NULL for whole table
  mouse_timepoint = "1WPI",

  # Outputs
  fig_dir = "./output/figures",
  csv_dir = "./output/tables"
)

# Create output dirs
for (d in c(PATHS$fig_dir, PATHS$csv_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# THEMES + COLORS
# =============================================================================

theme_pub <- function(base_size = 14) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 4,
                                   hjust = 0.5, margin = margin(b = 10)),
      plot.subtitle = element_text(size = base_size, hjust = 0.5,
                                   color = "gray30", margin = margin(b = 10)),
      axis.title    = element_text(face = "bold", size = base_size),
      axis.text     = element_text(size = base_size - 2, color = "black"),
      legend.title  = element_text(face = "bold", size = base_size),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin   = margin(15, 15, 15, 15)
    )
}

SPECIES_COLORS <- c("Mouse CMV" = "#C93312", "Human AD" = "#3B6FB6")

# =============================================================================
# IO HELPERS
# =============================================================================

save_fig <- function(name, plot, w = 10, h = 7) {
  out <- file.path(PATHS$fig_dir, paste0(name, ".png"))
  ggsave(out, plot, width = w, height = h, dpi = 300, bg = "white")
  message("[fig saved] ", basename(out))
}

save_csv <- function(name, df) {
  out <- file.path(PATHS$csv_dir, paste0(name, ".csv"))
  write.csv(df, out, row.names = FALSE)
  message("[csv saved] ", basename(out))
}

# =============================================================================
# ORTHOLOG MAPPING
# =============================================================================

# Curated 1:1 mouse->human ortholog map for the genes most relevant to
# cross-species vascular / glial AD comparison. For everything else, a
# simple uppercase conversion is used (correct for >90% of 1:1 orthologs).
ortho_map <- c(
  # BBB / endothelial identity
  "Cldn5"="CLDN5","Ocln"="OCLN","Tjp1"="TJP1","Cdh5"="CDH5","Ctnnb1"="CTNNB1",
  "Pecam1"="PECAM1","Flt1"="FLT1","Kdr"="KDR","Vwf"="VWF","Fn1"="FN1",
  # Transport
  "Slc2a1"="SLC2A1","Slc7a5"="SLC7A5","Abcb1"="ABCB1","Mfsd2a"="MFSD2A",
  # Pericyte
  "Pdgfrb"="PDGFRB","Rgs5"="RGS5","Cspg4"="CSPG4","Notch3"="NOTCH3",
  # ISGs
  "Stat1"="STAT1","Stat2"="STAT2","Ifit1"="IFIT1","Ifit3"="IFIT3",
  "Ifit3b"="IFIT3","Isg15"="ISG15","Ifitm3"="IFITM3","Bst2"="BST2",
  "Gbp3"="GBP3","Gbp5"="GBP5","Oasl2"="OASL","Irf2"="IRF2",
  "Ifitm1"="IFITM1","Tap2"="TAP2",
  # Microglial homeostatic
  "P2ry12"="P2RY12","Tmem119"="TMEM119","Cx3cr1"="CX3CR1","Hexb"="HEXB",
  # DAM
  "Tyrobp"="TYROBP","Apoe"="APOE","B2m"="B2M","Fth1"="FTH1",
  "Trem2"="TREM2","Cst7"="CST7","Lpl"="LPL","Cd9"="CD9",
  # Other glial
  "Ptprc"="PTPRC","Aif1"="AIF1","Gfap"="GFAP","Aqp4"="AQP4",
  "Mbp"="MBP","Mog"="MOG"
)

map_to_human <- function(mouse_genes) {
  mapped <- ortho_map[mouse_genes]
  unmapped <- is.na(mapped)
  mapped[unmapped] <- toupper(mouse_genes[unmapped])
  unname(mapped)
}

# =============================================================================
# LOAD SEA-AD NEBULA DE FILES
# =============================================================================

# Map filename prefix -> broad cell type label
.broad_for_file <- function(bn) {
  case_when(
    grepl("^Endo_Endo",            bn) ~ "Endothelial",
    grepl("^Micro-PVM_Micro-PVM",  bn) ~ "Microglia",
    grepl("^VLMC_Pericyte",        bn) ~ "Pericyte",
    grepl("^VLMC_VLMC",            bn) ~ "VLMC",
    grepl("^VLMC_SMC",             bn) ~ "SMC",
    grepl("^Astro_Astro",          bn) ~ "Astrocyte",
    grepl("^Oligo_Oligo",          bn) ~ "Oligodendrocyte",
    TRUE                               ~ "Other"
  )
}

read_nebula_file <- function(path) {
  d <- tryCatch(read.csv(path, row.names = 1, check.names = FALSE),
                error = function(e) NULL)
  if (is.null(d)) return(NULL)

  lfc_col <- grep("logFC_Continuous", colnames(d), value = TRUE)
  p_col   <- grep("p_Continuous",     colnames(d), value = TRUE)
  if (length(lfc_col) == 0 || length(p_col) == 0) return(NULL)

  out <- data.frame(
    gene   = rownames(d),
    logFC  = d[[lfc_col[1]]],
    pvalue = d[[p_col[1]]],
    stringsAsFactors = FALSE
  )
  out$FDR <- p.adjust(out$pvalue, method = "BH")
  out[!is.na(out$pvalue), , drop = FALSE]
}

#' Load SEA-AD Nebula DE results for one donor stratification.
#'
#' @param phase_dir directory containing the per-supertype CSVs
#' @param phase_label "All" | "Early" | "Late" — added as a column for tracking
#'
#' @return list with two elements:
#'   - all      : data.frame of every (gene, supertype) row
#'   - by_broad : named list of data.frames, one per broad cell type, with
#'                rows aggregated by taking the smallest-pvalue supertype
#'                per gene
load_sea_ad <- function(phase_dir, phase_label = "All") {
  if (!dir.exists(phase_dir)) {
    stop("SEA-AD directory not found: ", phase_dir)
  }

  csvs <- list.files(phase_dir, pattern = "_DE\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) {
    stop("No SEA-AD CSV files found in: ", phase_dir,
         "\nRun 00_download_sea_ad_data.sh first.")
  }
  message("Loading ", length(csvs), " SEA-AD CSVs from ", phase_dir)

  all_data <- bind_rows(lapply(csvs, function(f) {
    bn  <- basename(f)
    sup <- sub("_across_Continuous.*", "", bn)
    out <- read_nebula_file(f)
    if (is.null(out)) {
      warning("  Skipped (CPS columns not found): ", bn)
      return(NULL)
    }
    out$supertype <- sup
    out$broad     <- .broad_for_file(bn)
    out$phase     <- phase_label
    out
  }))

  # Aggregate: for each broad type and each gene, keep the supertype with
  # the smallest p-value
  by_broad <- list()
  for (bt in unique(all_data$broad)) {
    by_broad[[bt]] <- all_data %>%
      filter(.data$broad == bt) %>%
      group_by(.data$gene) %>%
      slice_min(.data$pvalue, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      as.data.frame()
  }

  message("Aggregated to ", length(by_broad), " broad cell types: ",
          paste(names(by_broad), collapse = ", "))
  list(all = all_data, by_broad = by_broad)
}

# =============================================================================
# LOAD MOUSE DE
# =============================================================================

#' Load the mouse differential expression table.
#'
#' Expected schema (CSV):
#'   gene         (mouse symbol; e.g. "Cldn5")
#'   cell_class   (e.g. "Vascular", "Immune", "Astro-Epen", "OPC-Oligo")
#'   timepoint    (e.g. "1WPI", "1MPI", "2MPI")
#'   logFC        (CMV vs Mock log2 fold change)
#'   pvalue
#'   FDR
#'
#' @param cell_class filter to one cell class (NULL = no filter)
#' @param timepoint  filter to one timepoint (NULL = no filter)
#'
#' @return data.frame with an added `human` column from map_to_human()
load_mouse_de <- function(cell_class = NULL,
                          timepoint  = PATHS$mouse_timepoint) {
  if (!file.exists(PATHS$mouse_de_file)) {
    stop("Mouse DE file not found: ", PATHS$mouse_de_file,
         "\nSee example_data/mouse_DE_example.csv for expected schema.")
  }
  de <- read.csv(PATHS$mouse_de_file, stringsAsFactors = FALSE)
  required <- c("gene","cell_class","timepoint","logFC","pvalue","FDR")
  missing  <- setdiff(required, colnames(de))
  if (length(missing) > 0) {
    stop("Mouse DE file is missing columns: ",
         paste(missing, collapse = ", "))
  }
  if (!is.null(timepoint))  de <- de[de$timepoint  == timepoint, , drop = FALSE]
  if (!is.null(cell_class)) de <- de[de$cell_class == cell_class, , drop = FALSE]
  de$human <- map_to_human(de$gene)
  de
}

# Quick sanity print when sourced
message("01_load_sea_ad.R loaded: ",
        length(ls(pattern = "^load_|^map_|^save_|^theme_|^read_")),
        " helpers available.")
