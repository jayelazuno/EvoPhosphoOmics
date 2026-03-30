# Author: Joshua Ayelazuno
# 2025-10-21
# Title: Time-Course LFC Summarization (per species) for -Pi Transcriptomes
# ==============================================================================
# Description:
#   This script compiles per-timepoint DESeq2 outputs into matrices that are easy
#   to visualize. For a given species, it scans:
#     06.results/<Species>/paired/
#   for per-contrast files named like:
#     <YYYYMMDD>__<Species>__paired__tXXXX_vs_t0000__DEGs_all__fdr0010_lfc1000.csv
#
#   From each file it extracts log2FoldChange values and constructs:
#     • LFC wide matrix  (rows = gene_id; columns = tXXXX; values = log2FC)
#     • LFC tidy-long table (gene_id, tp, log2FC) for ggplot/faceting
#   Missing values can be zero-filled (zero_fill = TRUE) so every gene has a
#   complete time grid.
#
# Output structure (per species):
#   06.results/<Species>/paired/combined/
#     <DATE>__<Species>__paired__LFC_matrix.csv
#     <DATE>__<Species>__paired__LFC_tidy_long.csv
#
# Inputs & assumptions:
#   • Per-timepoint DEGS_all results already exist under 06.results/<Species>/paired/
#   • Filenames contain a date tag (YYYYMMDD), species code, and tXXXX token.
#   • Files include at least: gene_id, log2FoldChange.
#   • If date_tag is NULL, the script accepts any 8-digit date in filenames.
#
# Parameters:
#   species        : Species code string (e.g., "Scer", "Cgla", "Calb", "Klac")
#   results_root   : Root results directory (default: here("06.results"))
#   date_tag       : Specific YYYYMMDD to filter files; NULL = any date
#   out_subdir     : Subfolder under paired/ for combined outputs ("combined")
#   zero_fill      : If TRUE, replace missing LFCs with 0 in outputs
#
# Dependencies:
#   readr, dplyr, tidyr, stringr, purrr, tibble
#
# Usage:
#   # Single species
#   lfc_Scer <- summarize_species_lfc_matrix("Scer",
#     results_root = here::here("06.results"),
#     date_tag     = "20251019",    # or NULL
#     out_subdir   = "combined",
#     zero_fill    = TRUE)
#
#   # All species
#   lfc_all <- summarize_all_species_lfc_matrices(
#     species_vec  = c("Scer","Cgla","Calb","Klac"),
#     results_root = here::here("06.results"),
#     date_tag     = "20251019",
#     out_subdir   = "combined",
#     zero_fill    = TRUE)
# ==============================================================================


summarize_species_lfc_matrix <- function(
    species,
    results_root = here::here("06.results"),
    date_tag = NULL,        # e.g. "20251019"; if NULL, accept any date
    out_subdir = "combined",
    zero_fill = TRUE        # fill NAs with 0
) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  req <- c("readr","dplyr","tidyr","stringr","tibble","purrr","here")
  
  miss <- req[!vapply(req, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(miss)) {
    message("Installing missing packages via BiocManager: ", paste(miss, collapse = ", "))
    BiocManager::install(miss, ask = FALSE, update = FALSE)
  }
  still_miss <- req[!vapply(req, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(still_miss)) {
    stop("These packages could not be installed/loaded: ",
         paste(still_miss, collapse = ", "),
         "\nPlease check your R setup and try again.")
  }
  
  suppressPackageStartupMessages({
    library(readr); library(dplyr); library(tidyr)
    library(stringr); library(purrr); library(tibble); library(here)
  })
  
  sp_dir <- file.path(results_root, species, "paired")
  if (!dir.exists(sp_dir)) stop("Missing directory: ", sp_dir)
  
  # Match DE_full files for this species/date
  date_pat <- if (is.null(date_tag)) "\\d{8}" else stringr::fixed(date_tag)
  pat_full <- paste0("^", date_pat, "__", species, "__paired__t\\d{4}_vs_t0000__DEGs_all__fdr0010_lfc1000\\.csv$")
  
  files <- list.files(sp_dir, pattern = pat_full, full.names = TRUE)
  if (!length(files)) stop("No DE_full files found for ", species, " (pattern: ", pat_full, ")")
  
  # read each DE_full -> (gene_id, tp, log2FC)
  tp_from_name <- function(x) stringr::str_extract(basename(x), "t\\d{4}")
  lfc_long <- purrr::map2_dfr(
    files, purrr::map_chr(files, tp_from_name),
    ~ readr::read_csv(.x, show_col_types = FALSE) |>
      transmute(gene_id, tp = .y, log2FC = log2FoldChange)
  )
  
  # wide: genes x timepoints
  lfc_wide <- lfc_long |>
    tidyr::pivot_wider(names_from = tp, values_from = log2FC) |>
    arrange(gene_id)
  
  # fill missing with 0 if requested
  if (zero_fill) {
    lfc_wide <- lfc_wide |>
      mutate(across(-gene_id, ~ dplyr::coalesce(., 0)))
  }
  
  # also provide tidy long with zero-fill for plotting (optional)
  lfc_long_out <- lfc_long
  if (zero_fill) {
    # complete ensures every gene_id x tp exists, then fill with 0
    all_tps <- sort(unique(lfc_long$tp))
    lfc_long_out <- lfc_long |>
      tidyr::complete(gene_id, tp = all_tps) |>
      mutate(log2FC = dplyr::coalesce(log2FC, 0)) |>
      arrange(gene_id, tp)
  }
  
  # write outputs
  out_dir <- file.path(sp_dir, out_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  tag <- ifelse(is.null(date_tag), "ALLDATES", date_tag)
  readr::write_csv(lfc_wide,     file.path(out_dir, sprintf("%s__%s__paired__LFC_matrix.csv", tag, species)))
  readr::write_csv(lfc_long_out, file.path(out_dir, sprintf("%s__%s__paired__LFC_tidy_long.csv", tag, species)))
  
  invisible(lfc_wide)
}

# Run for all species convenience
summarize_all_species_lfc_matrices <- function(
    species_vec = c("Scer","Cgla","Calb","Klac"),
    results_root = here::here("06.results"),
    date_tag = NULL,
    out_subdir = "combined",
    zero_fill = TRUE
) {
  out <- vector("list", length(species_vec)); names(out) <- species_vec
  for (sp in species_vec) {
    message("Summarizing LFC matrix for ", sp, " ...")
    out[[sp]] <- summarize_species_lfc_matrix(
      species = sp, results_root = results_root, date_tag = date_tag,
      out_subdir = out_subdir, zero_fill = zero_fill
    )
  }
  invisible(out)
}

# One species
#lfc_Scer <- summarize_species_lfc_matrix("Scer", results_root = here("06.results"), date_tag = "20251019")

# All four species
lfc_all <- summarize_all_species_lfc_matrices(
  species_vec = c("Scer","Cgla","Calb","Klac"),
  results_root = here("06.results"),
  date_tag = "20251019"  # or NULL to accept any date
)

