---
title: "Characterizing the evolutionary divergence oxidative stress response (OSR) structure and dynamics across four yeast species under -Pi stress"
author: "Joshua Ayelazuno"
date: "2025-11-11"
output: html_document
---
# ============================================================
# ESR pipeline (Gasch et al. ESR clusters) -> 4-species heatmaps
# - Parses NAME into: gene_id (ORF), gene_name (fallback ORF), description
# - Collapses diverse descriptions into major FUNCTION GROUPS (small legend)
# - Uses the SAME pipeline as OSR:
#     (1) read LFC tidy for all species
#     (2) standardize to 10 timepoints per species
#     (3) orthology map ESR genes (Scer->others) + choose best member
#     (4) write mapped CSVs per species, read back, build heatmaps
#     (5) save single-species PDFs + 4-species overview PDF

# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("patchwork", quietly = TRUE)) BiocManager::install("patchwork")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(readxl)
  library(here)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)      # gpar()
  library(patchwork)
  library(openxlsx)
})



`%||%` <- function(x, y) if (!is.null(x)) x else y

# ============================================================
# Paths

scer_dir  <- here("06.results", "Scer")
plots_dir <- here("08.plots", "Scer")

# Gasch ESR file (sheet2)
esr_path  <- file.path(scer_dir, "20250109_ESR_clusters_UPDATED_2017.xlsx")
esr_sheet <- 2
if (!file.exists(esr_path)) stop("ESR file not found: ", esr_path)

# Heatmap limits
lfc_limits <- c(-6, 6)
na_col     <- "grey85"

# Layout for ~800 rows
single_cell_w_mm <- 3.2
single_cell_h_mm <- 0.2
single_row_font  <- 5
single_col_font  <- 8

overview_cell_w_mm <- 2.6
overview_cell_h_mm <- 0.2
overview_col_font  <- 7
overview_block_rot <- 45
overview_show_row_names <- FALSE

# ============================================================
# Timepoint standardization (10 points per species)

.tp_sort_key <- function(tp) {
  tp_chr <- as.character(tp)
  num <- suppressWarnings(as.numeric(str_extract(tp_chr, "\\d+\\.?\\d*")))
  tibble(tp_chr = tp_chr, num = num) %>%
    mutate(ord_num = ifelse(is.na(num), Inf, num)) %>%
    arrange(ord_num, tp_chr) %>%
    pull(tp_chr)
}

make_tp_map_10 <- function(lfc_tidy_df, species_code) {
  stopifnot(all(c("gene_id", "tp", "log2FC") %in% names(lfc_tidy_df)))
  
  tp_chr <- unique(as.character(lfc_tidy_df$tp))
  
  # Your rule: Calb/Klac missing two late points
  if (species_code %in% c("Calb", "Klac")) {
    tp_chr <- tp_chr[!tp_chr %in% c("t0150", "t0210")]
  }
  
  tp_sorted <- .tp_sort_key(tp_chr)
  tp_keep <- tp_sorted[seq_len(min(10, length(tp_sorted)))]
  
  tibble(
    tp = tp_keep,
    time_std = paste0("T", sprintf("%02d", seq_along(tp_keep)))
  )
}

apply_tp_map <- function(lfc_tidy_df, tp_map) {
  lfc_tidy_df %>%
    mutate(
      gene_id = str_trim(as.character(gene_id)),
      tp      = str_trim(as.character(tp)),
      log2FC  = suppressWarnings(as.numeric(log2FC))
    ) %>%
    left_join(tp_map %>% mutate(tp = str_trim(as.character(tp))), by = "tp") %>%
    filter(!is.na(time_std)) %>%
    mutate(time_std = factor(time_std, levels = paste0("T", sprintf("%02d", 1:10))))
}

# Preserve species prefix (Cg-/Ca-/Kl-) but uppercase the rest
upper_preserve_prefix <- function(x) {
  x <- str_trim(as.character(x))
  m <- str_match(x, "^(Cg-|Ca-|Kl-)?(.*)$")
  pref <- m[, 2]
  rest <- m[, 3]
  paste0(ifelse(is.na(pref), "", pref), str_to_upper(rest))
}

# ============================================================
# Legends / label cleanup

make_na_legend <- function(label = "No ortholog (NA)", col = "grey85") {
  ComplexHeatmap::Legend(
    labels = label,
    legend_gp = grid::gpar(fill = col, col = NA),
    title = NULL
  )
}

# ============================================================
# ESR parsing + functional category collapse (ESR-specific)

parse_esr_NAME <- function(name_vec) {
  x <- str_squish(as.character(name_vec))
  
  gene_id <- str_extract(x, "^[A-Z0-9]{2,3}\\d{3}[A-Z]")  # YAL061W etc.
  rest1   <- str_trim(str_remove(x, "^[A-Z0-9]{2,3}\\d{3}[A-Z]\\s*"))
  
  gene_name_raw <- str_extract(rest1, "^[A-Za-z0-9_\\-]+")
  gene_name <- if_else(
    is.na(gene_name_raw) | gene_name_raw == "" | str_to_upper(gene_name_raw) == "UNKNOWN",
    NA_character_, gene_name_raw
  )
  
  rest2 <- if_else(
    !is.na(gene_name_raw) & gene_name_raw != "",
    str_trim(str_remove(rest1, paste0("^", fixed(gene_name_raw), "\\s*"))),
    rest1
  )
  
  tibble(
    gene_id   = str_to_upper(str_trim(gene_id)),
    gene_name = str_to_upper(str_trim(gene_name %||% "")),
    desc_raw  = rest2
  ) %>%
    mutate(gene_name = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name))
}

collapse_esr_function <- function(desc_raw) {
  d <- str_to_lower(str_squish(as.character(desc_raw)))
  
  case_when(
    str_detect(d, "ribosom|translation|trna|rrna|initiation|elongation") ~ "Translation / Ribosome",
    str_detect(d, "heat shock|hsp|chaperon|fold|proteas|ubiquitin|degrad") ~ "Protein folding / QC",
    str_detect(d, "mitochond|respirat|oxidative phosphorylation|atp synth") ~ "Mitochondria / Respiration",
    str_detect(d, "glycol|gluconeo|trehalose|carbon|metabolism|tca|ferment") ~ "Carbon metabolism",
    str_detect(d, "amino acid|nitrogen|biosynth|purine|pyrimid|nucleotide") ~ "Biosynthesis (AA/Nuc)",
    str_detect(d, "cell wall|membrane|ergosterol|lipid|transport|secret|vesicle|er\\b|golgi") ~ "Membrane / Trafficking",
    str_detect(d, "dna|repair|recomb|chromatin|histone|replication|cell cycle|mitosis") ~ "Genome / Cell cycle",
    str_detect(d, "autophag|vacuol|lysosom|peroxisom") ~ "Autophagy / Vacuole",
    str_detect(d, "stress|oxid|detox|redox|peroxid|glutathione|thioredoxin") ~ "Stress / Redox",
    str_detect(d, "transcription|rna polymer|splicing|mrna|rna binding") ~ "Gene expression (RNA)",
    str_detect(d, "\\bunknown\\b|uncharacterized|hypothetical") ~ "Unknown",
    TRUE ~ "Other"
  )
}

esr_category_order_tbl <- tribble(
  ~functional_category,             ~category_order,
  "Stress / Redox",                 1,
  "Protein folding / QC",           2,
  "Translation / Ribosome",         3,
  "Gene expression (RNA)",          4,
  "Genome / Cell cycle",            5,
  "Membrane / Trafficking",         6,
  "Mitochondria / Respiration",     7,
  "Carbon metabolism",              8,
  "Biosynthesis (AA/Nuc)",          9,
  "Autophagy / Vacuole",           10,
  "Other",                         11,
  "Unknown",                       12
)

esr_category_colors <- c(
  "Stress / Redox"             = "#E41A1C",
  "Protein folding / QC"       = "#377EB8",
  "Translation / Ribosome"     = "#4DAF4A",
  "Gene expression (RNA)"      = "#984EA3",
  "Genome / Cell cycle"        = "#FF7F00",
  "Membrane / Trafficking"     = "#A65628",
  "Mitochondria / Respiration" = "#F781BF",
  "Carbon metabolism"          = "#999999",
  "Biosynthesis (AA/Nuc)"      = "#66C2A5",
  "Autophagy / Vacuole"        = "#FC8D62",
  "Other"                      = "#8DA0CB",
  "Unknown"                    = "#BDBDBD"
)

# ============================================================
# Build complete matrix (ALL genes x ALL times)
# IMPORTANT: forces gene labels to UPPERCASE (prefix preserved)

make_complete_hm_matrix <- function(gene_df, lfc_std,
                                    time_levels = paste0("T", sprintf("%02d", 1:10)),
                                    keep_na_if_gene_id_missing = FALSE) {
  
  gene_df <- gene_df %>%
    mutate(
      gene_id        = na_if(str_trim(as.character(gene_id)), ""),
      scer_gene_name = upper_preserve_prefix(str_trim(as.character(scer_gene_name))),
      functional_category = str_trim(as.character(functional_category)),
      category_order = suppressWarnings(as.numeric(category_order))
    ) %>%
    mutate(
      functional_category = if_else(is.na(functional_category) | functional_category == "", "Other", functional_category),
      category_order = if_else(is.na(category_order), 99, category_order)
    ) %>%
    distinct(scer_gene_name, .keep_all = TRUE) %>%
    arrange(category_order, scer_gene_name)
  
  complete_grid <- tidyr::expand_grid(
    scer_gene_name = gene_df$scer_gene_name,
    time_std = factor(time_levels, levels = time_levels)
  )
  
  observed <- lfc_std %>%
    mutate(
      gene_id  = na_if(str_trim(as.character(gene_id)), ""),
      time_std = factor(as.character(time_std), levels = time_levels),
      log2FC   = suppressWarnings(as.numeric(log2FC))
    ) %>%
    select(gene_id, time_std, log2FC)
  
  id_to_name <- gene_df %>% select(gene_id, scer_gene_name) %>% distinct()
  
  complete_long <- complete_grid %>%
    left_join(id_to_name, by = "scer_gene_name") %>%
    left_join(observed, by = c("gene_id", "time_std"))
  
  if (keep_na_if_gene_id_missing) {
    complete_long <- complete_long %>%
      mutate(log2FC = if_else(is.na(gene_id), as.numeric(NA), tidyr::replace_na(log2FC, 0)))
  } else {
    complete_long <- complete_long %>%
      mutate(log2FC = tidyr::replace_na(log2FC, 0))
  }
  
  mat <- complete_long %>%
    select(scer_gene_name, time_std, log2FC) %>%
    pivot_wider(names_from = time_std, values_from = log2FC) %>%
    column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat <- mat[gene_df$scer_gene_name, , drop = FALSE]
  
  row_annot <- gene_df %>%
    filter(scer_gene_name %in% rownames(mat)) %>%
    arrange(match(scer_gene_name, rownames(mat))) %>%
    select(scer_gene_name, functional_category, category_order)
  
  list(mat = mat, row_annot = row_annot)
}

# ============================================================
# Heatmap plotter (row_labels allows prefixes in single-species plots)

plot_esr_heatmap <- function(mat, row_annot, category_colors,
                             row_labels = NULL,
                             show_row_names = TRUE,
                             lfc_limits = c(-6, 6),
                             heatmap_name = "log2FC",
                             na_col = "grey85",
                             cell_w_mm = 3.2,
                             cell_h_mm = 1.6,
                             row_name_max_mm = 55,
                             row_font = 5,
                             col_font = 8,
                             legend_font = 9) {
  
  desired_w_mm <- ncol(mat) * cell_w_mm
  desired_h_mm <- nrow(mat) * cell_h_mm
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = row_annot$functional_category,
    col = list(Function = category_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Function = list(
        title = "Function",
        title_gp  = grid::gpar(fontsize = legend_font + 1, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = legend_font, fontface = "bold")
      )
    )
  )
  
  if (is.null(row_labels)) row_labels <- rownames(mat)
  
  list(
    ht = ComplexHeatmap::Heatmap(
      mat,
      name = heatmap_name,
      na_col = na_col,
      col = circlize::colorRamp2(
        c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
        c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
      ),
      width  = grid::unit(desired_w_mm, "mm"),
      height = grid::unit(desired_h_mm, "mm"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      
      show_row_names = show_row_names,
      row_labels = row_labels,
      row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
      row_names_max_width = grid::unit(row_name_max_mm, "mm"),
      
      show_column_names = TRUE,
      column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold"),
      
      left_annotation = row_ha,
      
      heatmap_legend_param = list(
        title = heatmap_name,
        at = c(lfc_limits[1], 0, lfc_limits[2]),
        title_gp  = grid::gpar(fontsize = legend_font + 1, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = legend_font, fontface = "bold")
      ),
      border = FALSE,
      column_title = NULL,
      row_title = NULL
    ),
    na_legend = make_na_legend(col = na_col)
  )
}

# ============================================================
# Overview heatmap (4 blocks)

make_species_overview_heatmap <- function(ht_list,
                                          block_order,
                                          pretty_block_expr,
                                          category_colors,
                                          show_row_names = TRUE,
                                          lfc_limits = c(-6, 6),
                                          na_col = "grey85",
                                          row_font = 7,
                                          col_font = 8,
                                          cell_w_mm = 3.2,
                                          cell_h_mm = 2.6,
                                          row_name_max_mm = 70,
                                          block_label_rot = 45,
                                          block_label_font = 11,
                                          column_gap_mm = 2,
                                          top_padding_mm = 10,
                                          draw_now = TRUE) {
  
  ref_block <- block_order[1]
  row_order <- rownames(ht_list[[ref_block]]$mat)
  
  block_mats <- list()
  block_labels_for_cols <- character(0)
  col_labels_for_cols <- character(0)
  
  for (bk in block_order) {
    mat <- ht_list[[bk]]$mat
    
    miss <- setdiff(row_order, rownames(mat))
    if (length(miss)) {
      add <- matrix(NA_real_, nrow = length(miss), ncol = ncol(mat),
                    dimnames = list(miss, colnames(mat)))
      mat <- rbind(mat, add)
    }
    mat <- mat[row_order, , drop = FALSE]
    
    block_mats[[bk]] <- mat
    block_labels_for_cols <- c(block_labels_for_cols, rep(bk, ncol(mat)))
    col_labels_for_cols   <- c(col_labels_for_cols, colnames(mat))
  }
  
  big_mat <- do.call(cbind, block_mats)
  storage.mode(big_mat) <- "numeric"
  
  block_factor <- factor(block_labels_for_cols, levels = block_order)
  slice_titles <- pretty_block_expr[levels(block_factor)]
  slice_titles <- as.expression(slice_titles)
  
  # Functional annotation from reference block (already aligned)
  ref_row_annot <- ht_list[[ref_block]]$row_annot %>%
    select(scer_gene_name, functional_category) %>%
    distinct()
  
  ref_row_annot <- ref_row_annot %>%
    filter(scer_gene_name %in% row_order) %>%
    slice(match(row_order, scer_gene_name))
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = ref_row_annot$functional_category,
    col = list(Function = category_colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(
      Function = list(
        title = "Function",
        title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
      )
    )
  )
  
  ht <- ComplexHeatmap::Heatmap(
    big_mat,
    name = "log2FC",
    na_col = na_col,
    col = circlize::colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    column_split = block_factor,
    column_gap   = grid::unit(column_gap_mm, "mm"),
    
    column_title = slice_titles,
    column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = block_label_font),
    column_title_rot = block_label_rot,
    
    show_row_names = show_row_names,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    show_column_names = TRUE,
    column_labels = col_labels_for_cols,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold", rot = 90),
    
    width  = grid::unit(ncol(big_mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(big_mat) * cell_h_mm, "mm"),
    
    left_annotation = row_ha,
    
    heatmap_legend_param = list(
      title = "log2FC",
      at = c(lfc_limits[1], 0, lfc_limits[2]),
      title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
    ),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = big_mat)
  
  if (draw_now) {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(top_padding_mm, 2, 2, 6), "mm")
    )
  }
  
  invisible(out)
}

# ============================================================
# Read ESR (Gasch sheet2) + build Scer ESR backbone

esr_raw <- suppressMessages(
  readxl::read_excel(esr_path, sheet = esr_sheet, col_types = "text", .name_repair = "minimal")
)
if (!("NAME" %in% names(esr_raw))) stop("Expected a column named 'NAME' in sheet ", esr_sheet)

name_parsed <- parse_esr_NAME(esr_raw$NAME)

esr_backbone <- name_parsed %>%
  transmute(
    gene_id,
    scer_gene_name = gene_name,
    desc_raw
  ) %>%
  mutate(functional_category = collapse_esr_function(desc_raw)) %>%
  left_join(esr_category_order_tbl, by = "functional_category") %>%
  mutate(
    functional_category = if_else(is.na(functional_category) | functional_category == "", "Other", functional_category),
    category_order      = if_else(is.na(category_order), 99, category_order),
    scer_gene_name      = upper_preserve_prefix(scer_gene_name),
    gene_id             = str_to_upper(str_trim(gene_id))
  ) %>%
  distinct(scer_gene_name, .keep_all = TRUE) %>%
  arrange(category_order, scer_gene_name)

readr::write_csv(
  esr_backbone %>% select(gene_id, scer_gene_name, functional_category, category_order, desc_raw),
  file.path(scer_dir, "ESR_genes_Scer_backbone_parsed.csv")
)

# ============================================================
# Read LFC tidy for all species + build tp maps

read_lfc_species <- function(sp_code) {
  f <- list.files(here("06.results", sp_code, "paired", "combined"),
                  pattern = "LFC_tidy", full.names = TRUE)
  if (length(f) == 0) stop("No LFC_tidy file found for: ", sp_code)
  readr::read_csv(f[1], show_col_types = FALSE) %>%
    mutate(
      gene_id = str_trim(as.character(gene_id)),
      tp      = str_trim(as.character(tp)),
      log2FC  = suppressWarnings(as.numeric(log2FC))
    )
}

lfc_by_species <- list(
  Scer = read_lfc_species("Scer"),
  Cgla = read_lfc_species("Cgla"),
  Calb = read_lfc_species("Calb"),
  Klac = read_lfc_species("Klac")
)

tp_maps <- list(
  Scer = make_tp_map_10(lfc_by_species$Scer, "Scer"),
  Cgla = make_tp_map_10(lfc_by_species$Cgla, "Cgla"),
  Calb = make_tp_map_10(lfc_by_species$Calb, "Calb"),
  Klac = make_tp_map_10(lfc_by_species$Klac, "Klac")
)

# ============================================================
# Orthology mapping + robust best-member choice (OSR logic)

orthofinder_dir <- here("05.metadata", "Orthogroups")
orthogroups <- readr::read_tsv(file.path(orthofinder_dir, "Orthogroups.tsv"),
                               show_col_types = FALSE)

required_cols <- c("Orthogroup", "S_cerevisiae", "C_glabrata", "C_albicans", "K_lactis")
if (!all(required_cols %in% names(orthogroups))) {
  stop("Orthogroups.tsv missing required columns: ",
       paste(setdiff(required_cols, names(orthogroups)), collapse = ", "))
}

# ---- Make orthogroups long by Scer gene id: exact join, no regex ----
orth_long_scer <- orthogroups %>%
  select(Orthogroup, S_cerevisiae, C_glabrata, C_albicans, K_lactis) %>%
  tidyr::separate_rows(S_cerevisiae, sep = ",\\s*") %>%
  mutate(
    S_cerevisiae = str_extract(S_cerevisiae, "^[^\\s|]+"),
    S_cerevisiae = str_trim(S_cerevisiae)
  ) %>%
  filter(!is.na(S_cerevisiae), S_cerevisiae != "") %>%
  distinct(Orthogroup, S_cerevisiae, .keep_all = TRUE)

# Join ESR Scer genes to orthogroups (0/1 row per Scer gene id; no many-to-many warning)
esr_orthologs <- esr_backbone %>%
  transmute(
    scer_gene_id   = gene_id,
    scer_gene_name = scer_gene_name,
    functional_category,
    category_order
  ) %>%
  distinct(scer_gene_id, .keep_all = TRUE) %>%
  left_join(
    orth_long_scer %>%
      transmute(
        Orthogroup,
        scer_gene_id = S_cerevisiae,
        Scer_genes = S_cerevisiae,
        Cgla_genes = C_glabrata,
        Calb_genes = C_albicans,
        Klac_genes = K_lactis
      ),
    by = "scer_gene_id"
  ) %>%
  select(
    Orthogroup, scer_gene_name, scer_gene_id,
    functional_category, category_order,
    Scer_genes, Cgla_genes, Calb_genes, Klac_genes
  )

esr_orthologs_long <- esr_orthologs %>%
  pivot_longer(cols = ends_with("_genes"), names_to = "species", values_to = "gene_ids") %>%
  mutate(
    species = str_remove(species, "_genes"),
    species = case_when(
      species == "Scer" ~ "S_cerevisiae",
      species == "Cgla" ~ "C_glabrata",
      species == "Calb" ~ "C_albicans",
      species == "Klac" ~ "K_lactis",
      TRUE ~ species
    )
  ) %>%
  separate_rows(gene_ids, sep = ",\\s*") %>%
  mutate(
    gene_ids = if_else(!is.na(gene_ids) & gene_ids != "",
                       str_extract(gene_ids, "^[^\\s|]+"),
                       NA_character_)
  )

# ID remaps (existing files)
xp_klla_map <- readr::read_csv(here("06.results", "Klac", "xp_to_klla0_mapping.csv"),
                               show_col_types = FALSE)
qng_cagl_map <- readr::read_csv(here("06.results", "Cgla", "qng_gwk_cagl_threeway_map_complete.csv"),
                                show_col_types = FALSE)

esr_orthologs_long <- esr_orthologs_long %>%
  mutate(original_id = gene_ids) %>%
  left_join(xp_klla_map %>% select(xp_id, locus_tag), by = c("gene_ids" = "xp_id")) %>%
  left_join(qng_cagl_map %>% select(qng_id, gwk60_id, cagl_id), by = c("gene_ids" = "qng_id")) %>%
  mutate(
    gene_ids = case_when(
      species == "K_lactis"   & str_detect(original_id, "^XP_") & !is.na(locus_tag) ~ locus_tag,
      species == "C_glabrata" & str_detect(original_id, "^QNG") & !is.na(gwk60_id)  ~ gwk60_id,
      TRUE ~ gene_ids
    )
  ) %>%
  select(-locus_tag, -gwk60_id)

# ---- Robust best-member choice: peak ABS(|log2FC|) ----
choose_best_member <- function(lfc_std, candidate_ids) {
  candidate_ids <- unique(na.omit(str_trim(as.character(candidate_ids))))
  candidate_ids <- candidate_ids[candidate_ids != ""]
  if (length(candidate_ids) == 0) return(NA_character_)
  
  lfc_std <- lfc_std %>% mutate(gene_id = str_trim(as.character(gene_id)))
  
  peaks <- lfc_std %>%
    filter(gene_id %in% candidate_ids) %>%
    group_by(gene_id) %>%
    summarise(
      peak_abs = suppressWarnings(max(abs(log2FC), na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    mutate(peak_abs = ifelse(is.infinite(peak_abs), NA_real_, peak_abs))
  
  if (nrow(peaks) == 0 || all(is.na(peaks$peak_abs))) return(candidate_ids[1])
  peaks %>% arrange(desc(peak_abs)) %>% slice(1) %>% pull(gene_id)
}

build_best_member_mapped_table <- function(orth_long_df, lfc_std_by_species, species_full) {
  
  df_sp <- orth_long_df %>%
    filter(species == species_full) %>%
    group_by(Orthogroup, scer_gene_name, scer_gene_id, functional_category, category_order) %>%
    summarise(
      candidates   = list(unique(stats::na.omit(gene_ids))),
      original_ids = list(unique(stats::na.omit(original_id))),
      cagl_ids     = list(unique(stats::na.omit(cagl_id))),
      .groups = "drop"
    )
  
  df_sp %>%
    mutate(
      gene_id = purrr::map_chr(candidates, function(x) choose_best_member(lfc_std_by_species, x)),
      original_id = purrr::map_chr(original_ids, function(x) if (length(x) == 0) NA_character_ else x[[1]]),
      cagl_id     = purrr::map_chr(cagl_ids,     function(x) if (length(x) == 0) NA_character_ else x[[1]])
    ) %>%
    select(
      Orthogroup, gene_id, original_id, cagl_id,
      scer_gene_name, scer_gene_id, functional_category, category_order
    ) %>%
    arrange(category_order, scer_gene_name)
}

# ============================================================
# WRITE mapped files (Cgla/Calb/Klac) + symmetry file for Scer

for (sp_full in c("C_glabrata", "C_albicans", "K_lactis")) {
  
  species_code <- dplyr::case_when(
    sp_full == "C_glabrata" ~ "Cgla",
    sp_full == "C_albicans" ~ "Calb",
    sp_full == "K_lactis"   ~ "Klac"
  )
  
  lfc_std <- apply_tp_map(lfc_by_species[[species_code]], tp_maps[[species_code]])
  
  best_tbl <- build_best_member_mapped_table(
    orth_long_df = esr_orthologs_long,
    lfc_std_by_species = lfc_std,
    species_full = sp_full
  )
  
  # Keep ALL ESR rows even if no orthogroup/ortholog
  species_genes <- esr_backbone %>%
    transmute(scer_gene_name, functional_category, category_order) %>%
    distinct(scer_gene_name, .keep_all = TRUE) %>%
    left_join(best_tbl, by = c("scer_gene_name", "functional_category", "category_order")) %>%
    transmute(
      Orthogroup,
      gene_id,
      original_id,
      cagl_id,
      scer_gene_name,
      functional_category,
      category_order
    ) %>%
    arrange(category_order, scer_gene_name)
  
  out_csv <- here("06.results", species_code, paste0("ESR_genes_", species_code, "_mapped_best.csv"))
  readr::write_csv(species_genes, out_csv)
  message("Wrote: ", out_csv)
}

# Scer symmetry file
readr::write_csv(
  esr_backbone %>% transmute(
    Orthogroup = NA_character_,
    gene_id,
    original_id = NA_character_,
    cagl_id = NA_character_,
    scer_gene_name,
    functional_category,
    category_order
  ),
  here("06.results", "Scer", "ESR_genes_Scer_mapped_best.csv")
)

# ============================================================
# Build per-species heatmap objects (mat + ht)
# IMPORTANT: overview must share SAME rownames across blocks (Scer names, no prefixes).
# For single species, show prefix via row_labels only.

make_species_esr_ht_obj <- function(species_code) {
  
  prefix_map <- c(Scer = "", Cgla = "Cg-", Calb = "Ca-", Klac = "Kl-")
  prefix <- prefix_map[[species_code]] %||% ""
  
  if (species_code == "Scer") {
    gene_df <- esr_backbone %>%
      transmute(gene_id, scer_gene_name, functional_category, category_order)
  } else {
    mapped_path <- here("06.results", species_code, paste0("ESR_genes_", species_code, "_mapped_best.csv"))
    gene_df <- readr::read_csv(mapped_path, show_col_types = FALSE) %>%
      transmute(
        gene_id = str_trim(as.character(gene_id)),
        scer_gene_name = str_trim(as.character(scer_gene_name)),
        functional_category = str_trim(as.character(functional_category)),
        category_order = suppressWarnings(as.numeric(category_order))
      )
  }
  
  # stabilize categories
  gene_df <- gene_df %>%
    mutate(
      functional_category = if_else(is.na(functional_category) | functional_category == "", "Other", functional_category),
      functional_category = factor(functional_category, levels = esr_category_order_tbl$functional_category),
      category_order = if_else(is.na(category_order), 99, category_order),
      scer_gene_name = upper_preserve_prefix(scer_gene_name)
    ) %>%
    arrange(category_order, scer_gene_name) %>%
    mutate(functional_category = as.character(functional_category))
  
  lfc_std <- apply_tp_map(lfc_by_species[[species_code]], tp_maps[[species_code]])
  observed <- lfc_std %>% select(gene_id, time_std, log2FC)
  
  inputs <- make_complete_hm_matrix(
    gene_df = gene_df,
    lfc_std = observed,
    time_levels = as.character(tp_maps[[species_code]]$time_std),
    keep_na_if_gene_id_missing = (species_code != "Scer")
  )
  
  row_labels <- if (species_code == "Scer") {
    rownames(inputs$mat)
  } else {
    upper_preserve_prefix(paste0(prefix, rownames(inputs$mat)))
  }
  
  p <- plot_esr_heatmap(
    mat = inputs$mat,
    row_annot = inputs$row_annot,
    category_colors = esr_category_colors,
    row_labels = row_labels,
    show_row_names = FALSE,
    lfc_limits = lfc_limits,
    na_col = na_col,
    cell_w_mm = single_cell_w_mm,
    cell_h_mm = single_cell_h_mm,
    row_font  = single_row_font,
    col_font  = single_col_font
  )
  
  list(
    ht = p$ht,
    na_legend = p$na_legend,
    mat = inputs$mat,
    row_annot = inputs$row_annot
  )
}

scer_ht <- make_species_esr_ht_obj("Scer")
cgla_ht <- make_species_esr_ht_obj("Cgla")
calb_ht <- make_species_esr_ht_obj("Calb")
klac_ht <- make_species_esr_ht_obj("Klac")

# ============================================================
# 4-species OVERVIEW heatmap

esr_species_ht_list <- list(Scer = scer_ht, Cgla = cgla_ht, Calb = calb_ht, Klac = klac_ht)
block_order <- c("Scer", "Cgla", "Calb", "Klac")

pretty_block_expr <- list(
  Scer = bquote(bolditalic("S. cerevisiae")),
  Cgla = bquote(bolditalic("C. glabrata")),
  Calb = bquote(bolditalic("C. albicans")),
  Klac = bquote(bolditalic("K. lactis"))
)

overview <- make_species_overview_heatmap(
  ht_list = esr_species_ht_list,
  block_order = block_order,
  pretty_block_expr = pretty_block_expr,
  category_colors = esr_category_colors,
  show_row_names = overview_show_row_names,
  lfc_limits = lfc_limits,
  na_col = na_col,
  col_font = overview_col_font,
  cell_w_mm = overview_cell_w_mm,
  cell_h_mm = overview_cell_h_mm,
  block_label_rot = overview_block_rot,
  column_gap_mm = 1,
  draw_now = FALSE
)

# ============================================================
# SAVE PLOTS

na_leg <- make_na_legend(label = "No ortholog (NA)", col = na_col)

# Single-species PDFs
pdf(file.path(plots_dir, "ESR_Scer_heatmap.pdf"), width = 5, height = 9)
ComplexHeatmap::draw(scer_ht$ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf(file.path(plots_dir, "ESR_Cgla_heatmap.pdf"), width = 5, height = 9)
ComplexHeatmap::draw(
  cgla_ht$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(cgla_ht$na_legend),
  merge_legends = TRUE
)
dev.off()

pdf(file.path(plots_dir, "ESR_Calb_heatmap.pdf"), width = 5, height = 9)
ComplexHeatmap::draw(
  calb_ht$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(calb_ht$na_legend),
  merge_legends = TRUE
)
dev.off()

pdf(file.path(plots_dir, "ESR_Klac_heatmap.pdf"), width = 5, height = 9)
ComplexHeatmap::draw(
  klac_ht$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(klac_ht$na_legend),
  merge_legends = TRUE
)
dev.off()

# Combined overview
pdf(file.path(plots_dir, "ESR_4species_overview_heatmap.pdf"), width = 7.5, height = 9)
ComplexHeatmap::draw(
  overview$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(na_leg),
  merge_legends = TRUE,
  padding = grid::unit(c(12, 2, 2, 6), "mm")
)
dev.off()

