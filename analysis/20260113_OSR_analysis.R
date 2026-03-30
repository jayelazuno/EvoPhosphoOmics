---
title: "Characterizing the evolutionary divergence oxidative stress response (OSR) structure and dynamics across four yeast species under -Pi stress"
author: "Joshua Ayelazuno"
date: "2025-11-11"
output: html_document
---
  
  # ============================================================
# - Read OSR gene set (Scer), map to gene_id via 20250910_Scer_gene_list.csv
# - Add category_order (given), keep requested columns
# - Overwrite OSR file as Excel (same base name)
# - Make Scer heatmap
# - Orthogroup mapping + XP/QNG conversions
# - If expanded family: choose "best" member (peak log2FC) for heatmap
# - Make Cgla/Calb/Klac heatmaps + 2x2 shared legend draw
# ============================================================

# -------------------- Packages ------------------------------
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

standardize_function <- function(x) {
  x <- str_trim(as.character(x))
  x <- if_else(is.na(x) | x == "", "Unknown Function", x)
  str_to_title(x)
}

# ============================================================
# Build complete matrix (ALL genes x ALL times)
# IMPORTANT: forces gene labels to UPPERCASE (prefix preserved)

make_complete_hm_matrix <- function(gene_df, lfc_std,
                                    time_levels = paste0("T", sprintf("%02d", 1:10)),
                                    keep_na_if_gene_id_missing = FALSE) {
  
  gene_df <- gene_df %>%
    mutate(
      gene_id        = na_if(str_trim(as.character(gene_id)), ""),
      scer_gene_name = upper_preserve_prefix(str_trim(as.character(scer_gene_name))),  # <-- UPPERCASE
      functional_category = standardize_function(functional_category),
      category_order = suppressWarnings(as.numeric(category_order))
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
  
  # enforce row order
  mat <- mat[gene_df$scer_gene_name, , drop = FALSE]
  
  row_annot <- gene_df %>%
    filter(scer_gene_name %in% rownames(mat)) %>%
    arrange(match(scer_gene_name, rownames(mat))) %>%
    select(scer_gene_name, functional_category, category_order)
  
  list(mat = mat, row_annot = row_annot)
}

# ============================================================
# Heatmap plotter (row_labels lets you show prefixes in single-species plots)

plot_osr_heatmap <- function(mat, row_annot, category_colors,
                             row_labels = NULL,
                             lfc_limits = c(-6, 6),
                             heatmap_name = "log2FC",
                             na_col = "grey85",
                             cell_w_mm = 6,
                             cell_h_mm = 3,
                             row_name_max_mm = 70,
                             row_font = 7,
                             col_font = 9,
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
        c(lfc_limits[1], lfc_limits[1] / 2, 0, lfc_limits[2] / 2, lfc_limits[2]),
        c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
      ),
      width  = grid::unit(desired_w_mm, "mm"),
      height = grid::unit(desired_h_mm, "mm"),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      
      show_row_names = TRUE,
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
# Overview heatmap (Gasch/Gurvich style)
# - species titles are plotmath expressions (bold+italic) at 45°
# ============================================================
make_osr_species_overview_heatmap <- function(ht_list,
                                              block_order,
                                              pretty_block_expr,
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
    
    # enforce same row order; missing rows become NA
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
  
  # ---- NEW: pull functional categories from reference block + add left annotation ----
  ref_row_annot <- ht_list[[ref_block]]$row_annot %>%
    mutate(
      scer_gene_name = upper_preserve_prefix(str_trim(as.character(scer_gene_name))),
      functional_category = standardize_function(functional_category)
    ) %>%
    distinct(scer_gene_name, .keep_all = TRUE) %>%
    filter(scer_gene_name %in% row_order) %>%
    arrange(match(scer_gene_name, row_order))
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = ref_row_annot$functional_category,
    col = list(Function = osr_category_colors),
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
    
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    show_column_names = TRUE,
    column_labels = col_labels_for_cols,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold", rot = 90),
    
    left_annotation = row_ha,   
    
    
    width  = grid::unit(ncol(big_mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(big_mat) * cell_h_mm, "mm"),
    
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
      padding = grid::unit(c(top_padding_mm, 2, 2, 6), "mm")
    )
  }
  
  invisible(out)
}

# ============================================================
# Paths

scer_dir  <- here("06.results", "Scer")
plots_dir <- here("08.plots", "Scer")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ============================================================
# OSR category order + colors

osr_category_order <- tribble(
  ~functional_category,        ~category_order,
  "Antioxidant",               1,
  "Chaperon",                  2,
  "Amino Acid Metabolism",     3,
  "Carbon Metabolism",         4,
  "Protein Degradation",       5,
  "Not Classified",            6,
  "Unknown Function",          7
)

osr_category_colors <- c(
  "Antioxidant"           = "#E41A1C",
  "Chaperon"              = "#377EB8",
  "Amino Acid Metabolism" = "#4DAF4A",
  "Carbon Metabolism"     = "#984EA3",
  "Protein Degradation"   = "#FF7F00",
  "Not Classified"        = "#999999",
  "Unknown Function"      = "#A65628"
)

# ============================================================
# Read OSR gene set + map to Scer gene_id

guess_osr_path <- function(scer_dir, base_name = "20260113_Sc_OSR_genes") {
  candidates <- c(
    file.path(scer_dir, base_name),
    file.path(scer_dir, paste0(base_name, ".xlsx")),
    file.path(scer_dir, paste0(base_name, ".csv")),
    file.path(scer_dir, paste0(base_name, ".tsv"))
  )
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit)) stop("OSR file not found. Looked for: ", paste(candidates, collapse = ", "))
  hit
}

read_osr_gene_set <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("xlsx", "xls")) {
    df <- readxl::read_excel(path, sheet = 1, col_types = "text")
  } else if (ext == "csv") {
    df <- readr::read_csv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  } else if (ext %in% c("tsv", "txt")) {
    df <- readr::read_tsv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  } else {
    df <- readr::read_csv(path, show_col_types = FALSE, col_types = cols(.default = "c"))
  }
  
  df %>%
    transmute(
      originalName = str_trim(as.character(originalName)),
      gene_name = str_trim(as.character(gene_name)),
      functional_category = standardize_function(functional_category),
      reference = str_trim(as.character(reference)),
      comment   = str_trim(as.character(comment))
    ) %>%
    left_join(osr_category_order, by = "functional_category") %>%
    mutate(category_order = as.numeric(category_order))
}

read_scer_gene_list <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
    transmute(
      gene_id   = str_trim(as.character(gene_id)),
      gene_name = str_trim(as.character(gene_name))
    ) %>%
    mutate(
      gene_name = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name),
      gene_key  = str_to_upper(gene_name)
    ) %>%
    distinct(gene_key, .keep_all = TRUE)
}

map_osr_to_gene_ids <- function(osr_df, scer_gene_list_df) {
  scer_keys <- scer_gene_list_df %>% select(gene_id, gene_name_ref = gene_name, gene_key)
  
  by_gene_name <- osr_df %>%
    mutate(gene_key = str_to_upper(str_trim(as.character(gene_name)))) %>%
    left_join(scer_keys, by = "gene_key") %>%
    mutate(mapped_by = if_else(!is.na(gene_id), "gene_name", NA_character_))
  
  out <- by_gene_name %>%
    mutate(orig_key = str_to_upper(str_trim(as.character(originalName)))) %>%
    left_join(
      scer_keys %>% rename(gene_id_orig = gene_id, gene_name_orig = gene_name_ref, orig_key = gene_key),
      by = "orig_key"
    ) %>%
    mutate(
      gene_id = if_else(is.na(gene_id) & !is.na(gene_id_orig), gene_id_orig, gene_id),
      mapped_by = if_else(!is.na(mapped_by), mapped_by,
                          if_else(!is.na(gene_id_orig), "originalName", NA_character_))
    ) %>%
    select(-gene_id_orig, -gene_name_orig, -gene_key, -orig_key)
  
  out %>%
    mutate(
      scer_gene_name = if_else(!is.na(gene_name) & gene_name != "", gene_name, originalName),
      scer_gene_name = str_trim(as.character(scer_gene_name))
    ) %>%
    distinct(scer_gene_name, .keep_all = TRUE) %>%
    arrange(category_order, scer_gene_name)
}

osr_path <- guess_osr_path(scer_dir, base_name = "20260113_Sc_OSR_genes1")
osr_raw  <- read_osr_gene_set(osr_path)

scer_gene_list <- read_scer_gene_list(here("06.results", "Scer", "20250910_Scer_gene_list.csv"))

osr_mapped <- map_osr_to_gene_ids(osr_raw, scer_gene_list) %>%
  transmute(
    gene_id,
    scer_gene_name,
    functional_category,
    category_order,
    reference,
    comment,
    mapped_by
  )
#write.csv(osr_mapped, file.path(scer_dir, "OSR_Scer_mapped.csv"), row.names = FALSE)


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
# Orthology mapping + best-member choice (Cgla/Calb/Klac)

orthofinder_dir <- here("05.metadata", "Orthogroups")
orthogroups <- readr::read_tsv(file.path(orthofinder_dir, "Orthogroups.tsv"),
                               show_col_types = FALSE)

required_cols <- c("Orthogroup", "S_cerevisiae", "C_glabrata", "C_albicans", "K_lactis")
if (!all(required_cols %in% names(orthogroups))) {
  stop("Orthogroups.tsv missing required columns: ",
       paste(setdiff(required_cols, names(orthogroups)), collapse = ", "))
}

scer_ids_regex <- paste(na.omit(osr_mapped$gene_id), collapse = "|")

# Build a long table of orthogroup memberships for OSR genes
osr_orthologs <- osr_mapped %>%
  transmute(
    scer_gene_id   = gene_id,
    scer_gene_name = scer_gene_name,
    functional_category,
    category_order
  ) %>%
  left_join(
    orthogroups %>%
      filter(str_detect(S_cerevisiae, scer_ids_regex)) %>%
      mutate(
        scer_gene_id_hit = map_chr(S_cerevisiae, function(genes) {
          matches <- str_extract_all(genes, scer_ids_regex)[[1]]
          if (length(matches) > 0) matches[1] else NA_character_
        })
      ) %>%
      select(
        Orthogroup,
        scer_gene_id_hit,
        Scer_genes = S_cerevisiae,
        Cgla_genes = C_glabrata,
        Calb_genes = C_albicans,
        Klac_genes = K_lactis
      ),
    by = c("scer_gene_id" = "scer_gene_id_hit")
  ) %>%
  select(
    Orthogroup, scer_gene_name, scer_gene_id,
    functional_category, category_order,
    Scer_genes, Cgla_genes, Calb_genes, Klac_genes
  )

osr_orthologs_long <- osr_orthologs %>%
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

# ID remaps (your existing files)
xp_klla_map <- readr::read_csv(here("06.results", "Klac", "xp_to_klla0_mapping.csv"),
                               show_col_types = FALSE)
qng_cagl_map <- readr::read_csv(here("06.results", "Cgla", "qng_gwk_cagl_threeway_map_complete.csv"),
                                show_col_types = FALSE)

osr_orthologs_long <- osr_orthologs_long %>%
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

choose_best_member <- function(lfc_std, candidate_ids) {
  candidate_ids <- unique(na.omit(str_trim(as.character(candidate_ids))))
  candidate_ids <- candidate_ids[candidate_ids != ""]
  if (length(candidate_ids) == 0) return(NA_character_)
  
  lfc_std <- lfc_std %>% mutate(gene_id = str_trim(as.character(gene_id)))
  
  peaks <- lfc_std %>%
    filter(gene_id %in% candidate_ids) %>%
    group_by(gene_id) %>%
    summarise(peak = suppressWarnings(max(log2FC, na.rm = TRUE)), .groups = "drop") %>%
    mutate(peak = ifelse(is.infinite(peak), NA_real_, peak))
  
  if (nrow(peaks) == 0 || all(is.na(peaks$peak))) return(candidate_ids[1])
  peaks %>% arrange(desc(peak)) %>% slice(1) %>% pull(gene_id)
}

build_best_member_mapped_table <- function(orth_long_df, lfc_std_by_species, species_full) {
  
  df_sp <- orth_long_df %>%
    dplyr::filter(species == species_full) %>%
    dplyr::group_by(Orthogroup, scer_gene_name, scer_gene_id, functional_category, category_order) %>%
    dplyr::summarise(
      candidates   = list(unique(stats::na.omit(gene_ids))),
      original_ids = list(unique(stats::na.omit(original_id))),
      cagl_ids     = list(unique(stats::na.omit(cagl_id))),
      .groups = "drop"
    )
  
  out <- df_sp %>%
    dplyr::mutate(
      gene_id = purrr::map_chr(candidates, function(x) choose_best_member(lfc_std_by_species, x)),
      original_id = purrr::map_chr(original_ids, function(x) if (length(x) == 0) NA_character_ else x[[1]]),
      cagl_id     = purrr::map_chr(cagl_ids,     function(x) if (length(x) == 0) NA_character_ else x[[1]])
    ) %>%
    dplyr::select(
      Orthogroup, gene_id, original_id, cagl_id,
      scer_gene_name, scer_gene_id, functional_category, category_order
    ) %>%
    dplyr::arrange(category_order, scer_gene_name)
  
  out
}

# ============================================================
# WRITE mapped files and read them back for plotting

for (sp_full in c("C_glabrata", "C_albicans", "K_lactis")) {
  
  species_code <- dplyr::case_when(
    sp_full == "C_glabrata" ~ "Cgla",
    sp_full == "C_albicans" ~ "Calb",
    sp_full == "K_lactis"   ~ "Klac"
  )
  
  lfc_std <- apply_tp_map(lfc_by_species[[species_code]], tp_maps[[species_code]])
  
  best_tbl <- build_best_member_mapped_table(
    orth_long_df = osr_orthologs_long,
    lfc_std_by_species = lfc_std,
    species_full = sp_full
  )
  
  # Keep ALL OSR rows even if no orthogroup/ortholog
  species_genes <- osr_mapped %>%
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
  
  out_csv <- here("06.results", species_code, paste0("OSR_genes_", species_code, "_mapped_best.csv"))
  readr::write_csv(species_genes, out_csv)
  message("Wrote: ", out_csv)
}

# ============================================================
# Build per-species heatmap objects (mat + ht)
# IMPORTANT: overview must share SAME rownames across blocks (Scer gene names, no prefixes).
# For single species, show prefix via row_labels only.

make_species_osr_ht_obj <- function(species_code) {
  
  prefix_map <- c(Scer = "", Cgla = "Cg-", Calb = "Ca-", Klac = "Kl-")
  prefix <- prefix_map[[species_code]] %||% ""
  
  if (species_code == "Scer") {
    gene_df <- osr_mapped %>%
      transmute(gene_id, scer_gene_name, functional_category, category_order)
  } else {
    mapped_path <- here("06.results", species_code, paste0("OSR_genes_", species_code, "_mapped_best.csv"))
    gene_df <- readr::read_csv(mapped_path, show_col_types = FALSE) %>%
      transmute(
        gene_id = str_trim(as.character(gene_id)),
        scer_gene_name = str_trim(as.character(scer_gene_name)),
        functional_category = standardize_function(functional_category),
        category_order = suppressWarnings(as.numeric(category_order))
      )
  }
  
  lfc_std <- apply_tp_map(lfc_by_species[[species_code]], tp_maps[[species_code]])
  observed <- lfc_std %>% select(gene_id, time_std, log2FC)
  
  inputs <- make_complete_hm_matrix(
    gene_df = gene_df,
    lfc_std = observed,
    time_levels = as.character(tp_maps[[species_code]]$time_std),
    keep_na_if_gene_id_missing = (species_code != "Scer")
  )
  
  # Show prefixes in single-species plots (uppercase after prefix)
  row_labels <- if (species_code == "Scer") {
    rownames(inputs$mat)
  } else {
    upper_preserve_prefix(paste0(prefix, rownames(inputs$mat)))
  }
  
  p <- plot_osr_heatmap(
    inputs$mat,
    inputs$row_annot,
    osr_category_colors,
    row_labels = row_labels
  )
  
  list(
    ht = p$ht,
    na_legend = p$na_legend,
    mat = inputs$mat,
    row_annot = inputs$row_annot
  )
}

scer_ht <- make_species_osr_ht_obj("Scer")
cgla_ht <- make_species_osr_ht_obj("Cgla")
calb_ht <- make_species_osr_ht_obj("Calb")
klac_ht <- make_species_osr_ht_obj("Klac")

# ============================================================
# 4-species OVERVIEW heatmap (titles bold+italic at 45°)

osr_species_ht_list <- list(Scer = scer_ht, Cgla = cgla_ht, Calb = calb_ht, Klac = klac_ht)
block_order <- c("Scer", "Cgla", "Calb", "Klac")

pretty_block_expr <- list(
  Scer = bquote(bolditalic("S. cerevisiae")),
  Cgla = bquote(bolditalic("C. glabrata")),
  Calb = bquote(bolditalic("C. albicans")),
  Klac = bquote(bolditalic("K. lactis"))
)

osr_overview <- make_osr_species_overview_heatmap(
  ht_list = osr_species_ht_list,
  block_order = block_order,
  pretty_block_expr = pretty_block_expr,
  column_gap_mm = 1,
  block_label_rot = 45,
  draw_now = FALSE
)

# ============================================================
# SAVE PLOTS

# Single-species
pdf(file.path(plots_dir, "OSR_Scer_heatmap.pdf"), width = 8.5, height = 10)
ComplexHeatmap::draw(scer_ht$ht, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

pdf(file.path(plots_dir, "OSR_Cgla_heatmap.pdf"), width = 8.5, height = 10)
ComplexHeatmap::draw(cgla_ht$ht, heatmap_legend_side = "right", annotation_legend_side = "right",
                     heatmap_legend_list = list(cgla_ht$na_legend), merge_legends = TRUE)
dev.off()

pdf(file.path(plots_dir, "OSR_Calb_heatmap.pdf"), width = 8.5, height = 10)
ComplexHeatmap::draw(calb_ht$ht, heatmap_legend_side = "right", annotation_legend_side = "right",
                     heatmap_legend_list = list(calb_ht$na_legend), merge_legends = TRUE)
dev.off()

pdf(file.path(plots_dir, "OSR_Klac_heatmap.pdf"), width = 8.5, height = 10)
ComplexHeatmap::draw(klac_ht$ht, heatmap_legend_side = "right", annotation_legend_side = "right",
                     heatmap_legend_list = list(klac_ht$na_legend), merge_legends = TRUE)
dev.off()

# Combined overview
na_leg <- make_na_legend(label = "No ortholog (NA)", col = "grey85")

pdf(file.path(plots_dir, "OSR_4species_overview_heatmap.pdf"), width = 13, height = 10)
ComplexHeatmap::draw(
  osr_overview$ht,
  heatmap_legend_side = "right",
  heatmap_legend_list = na_leg,
  merge_legends = TRUE,
  padding = grid::unit(c(12, 2, 2, 6), "mm")
)
dev.off()




# ============================================================
# Panel B (OSR): Induction time is defined as the time point at which an OSR gene
# reaches its maximum log2FC (peak). This compares temporal execution order
# independent of response magnitude.
# Induction time = PEAK log2FC across T01–T10
# ============================================================

cgla_inputs <- make_complete_hm_matrix(
  gene_df = read_csv(
    here("06.results", "Cgla", "OSR_genes_Cgla_mapped_best.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      gene_id = str_trim(as.character(gene_id)),
      scer_gene_name = paste0("Cg-", scer_gene_name),
      functional_category = standardize_function(functional_category),
      category_order = as.numeric(category_order)
    ),
  lfc_std = apply_tp_map(lfc_by_species$Cgla, tp_maps$Cgla),
  time_levels = tp_maps$Cgla$time_std %>% as.character(),
  keep_na_if_gene_id_missing = TRUE
)

calb_inputs <- make_complete_hm_matrix(
  gene_df = read_csv(
    here("06.results", "Calb", "OSR_genes_Calb_mapped_best.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      gene_id = str_trim(as.character(gene_id)),
      scer_gene_name = paste0("Ca-", scer_gene_name),
      functional_category = standardize_function(functional_category),
      category_order = as.numeric(category_order)
    ),
  lfc_std = apply_tp_map(lfc_by_species$Calb, tp_maps$Calb),
  time_levels = tp_maps$Calb$time_std %>% as.character(),
  keep_na_if_gene_id_missing = TRUE
)

klac_inputs <- make_complete_hm_matrix(
  gene_df = read_csv(
    here("06.results", "Klac", "OSR_genes_Klac_mapped_best.csv"),
    show_col_types = FALSE
  ) %>%
    transmute(
      gene_id = str_trim(as.character(gene_id)),
      scer_gene_name = paste0("Kl-", scer_gene_name),
      functional_category = standardize_function(functional_category),
      category_order = as.numeric(category_order)
    ),
  lfc_std = apply_tp_map(lfc_by_species$Klac, tp_maps$Klac),
  time_levels = tp_maps$Klac$time_std %>% as.character(),
  keep_na_if_gene_id_missing = TRUE
)

#  helpers 
# remove species prefixes added for Panel A visualization
strip_species_prefix <- function(x) {
  sub("^(Cg-|Ca-|Kl-)", "", x)
}

# ensure T01..T10 order
order_time_cols <- function(mat) {
  tnum <- suppressWarnings(as.integer(sub("^T", "", colnames(mat))))
  mat[, order(tnum), drop = FALSE]
}

# compute PEAK induction time (index)
compute_peak_time <- function(mat) {
  mat <- order_time_cols(mat)
  
  mat_num <- apply(mat, 2, as.numeric)
  rownames(mat_num) <- rownames(mat)
  tp <- colnames(mat_num)
  
  idx <- apply(mat_num, 1, function(x) {
    if (all(is.na(x))) return(NA_integer_)
    m <- max(x, na.rm = TRUE)
    which(x == m)[1]   # earliest peak if tied
  })
  
  tibble(
    gene = strip_species_prefix(rownames(mat_num)),
    peak_tp  = ifelse(is.na(idx), NA_character_, tp[idx]),
    peak_idx = idx
  )
}

# rank genes by peak time (earliest = rank 1)
rank_peak_time <- function(tbl) {
  tbl %>%
    mutate(
      rank = if_else(
        is.na(peak_idx),
        NA_real_,
        rank(peak_idx, ties.method = "min")
      )
    )
}

# pairwise Spearman rank correlation
pairwise_rank_cor <- function(rank_tbl, sp1, sp2, min_genes = 8) {
  x <- rank_tbl %>% filter(species == sp1) %>% select(gene, r1 = rank)
  y <- rank_tbl %>% filter(species == sp2) %>% select(gene, r2 = rank)
  
  m <- inner_join(x, y, by = "gene") %>%
    filter(!is.na(r1) & !is.na(r2))
  
  if (nrow(m) < min_genes) {
    return(tibble(
      sp1 = sp1, sp2 = sp2,
      n_genes = nrow(m),
      spearman_rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  ct <- suppressWarnings(cor.test(m$r1, m$r2, method = "spearman", exact = FALSE))
  
  tibble(
    sp1 = sp1, sp2 = sp2,
    n_genes = nrow(m),
    spearman_rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

# Panel A OSR matrices 
# MUST be the same matrices used for Panel A OSR heatmaps
mats_osr <- list(
  Scer = scer_inputs$mat,
  Cgla = cgla_inputs$mat,
  Calb = calb_inputs$mat,
  Klac = klac_inputs$mat
)

# sanity check
sapply(mats_osr, nrow)

#  compute peak-time ranks 
rank_tbl_osr <- purrr::imap_dfr(mats_osr, function(mat, sp) {
  compute_peak_time(mat) %>%
    mutate(species = sp) %>%
    rank_peak_time()
})

# Inspect earliest OSR genes per species (optional)
rank_tbl_osr %>%
  filter(!is.na(rank)) %>%
  arrange(species, rank) %>%
  group_by(species) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  print(n = 40)

#  pairwise correlations
species <- names(mats_osr)

pair_tbl_osr <- purrr::map_dfr(
  combn(species, 2, simplify = FALSE),
  \(p) pairwise_rank_cor(rank_tbl_osr, p[1], p[2])
)

pair_tbl_osr %>% arrange(desc(spearman_rho)) %>% print(n = Inf)

#  correlation heatmap 
cor_mat_osr <- pair_tbl_osr %>%
  select(sp1, sp2, spearman_rho) %>%
  bind_rows(pair_tbl_osr %>% transmute(sp1 = sp2, sp2 = sp1, spearman_rho)) %>%
  bind_rows(tibble(sp1 = species, sp2 = species, spearman_rho = 1)) %>%
  distinct() %>%
  mutate(
    sp1 = factor(sp1, levels = species),
    sp2 = factor(sp2, levels = species)
  )

ggplot(cor_mat_osr, aes(sp1, sp2, fill = spearman_rho)) +
  geom_tile(color = "white") +
  coord_equal() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman \u03c1 (OSR genes)"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    axis.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box.margin = margin(b = -25.5)
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      barwidth = unit(10, "cm"),
      barheight = unit(0.4, "cm")
    )
  )

# ============================================================
# Rank–rank scatter plots (OSR)
# computes peak-time ranks (earliest peak = 1)
# ============================================================

# ---- italic species labels (plotmath) ----
species_title_expr <- function(other) {
  switch(
    other,
    "Cgla" = expression(italic("S. cerevisiae")~"vs"~italic("C. glabrata")),
    "Calb" = expression(italic("S. cerevisiae")~"vs"~italic("C. albicans")),
    "Klac" = expression(italic("S. cerevisiae")~"vs"~italic("K. lactis")),
    expression(italic("S. cerevisiae")~"vs"~other)
  )
}

species_y_expr <- function(other) {
  switch(
    other,
    "Cgla" = expression(bold("Peak-time rank in " * italic("C. glabrata"))),
    "Calb" = expression(bold("Peak-time rank in " * italic("C. albicans"))),
    "Klac" = expression(bold("Peak-time rank in " * italic("K. lactis"))),
    expression(bold("Peak-time rank"))
  )
}

make_rankrank_df <- function(rank_tbl, other_species, x_species = "Scer") {
  x <- rank_tbl %>%
    filter(species == x_species) %>%
    select(gene, x_rank = rank)
  
  y <- rank_tbl %>%
    filter(species == other_species) %>%
    select(gene, y_rank = rank)
  
  inner_join(x, y, by = "gene") %>%
    drop_na(x_rank, y_rank)
}

plot_rankrank <- function(rank_tbl, other_species, point_size = 2.4) {
  
  df <- make_rankrank_df(rank_tbl, other_species)
  
  ct <- suppressWarnings(cor.test(df$x_rank, df$y_rank, method = "spearman", exact = FALSE))
  
  rho <- unname(ct$estimate)
  pvl <- ct$p.value
  
  ggplot(df, aes(x = x_rank, y = y_rank)) +
    geom_point(size = point_size) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    coord_equal() +
    labs(
      title = species_title_expr(other_species),
      subtitle = paste0("Spearman \u03c1 = ", sprintf("%.2f", rho),
                        "   (n = ", nrow(df), ", p = ", signif(pvl, 2), ")"),
      x = expression(bold("Peak-time rank in " * italic("S. cerevisiae"))),
      y = species_y_expr(other_species)
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text  = element_text(face = "bold", colour = "black"),
      plot.title = element_text(face = "bold", colour = "black"),
      plot.subtitle = element_text(face = "bold", colour = "black")
    )
}

# ---- OSR plots (Scer vs each species) ----
pC_cgla_osr <- plot_rankrank(rank_tbl_osr, "Cgla")
pC_calb_osr <- plot_rankrank(rank_tbl_osr, "Calb")
pC_klac_osr <- plot_rankrank(rank_tbl_osr, "Klac")

(pC_cgla_osr | pC_calb_osr | pC_klac_osr)

# ============================================================
# Directionality of rewiring (OSR): Δrank vs Scer
# Δrank = rank_other - rank_Scer
# Positive = later than Scer, Negative = earlier than Scer
# ============================================================

make_delta_rank <- function(rank_tbl, other_species, ref_species = "Scer") {
  
  ref <- rank_tbl %>%
    filter(species == ref_species) %>%
    select(gene, rank_ref = rank)
  
  oth <- rank_tbl %>%
    filter(species == other_species) %>%
    select(gene, rank_other = rank)
  
  inner_join(ref, oth, by = "gene") %>%
    drop_na(rank_ref, rank_other) %>%
    mutate(
      other = other_species,
      delta_rank = rank_other - rank_ref,
      abs_shift  = abs(delta_rank)
    )
}

plot_delta_rank <- function(delta_df, other_species) {
  
  delta_df <- delta_df %>%
    arrange(delta_rank) %>%
    mutate(gene = factor(gene, levels = gene))
  
  title_expr <- switch(
    other_species,
    "Cgla" = expression(Delta~"rank (" * italic("C. glabrata") - italic("S. cerevisiae") * ")"),
    "Calb" = expression(Delta~"rank (" * italic("C. albicans") - italic("S. cerevisiae") * ")"),
    "Klac" = expression(Delta~"rank (" * italic("K. lactis") - italic("S. cerevisiae") * ")"),
    expression(Delta~"rank (other - S. cerevisiae)")
  )
  
  ggplot(delta_df, aes(x = gene, y = delta_rank)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_col() +
    coord_flip() +
    labs(
      title = title_expr,
      x = NULL,
      y = expression(bold(Delta~"peak-time rank"))
    ) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.y = element_text(face = "bold", colour = "black"),
      axis.text.x = element_text(face = "bold", colour = "black"),
      axis.title  = element_text(face = "bold", colour = "black"),
      plot.title  = element_text(face = "bold", colour = "black")
    )
}

# ---- build Δrank dataframes (OSR) ----
d_cgla_osr <- make_delta_rank(rank_tbl_osr, "Cgla")
d_calb_osr <- make_delta_rank(rank_tbl_osr, "Calb")
d_klac_osr <- make_delta_rank(rank_tbl_osr, "Klac")

# ---- call individually ----
pD_cgla_osr <- plot_delta_rank(d_cgla_osr, "Cgla")
pD_calb_osr <- plot_delta_rank(d_calb_osr, "Calb")
pD_klac_osr <- plot_delta_rank(d_klac_osr, "Klac")

(pD_cgla_osr | pD_calb_osr | pD_klac_osr)
