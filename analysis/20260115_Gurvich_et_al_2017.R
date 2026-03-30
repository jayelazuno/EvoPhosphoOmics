---
title: "Gurvich et al. 2017 S. cer"
author: "Joshua Ayelazuno"
date: "2026-01-14"
output: html_document
---

# ============================================================
# Gurvich et al. (2017) goal: define the dynamics and regulatory programs of the
# phosphate-starvation response in S. cerevisiae across graded Pi limitation.
# Experimental design: cells were grown in phosphate-replete media (7.3 mM Pi) and
# shifted to 0, 0.06, 0.2, or 0.5 mM Pi; gene-level log2 fold-changes were reported
# over a common time course (0, 1, 2, 3.5, 5, 6.5, 8, 9.5, 24.67 h).For our
# analysis, we extracted only the phosphate-starvation portion of their dataset, 
# parsed the two-row headers to assign each column to its correct Pi block (including 
# the 0 h baseline within each block), tidied the data into a long format, and then 
# (i) visualized Gurvich “stress genes” directly as heatmaps and (ii) mapped Gurvich 
# gene symbols onto our curated OSR backbone (using SCER gene name first, then ORF ID) 
# to generate comparable OSR-focused heatmaps that highlight overlap vs. missing (NA) measurements.
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE)) BiocManager::install("circlize")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(here)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(gridExtra)
})

#  Helpers functions 
`%||%` <- function(x, y) if (!is.null(x)) x else y

.keep_times_tol <- function(x, times, tol = 1e-6) {
  x <- as.numeric(x)
  keep <- rep(FALSE, length(x))
  for (t in times) keep <- keep | (abs(x - t) < tol)
  keep
}

standardize_function <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- dplyr::if_else(is.na(x) | x == "", "Unknown Function", x)
  stringr::str_to_title(x)
}

# ============================================================
# Inputs / Paths

scer_dir  <- here("06.results", "Scer")
plots_dir <- here("08.plots", "Scer")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

times_keep <- c(0, 1, 2, 3.5, 5, 6.5, 8, 9.5, 24.67)

# OSR colors
osr_category_colors <- c(
  "Antioxidant"           = "#E41A1C",
  "Chaperon"              = "#377EB8",
  "Amino Acid Metabolism" = "#4DAF4A",
  "Carbon Metabolism"     = "#984EA3",
  "Protein Degradation"   = "#FF7F00",
  "Not Classified"        = "#999999",
  "Unknown Function"      = "#A65628"
)

# OSR backbone
osr_backbone_path <- file.path(scer_dir, "20260113_Sc_OSR_genes.xlsx")
if (!file.exists(osr_backbone_path)) stop("OSR backbone not found: ", osr_backbone_path)

osr_backbone <- readxl::read_excel(osr_backbone_path, sheet = 1, col_types = "text") %>%
  transmute(
    gene_id = stringr::str_to_upper(stringr::str_trim(as.character(gene_id))),
    scer_gene_name = stringr::str_to_upper(stringr::str_trim(as.character(scer_gene_name))),
    functional_category = standardize_function(functional_category),
    category_order = suppressWarnings(as.numeric(category_order))
  ) %>%
  mutate(
    scer_gene_name = if_else(is.na(scer_gene_name) | scer_gene_name == "", gene_id, scer_gene_name),
    functional_category = factor(functional_category, levels = names(osr_category_colors))
  ) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(category_order, scer_gene_name)

# ============================================================
# Gurvich et al. 2017 read + tidy (Sheet 1)
# - keep 0h for each Pi block (baseline 7.3mM column is assigned to next Pi)

gurvich_path <- file.path(scer_dir, "Gurvich-2017-S2Data.xlsx")
if (!file.exists(gurvich_path)) stop("Gurvich file not found: ", gurvich_path)

g_raw <- readxl::read_excel(gurvich_path, sheet = 1, col_names = FALSE)

hdr_pi   <- g_raw[1, ] %>% unlist(use.names = FALSE) %>% as.character()
hdr_time <- g_raw[2, ] %>% unlist(use.names = FALSE) %>% as.character()
hdr_time_num <- suppressWarnings(as.numeric(hdr_time))

data_cols <- which(!is.na(hdr_time_num))

pi_tokens <- hdr_pi[data_cols] %>%
  stringr::str_squish() %>%
  dplyr::na_if("") %>%
  as.character()

pi_tokens_fill <- tidyr::fill(tibble(tok = pi_tokens), tok, .direction = "down")$tok

col_meta_all <- tibble(
  col_index = data_cols,
  tok = stringr::str_squish(pi_tokens_fill),
  time_h = as.numeric(hdr_time_num[data_cols])
)

# For each column, assign pi_condition as:
# - if tok != 7.3mM, it's already the target Pi column
# - if tok == 7.3mM and next token is target Pi, treat that 7.3mM (baseline) as belonging to that target Pi block

col_meta <- col_meta_all %>%
  mutate(
    tok_next = dplyr::lead(tok),
    pi_condition = dplyr::case_when(
      !is.na(tok) & tok != "7.3mM" ~ tok,
      tok == "7.3mM" & !is.na(tok_next) & tok_next != "7.3mM" ~ tok_next,
      TRUE ~ NA_character_
    )
  ) %>%
  transmute(
    col_index,
    pi_condition = stringr::str_squish(pi_condition),
    time_h = as.numeric(time_h)
  ) %>%
  filter(!is.na(pi_condition)) %>%
  filter(.keep_times_tol(time_h, times_keep))

# gene rows start at row 3
g_dat <- g_raw[-c(1, 2), , drop = FALSE]

gene_sym <- g_dat[[1]] %>% as.character() %>% stringr::str_trim()

wide_vals <- g_dat[, col_meta$col_index, drop = FALSE] %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.x)))))

tmp <- wide_vals
names(tmp) <- paste0("V", col_meta$col_index)

gurvich_pi_tidy <- tmp %>%
  mutate(gurvich_gene = gene_sym) %>%
  pivot_longer(
    cols = starts_with("V"),
    names_to = "Vcol",
    values_to = "log2FC"
  ) %>%
  mutate(col_index = as.integer(stringr::str_remove(Vcol, "^V"))) %>%
  left_join(col_meta, by = "col_index") %>%
  transmute(
    gurvich_gene = stringr::str_to_upper(stringr::str_trim(as.character(gurvich_gene))),
    pi_condition = stringr::str_squish(pi_condition),
    time_h = as.numeric(time_h),
    time_label = paste0(sprintf("%g", as.numeric(time_h)), "h"),
    log2FC = as.numeric(log2FC)
  ) %>%
  filter(!is.na(pi_condition), !is.na(time_h)) %>%
  filter(.keep_times_tol(time_h, times_keep)) %>%
  filter(!is.na(gurvich_gene) & gurvich_gene != "")

# ============================================================
# Map Gurvich genes to OSR backbone (prefer scer_gene_name; fallback gene_id)

osr_keymap <- bind_rows(
  osr_backbone %>%
    transmute(key = stringr::str_to_upper(stringr::str_trim(scer_gene_name)),
              gene_id, scer_gene_name, functional_category, category_order),
  osr_backbone %>%
    transmute(key = stringr::str_to_upper(stringr::str_trim(gene_id)),
              gene_id, scer_gene_name, functional_category, category_order)
) %>%
  filter(!is.na(key), key != "") %>%
  distinct(key, .keep_all = TRUE)

gurvich_osr_long <- gurvich_pi_tidy %>%
  left_join(osr_keymap, by = c("gurvich_gene" = "key")) %>%
  mutate(functional_category = factor(functional_category, levels = names(osr_category_colors)))

gurvich_osr_long <- gurvich_osr_long %>%
  filter(!is.na(scer_gene_name) & scer_gene_name != "")

# ============================================================
# Heatmap inputs: stress genes (no functional annotation)
# NOTE: NA means "not quantified / missing in Gurvich table" -> grey

make_gurvich_stress_inputs <- function(pi_keep, gurvich_pi_tidy, times_keep) {
  df <- gurvich_pi_tidy %>%
    filter(pi_condition == pi_keep) %>%
    filter(.keep_times_tol(time_h, times_keep)) %>%
    mutate(time_label = paste0(sprintf("%g", time_h), "h"))
  
  time_levels <- df %>%
    distinct(time_h, time_label) %>%
    arrange(time_h) %>%
    pull(time_label)
  
  sum_df <- df %>%
    group_by(gurvich_gene, time_label) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    mutate(time_label = factor(time_label, levels = time_levels))
  
  mat <- sum_df %>%
    pivot_wider(names_from = time_label, values_from = log2FC) %>%
    column_to_rownames("gurvich_gene") %>%
    as.matrix()
  
  mat <- mat[, time_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  
  list(mat = mat, time_levels = time_levels)
}

make_gurvich_stress_heatmap <- function(pi_keep,
                                        gurvich_pi_tidy,
                                        times_keep,
                                        lfc_limits = c(-4, 4),
                                        na_col = "grey85",
                                        row_font = 7,
                                        col_font = 9,
                                        cell_w_mm = 5,
                                        cell_h_mm = 2.8,
                                        row_name_max_mm = 70,
                                        draw_now = TRUE) {
  
  inputs <- make_gurvich_stress_inputs(pi_keep, gurvich_pi_tidy, times_keep)
  mat <- inputs$mat
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "log2FC",
    na_col = na_col,
    col = circlize::colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold"),
    
    width  = grid::unit(ncol(mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_h_mm, "mm"),
    
    column_title = pi_keep,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    
    heatmap_legend_param = list(
      title = "log2FC",
      at = c(lfc_limits[1], 0, lfc_limits[2]),
      title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
    ),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat)
  if (draw_now) ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                                     padding = grid::unit(c(2,2,2,6), "mm"))
  invisible(out)
}

# ============================================================
# Heatmap inputs: OSR genes (with functional annotation)

make_gurvich_osr_inputs <- function(pi_keep, gurvich_osr_long, osr_backbone, times_keep) {
  
  df <- gurvich_osr_long %>%
    filter(pi_condition == pi_keep) %>%
    filter(.keep_times_tol(time_h, times_keep)) %>%
    mutate(time_label = paste0(sprintf("%g", time_h), "h"))
  
  time_levels <- df %>%
    distinct(time_h, time_label) %>%
    arrange(time_h) %>%
    pull(time_label)
  
  sum_df <- df %>%
    group_by(scer_gene_name, time_label) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    mutate(time_label = factor(time_label, levels = time_levels))
  
  backbone2 <- osr_backbone %>%
    arrange(category_order, scer_gene_name) %>%
    mutate(scer_gene_name = make.unique(scer_gene_name)) %>%
    distinct(scer_gene_name, .keep_all = TRUE)
  
  joined <- backbone2 %>%
    left_join(sum_df, by = "scer_gene_name") %>%
    mutate(time_label = factor(as.character(time_label), levels = time_levels))
  
  mat <- joined %>%
    select(scer_gene_name, time_label, log2FC) %>%
    pivot_wider(names_from = time_label, values_from = log2FC) %>%
    column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  mat <- mat[, time_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(rownames(mat))
  
  row_annot <- backbone2 %>%
    filter(scer_gene_name %in% rownames(mat)) %>%
    arrange(match(scer_gene_name, rownames(mat))) %>%
    select(scer_gene_name, functional_category, category_order)
  
  list(mat = mat, row_annot = row_annot, time_levels = time_levels)
}

make_gurvich_osr_heatmap <- function(pi_keep,
                                     gurvich_osr_long,
                                     osr_backbone,
                                     osr_category_colors,
                                     times_keep,
                                     lfc_limits = c(-4, 4),
                                     na_col = "grey85",
                                     row_font = 7,
                                     col_font = 9,
                                     cell_w_mm = 5,
                                     cell_h_mm = 2.8,
                                     row_name_max_mm = 70,
                                     draw_now = TRUE) {
  
  inputs <- make_gurvich_osr_inputs(pi_keep, gurvich_osr_long, osr_backbone, times_keep)
  mat <- inputs$mat
  row_annot <- inputs$row_annot
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = row_annot$functional_category,
    col = list(Function = osr_category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    annotation_legend_param = list(
      Function = list(
        title = "Function",
        title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
      )
    )
  )
  
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "log2FC",
    na_col = na_col,
    col = circlize::colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = row_ha,
    
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold"),
    
    width  = grid::unit(ncol(mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_h_mm, "mm"),
    
    column_title = pi_keep,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    
    heatmap_legend_param = list(
      title = "log2FC",
      at = c(lfc_limits[1], 0, lfc_limits[2]),
      title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
    ),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat)
  if (draw_now) ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                                     annotation_legend_side = "right",
                                     padding = grid::unit(c(2,2,2,6), "mm"))
  invisible(out)
}

# ============================================================
# Overview heatmaps (generic)

make_overview_heatmap_no_row_anno <- function(ht_list,
                                              block_order,
                                              pretty_block,
                                              lfc_limits = c(-4, 4),
                                              na_col = "grey85",
                                              row_font = 7,
                                              col_font = 8,
                                              cell_w_mm = 3.2,
                                              cell_h_mm = 2.6,
                                              row_name_max_mm = 70,
                                              block_label_rot = 45,
                                              block_label_font = 10,
                                              column_gap_mm = 1,
                                              top_padding_mm = 10,
                                              draw_now = TRUE) {
  
  all_genes <- Reduce(union, lapply(ht_list, function(x) rownames(x$mat)))
  all_genes <- sort(all_genes)
  
  block_mats <- list()
  block_labels_for_cols <- character(0)
  col_labels_for_cols <- character(0)
  
  for (bk in block_order) {
    if (!bk %in% names(ht_list)) stop("Missing block in ht_list: ", bk)
    mat <- ht_list[[bk]]$mat
    
    miss <- setdiff(all_genes, rownames(mat))
    if (length(miss)) {
      add <- matrix(NA_real_, nrow = length(miss), ncol = ncol(mat),
                    dimnames = list(miss, colnames(mat)))
      mat <- rbind(mat, add)
    }
    mat <- mat[all_genes, , drop = FALSE]
    
    block_mats[[bk]] <- mat
    block_labels_for_cols <- c(block_labels_for_cols, rep(bk, ncol(mat)))
    col_labels_for_cols   <- c(col_labels_for_cols, colnames(mat))
  }
  
  big_mat <- do.call(cbind, block_mats)
  storage.mode(big_mat) <- "numeric"
  
  block_factor <- factor(block_labels_for_cols, levels = block_order)
  block_levels <- levels(block_factor)
  
  block_pretty <- pretty_block[block_levels]
  block_pretty[is.na(block_pretty) | block_pretty == ""] <- block_levels
  slice_titles <- unname(block_pretty)
  
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
    column_title_gp = grid::gpar(fontsize = block_label_font, fontface = "bold"),
    column_title_rot = block_label_rot,
    
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    show_column_names = TRUE,
    column_labels = col_labels_for_cols,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold", rot = 90),
    
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
  if (draw_now) ComplexHeatmap::draw(ht, heatmap_legend_side = "right",
                                     padding = grid::unit(c(top_padding_mm,2,2,6), "mm"))
  invisible(out)
}

make_overview_heatmap_with_row_anno <- function(ht_list,
                                                block_order,
                                                pretty_block,
                                                backbone,
                                                category_colors,
                                                lfc_limits = c(-4, 4),
                                                na_col = "grey85",
                                                row_font = 7,
                                                col_font = 8,
                                                cell_w_mm = 3.2,
                                                cell_h_mm = 2.6,
                                                row_name_max_mm = 70,
                                                block_label_rot = 45,
                                                block_label_font = 10,
                                                column_gap_mm = 1,
                                                top_padding_mm = 10,
                                                draw_now = TRUE) {
  
  backbone2 <- backbone %>%
    arrange(category_order, scer_gene_name) %>%
    mutate(scer_gene_name = make.unique(scer_gene_name)) %>%
    distinct(scer_gene_name, .keep_all = TRUE)
  
  row_order <- backbone2$scer_gene_name
  
  block_mats <- list()
  block_labels_for_cols <- character(0)
  col_labels_for_cols <- character(0)
  
  for (bk in block_order) {
    if (!bk %in% names(ht_list)) stop("Missing block in ht_list: ", bk)
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
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = backbone2$functional_category,
    col = list(Function = category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 9, fontface = "bold"),
    annotation_legend_param = list(
      Function = list(
        title = "Function",
        title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
        labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
      )
    )
  )
  
  block_factor <- factor(block_labels_for_cols, levels = block_order)
  block_levels <- levels(block_factor)
  
  block_pretty <- pretty_block[block_levels]
  block_pretty[is.na(block_pretty) | block_pretty == ""] <- block_levels
  slice_titles <- unname(block_pretty)
  
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
    
    left_annotation = row_ha,
    
    column_split = block_factor,
    column_gap   = grid::unit(column_gap_mm, "mm"),
    
    column_title = slice_titles,
    column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = block_label_font, fontface = "bold"),
    column_title_rot = block_label_rot,
    
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    show_column_names = TRUE,
    column_labels = col_labels_for_cols,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold", rot = 90),
    
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
  if (draw_now) ComplexHeatmap::draw(ht,
                                     heatmap_legend_side = "right",
                                     annotation_legend_side = "right",
                                     merge_legends = TRUE,
                                     padding = grid::unit(c(top_padding_mm,2,2,6), "mm"))
  invisible(out)
}

# ============================================================
# Build  plots and SAVE PDFs

gurvich_pi_conditions <- c("0mM Pi", "0.06mM Pi", "0.2mM Pi", "0.5mM Pi")
gurvich_pi_conditions <- gurvich_pi_conditions[gurvich_pi_conditions %in% unique(gurvich_pi_tidy$pi_condition)]

gurvich_pretty_block <- c(
  "0mM Pi"    = "0 mM Pi",
  "0.06mM Pi" = "0.06 mM Pi",
  "0.2mM Pi"  = "0.2 mM Pi",
  "0.5mM Pi"  = "0.5 mM Pi"
)

#- Plot 1: stress genes only 
gurvich_stress_ht_list <- setNames(vector("list", length(gurvich_pi_conditions)), gurvich_pi_conditions)
for (pc in gurvich_pi_conditions) {
  gurvich_stress_ht_list[[pc]] <- make_gurvich_stress_heatmap(
    pi_keep = pc,
    gurvich_pi_tidy = gurvich_pi_tidy,
    times_keep = times_keep,
    draw_now = FALSE
  )
}

gurvich_stress_overview <- make_overview_heatmap_no_row_anno(
  ht_list = gurvich_stress_ht_list,
  block_order = gurvich_pi_conditions,
  pretty_block = gurvich_pretty_block,
  column_gap_mm = 1,
  top_padding_mm = 12,
  draw_now = FALSE
)

na_leg <- ComplexHeatmap::Legend(
  labels = "NA",
  title  = NULL,
  legend_gp = grid::gpar(fill = "grey85", col = NA)
)

out_pdf1 <- file.path(plots_dir, "Gurvich2017_PiGradient_stressGenes_heatmap.pdf")
pdf(out_pdf1, width = 8, height = 7)
ComplexHeatmap::draw(
  gurvich_stress_overview$ht,
  heatmap_legend_side = "right",
  heatmap_legend_list = na_leg,
  padding = grid::unit(c(12,2,2,6), "mm")
)
dev.off()

# Plot 2: overlap-only OSR genes
overlap_genes <- gurvich_osr_long %>% distinct(scer_gene_name) %>% pull(scer_gene_name)

osr_backbone_overlap <- osr_backbone %>%
  filter(scer_gene_name %in% overlap_genes) %>%
  arrange(category_order, scer_gene_name)

gurvich_ht_list_overlap <- setNames(vector("list", length(gurvich_pi_conditions)), gurvich_pi_conditions)
for (pc in gurvich_pi_conditions) {
  gurvich_ht_list_overlap[[pc]] <- make_gurvich_osr_heatmap(
    pi_keep = pc,
    gurvich_osr_long = gurvich_osr_long,
    osr_backbone = osr_backbone_overlap,
    osr_category_colors = osr_category_colors,
    times_keep = times_keep,
    draw_now = FALSE
  )
}

gurvich_overview_overlap <- make_overview_heatmap_with_row_anno(
  ht_list = gurvich_ht_list_overlap,
  block_order = gurvich_pi_conditions,
  pretty_block = gurvich_pretty_block,
  backbone = osr_backbone_overlap,
  category_colors = osr_category_colors,
  column_gap_mm = 1,
  top_padding_mm = 12,
  draw_now = FALSE
)

out_pdf2 <- file.path(plots_dir, "Gurvich2017_PiGradient_OSR_overlapOnly_heatmap.pdf")
pdf(out_pdf2, width = 8, height = 4)
ComplexHeatmap::draw(
  gurvich_overview_overlap$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
  padding = grid::unit(c(12,2,2,6), "mm")
)
dev.off()

#Plot 3: full OSR backbone + explicit Missing legend -
gurvich_ht_list_full <- setNames(vector("list", length(gurvich_pi_conditions)), gurvich_pi_conditions)
for (pc in gurvich_pi_conditions) {
  gurvich_ht_list_full[[pc]] <- make_gurvich_osr_heatmap(
    pi_keep = pc,
    gurvich_osr_long = gurvich_osr_long,
    osr_backbone = osr_backbone,
    osr_category_colors = osr_category_colors,
    times_keep = times_keep,
    draw_now = FALSE
  )
}

gurvich_overview_full <- make_overview_heatmap_with_row_anno(
  ht_list = gurvich_ht_list_full,
  block_order = gurvich_pi_conditions,
  pretty_block = gurvich_pretty_block,
  backbone = osr_backbone,
  category_colors = osr_category_colors,
  column_gap_mm = 1,
  top_padding_mm = 12,
  draw_now = FALSE
)

missing_leg <- ComplexHeatmap::Legend(
  labels = "Missing",
  title  = NULL,
  legend_gp = grid::gpar(fill = "grey85", col = NA)
)

out_pdf3 <- file.path(plots_dir, "Gurvich2017_PiGradient_OSR_full_heatmap.pdf")
pdf(out_pdf3, width = 8, height = 8)
ComplexHeatmap::draw(
  gurvich_overview_full$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  merge_legends = TRUE,
  heatmap_legend_list = missing_leg,
  padding = grid::unit(c(12,2,2,6), "mm")
)
dev.off()


