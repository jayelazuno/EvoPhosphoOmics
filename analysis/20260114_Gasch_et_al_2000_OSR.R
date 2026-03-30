---
title: "Gasch et al. 2000 S. cer OSR "
author: "Joshua Ayelazuno"
date: "2026-01-14"
output: html_document
---
  
# ============================================================
# # Gasch et al. (2000) goal: build a comprehensive, time-resolved atlas of the
# environmental stress response (ESR) by measuring genome-wide transcriptional
# dynamics across many acute stresses and nutrient perturbations in S. cerevisiae.
# Experimental design: cultures were exposed to diverse stress conditions (e.g.,
# heat shock, oxidants, osmotic shock, and starvation cues) and transcript changes
# were measured across stress-specific time courses, reported as processed log2
# expression changes relative to matched baselines.
# What we did here: parsed the processed microarray table, extracted condition + time metadata
# from column headers, converted all measurements into a tidy long format, and subset to our
# curated OSR backbone using systematic ORF IDs (UID). We then generated (i) per-condition
# OSR heatmaps and (ii) an overview heatmap that concatenates selected stresses side-by-side
# with a consistent color scale and OSR functional annotation for cross-condition comparison.
# ============================================================
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE)) BiocManager::install("circlize")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
if (!requireNamespace("scales", quietly = TRUE)) install.packages("scales")
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

#  Helpers 
`%||%` <- function(x, y) if (!is.null(x)) x else y

standardize_function <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- dplyr::if_else(is.na(x) | x == "", "Unknown Function", x)
  stringr::str_to_title(x)
}

.clean_col_text <- function(x) {
  x %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("030inutes", "030 minutes") %>%  
    stringr::str_replace_all("\\(\\s*", "(") %>%
    stringr::str_replace_all("\\s*\\)", ")")
}

.standardize_time_unit <- function(u) {
  u <- stringr::str_to_lower(stringr::str_trim(u))
  dplyr::case_when(
    u %in% c("min", "mins", "minute", "minutes", "min.") ~ "min",
    u %in% c("h", "hr", "hrs", "hour", "hours") ~ "h",
    u %in% c("d", "day", "days") ~ "d",
    TRUE ~ NA_character_
  )
}

.time_to_minutes <- function(v, u) {
  u <- .standardize_time_unit(u)
  dplyr::case_when(
    is.na(v) | is.na(u) ~ NA_real_,
    u == "min" ~ as.numeric(v),
    u == "h"   ~ as.numeric(v) * 60,
    u == "d"   ~ as.numeric(v) * 1440,
    TRUE ~ NA_real_
  )
}

.clean_condition <- function(x) {
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all("redo|rescan", "") %>%
    stringr::str_replace_all("\\*.*?\\)", ")") %>%       
    stringr::str_replace_all("\\(.*?problem.*?\\)", "") %>%
    stringr::str_replace_all("reference pool", "") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()
}

# Inputs 
scer_dir <- here("06.results", "Scer")

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

# Read OSR backbone file that will be used to subset OSR genes
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

#  Read Gasch sheet2 
gasch_path <- file.path(scer_dir, "20260109_Gasch2000_ProcessedLog2Data.xlsx")
if (!file.exists(gasch_path)) stop("Gasch file not found: ", gasch_path)

gasch_raw <- readxl::read_excel(gasch_path, sheet = 2, col_types = "text")

#Clean Gasch -> tidy long 
gasch_df <- gasch_raw %>%
  rename(
    gene_id   = UID,
    gene_name = NAME,
    gweight   = GWEIGHT
  ) %>%
  mutate(
    gene_id   = stringr::str_to_upper(stringr::str_trim(as.character(gene_id))),
    gene_name = stringr::str_trim(as.character(gene_name))
  ) %>%
  select(-matches("^\\.\\.\\.")) %>%
  mutate(across(-c(gene_id, gene_name, gweight), ~ {
    x <- stringr::str_trim(as.character(.x))
    x <- na_if(x, "")
    suppressWarnings(as.numeric(x))
  }))

gasch_long <- gasch_df %>%
  pivot_longer(
    cols = -c(gene_id, gene_name, gweight),
    names_to = "condition_time_raw",
    values_to = "log2FC"
  ) %>%
  mutate(condition_time_raw = .clean_col_text(condition_time_raw))

gasch_long <- gasch_long %>%
  mutate(
    time_in_paren = stringr::str_match(
      condition_time_raw,
      "\\((\\s*[0-9]+\\.?[0-9]*)\\s*(min|mins|minutes|h|hr|hrs|hours|d|day|days)\\s*\\)"
    ),
    time_value_paren = suppressWarnings(as.numeric(time_in_paren[, 2])),
    time_unit_paren  = time_in_paren[, 3],
    
    time_trail = stringr::str_match(
      condition_time_raw,
      "\\s([0-9]+\\.?[0-9]*)\\s*(minutes|minute|mins|min\\.?|hours|hour|hrs|hr|h|days|day|d)\\.?$"
    ),
    time_value_trail = suppressWarnings(as.numeric(time_trail[, 2])),
    time_unit_trail  = time_trail[, 3],
    
    time_value = coalesce(time_value_paren, time_value_trail),
    time_unit  = coalesce(time_unit_paren,  time_unit_trail),
    time_unit  = .standardize_time_unit(time_unit),
    time_min   = .time_to_minutes(time_value, time_unit),
    
    condition_raw = case_when(
      !is.na(time_value_paren) ~ stringr::str_remove(
        condition_time_raw,
        "\\s*\\(\\s*[0-9]+\\.?[0-9]*\\s*(min|mins|minutes|h|hr|hrs|hours|d|day|days)\\s*\\)\\s*.*$"
      ),
      !is.na(time_value_trail) ~ stringr::str_remove(
        condition_time_raw,
        "\\s*[0-9]+\\.?[0-9]*\\s*(minutes|minute|mins|min\\.?|hours|hour|hrs|hr|h|days|day|d)\\.?$"
      ),
      TRUE ~ condition_time_raw
    ),
    condition = .clean_condition(condition_raw),
    
    time_label = case_when(
      time_unit == "min" ~ paste0(sprintf("%g", time_value), "m"),
      time_unit == "h"   ~ paste0(sprintf("%g", time_value), "h"),
      time_unit == "d"   ~ paste0(sprintf("%g", time_value), "d"),
      TRUE ~ NA_character_
    )
  ) %>%
  select(gene_id, gene_name, gweight, condition, time_value, time_unit, time_min, time_label, log2FC) %>%
  filter(!is.na(time_min))

gasch_conditions <- gasch_long %>%
  distinct(condition) %>%
  arrange(condition) %>%
  pull(condition)

#  Join OSR backbone 
gasch_osr_long <- gasch_long %>%
  inner_join(osr_backbone, by = "gene_id") %>%
  mutate(functional_category = factor(functional_category, levels = names(osr_category_colors)))

#  Build matrix + annotation 
make_gasch_osr_inputs <- function(condition_keep, gasch_osr_long, osr_backbone) {
  
  df <- gasch_osr_long %>% filter(condition == condition_keep)
  
  time_levels <- df %>%
    distinct(time_label, time_min) %>%
    arrange(time_min, time_label) %>%
    pull(time_label)
  
  sum_df <- df %>%
    group_by(gene_id, time_label) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    mutate(time_label = factor(time_label, levels = time_levels))
  
  backbone2 <- osr_backbone %>%
    arrange(category_order, scer_gene_name) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  joined <- backbone2 %>%
    left_join(sum_df, by = "gene_id") %>%
    mutate(time_label = factor(as.character(time_label), levels = time_levels))
  
  mat <- joined %>%
    select(scer_gene_name, time_label, log2FC) %>%
    pivot_wider(names_from = time_label, values_from = log2FC) %>%
    column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  mat <- mat[, time_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  
  rn <- rownames(mat)
  rn[rn == "" | is.na(rn)] <- "UNLABELED"
  rownames(mat) <- make.unique(rn)
  
  row_annot <- backbone2 %>%
    mutate(scer_gene_name = make.unique(scer_gene_name)) %>%
    filter(scer_gene_name %in% rownames(mat)) %>%
    arrange(match(scer_gene_name, rownames(mat))) %>%
    select(scer_gene_name, functional_category, category_order)
  
  list(mat = mat, row_annot = row_annot, time_levels = time_levels)
}

#  Plot (ONE condition)
make_gasch_osr_heatmap <- function(condition_keep,
                                   gasch_osr_long,
                                   osr_backbone,
                                   osr_category_colors,
                                   lfc_limits = c(-4, 4),
                                   na_col = "grey85",
                                   row_font = 5,
                                   col_font = 9,
                                   cell_w_mm = 5,
                                   cell_h_mm = 0.8,
                                   row_name_max_mm = 70,
                                   draw_now = TRUE) {
  
  inputs <- make_gasch_osr_inputs(condition_keep, gasch_osr_long, osr_backbone)
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
    
    show_row_names = TRUE,
    row_labels = rownames(mat),
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    left_annotation = row_ha,
    
    show_column_names = TRUE,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold"),
    
    width  = grid::unit(ncol(mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_h_mm, "mm"),
    
    column_title = stringr::str_to_title(condition_keep),
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    
    heatmap_legend_param = list(
      title = "log2FC",
      at = c(lfc_limits[1], 0, lfc_limits[2]),
      title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
    ),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat, row_annot = row_annot, inputs = inputs)
  
  if (draw_now) {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = grid::unit(c(2, 2, 2, 6), "mm")
    )
  }
  
  invisible(out)
}

# Heat-shock TEMP-SERIES combo heatmap 
make_gasch_osr_heatmap_multi_condition <- function(conditions_keep,
                                                   gasch_osr_long,
                                                   osr_backbone,
                                                   osr_category_colors,
                                                   lfc_limits = c(-4, 4),
                                                   na_col = "grey85",
                                                   row_font = 5,
                                                   col_font = 9,
                                                   cell_w_mm = 7,
                                                   cell_h_mm = 0.8,
                                                   row_name_max_mm = 70,
                                                   column_title = "Heat shock (20 min): start temp \u2192 37\u00b0C",
                                                   draw_now = TRUE) {
  
  stopifnot(length(conditions_keep) >= 2)
  
  backbone2 <- osr_backbone %>%
    arrange(category_order, scer_gene_name) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  df <- gasch_osr_long %>%
    filter(condition %in% conditions_keep) %>%
    mutate(
      gene_id = stringr::str_to_upper(stringr::str_trim(as.character(gene_id))),
      condition = as.character(condition)
    )
  
  sum_df <- df %>%
    group_by(gene_id, condition) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    mutate(log2FC = ifelse(is.nan(log2FC), NA_real_, log2FC))
  
  joined <- backbone2 %>%
    left_join(sum_df, by = "gene_id")
  
  mat <- joined %>%
    select(scer_gene_name, condition, log2FC) %>%
    tidyr::pivot_wider(names_from = condition, values_from = log2FC) %>%
    tibble::column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  missing_cols <- setdiff(conditions_keep, colnames(mat))
  if (length(missing_cols) > 0) {
    stop("These conditions were not found as columns in the matrix: ",
         paste(missing_cols, collapse = ", "))
  }
  mat <- mat[, conditions_keep, drop = FALSE]
  storage.mode(mat) <- "numeric"
  
  rn <- rownames(mat)
  rn[rn == "" | is.na(rn)] <- "UNLABELED"
  rownames(mat) <- make.unique(rn)
  
  row_annot <- backbone2 %>%
    mutate(scer_gene_name = make.unique(scer_gene_name)) %>%
    filter(scer_gene_name %in% rownames(mat)) %>%
    arrange(match(scer_gene_name, rownames(mat))) %>%
    select(scer_gene_name, functional_category, category_order)
  
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
  
  pretty_cols <- stringr::str_match(colnames(mat), "heat shock\\s+([0-9]+)\\s+to\\s+37")[, 2]
  pretty_cols <- ifelse(is.na(pretty_cols), colnames(mat), paste0(pretty_cols, "\u2192", "37"))
  
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
    column_labels = pretty_cols,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold"),
    
    width  = grid::unit(ncol(mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_h_mm, "mm"),
    
    column_title = column_title,
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    
    heatmap_legend_param = list(
      title = "log2FC",
      at = c(lfc_limits[1], 0, lfc_limits[2]),
      title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
    ),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat, row_annot = row_annot, conditions_keep = conditions_keep)
  
  if (draw_now) {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = grid::unit(c(6, 2, 2, 6), "mm")
    )
  }
  
  invisible(out)
}

# ============================================================
# OVERVIEW HEATMAP

make_gasch_osr_overview_heatmap <- function(gasch_ht_list,
                                            heat_shock_temp_series_ht,
                                            block_order,
                                            pretty_block,
                                            osr_backbone,
                                            osr_category_colors,
                                            lfc_limits = c(-4, 4),
                                            na_col = "grey85",
                                            row_font = 2.5,
                                            col_font = 7,
                                            cell_w_mm = 3.2,
                                            cell_h_mm = 0.8,
                                            row_name_max_mm = 70,
                                            title = " ",
                                            block_label_rot = 45,
                                            block_label_font = 10,
                                            column_gap_mm = 1,
                                            top_padding_mm = 10,
                                            draw_now = TRUE) {
  
  backbone2 <- osr_backbone %>%
    arrange(category_order, scer_gene_name) %>%
    mutate(scer_gene_name = make.unique(scer_gene_name)) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  row_order <- backbone2$scer_gene_name
  
  block_mats <- list()
  block_labels_for_cols <- character(0)
  col_labels_for_cols <- character(0)
  
  for (bk in block_order) {
    
    if (bk == "heat shock temp series") {
      mat <- heat_shock_temp_series_ht$mat
      mat <- mat[row_order, , drop = FALSE]
      block_mats[[bk]] <- mat
      
      block_labels_for_cols <- c(block_labels_for_cols, rep(bk, ncol(mat)))
      col_labels_for_cols   <- c(col_labels_for_cols, colnames(mat))
      
    } else {
      if (!bk %in% names(gasch_ht_list)) stop("Missing block in gasch_ht_list: ", bk)
      mat <- gasch_ht_list[[bk]]$mat
      mat <- mat[row_order, , drop = FALSE]
      block_mats[[bk]] <- mat
      
      block_labels_for_cols <- c(block_labels_for_cols, rep(bk, ncol(mat)))
      col_labels_for_cols   <- c(col_labels_for_cols, colnames(mat))
    }
  }
  
  big_mat <- do.call(cbind, block_mats)
  storage.mode(big_mat) <- "numeric"
  
  row_annot <- backbone2 %>%
    select(scer_gene_name, functional_category, category_order) %>%
    arrange(match(scer_gene_name, rownames(big_mat)))
  
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
  
  # split columns by condition block
  block_factor <- factor(block_labels_for_cols, levels = block_order)
  block_levels <- levels(block_factor)
  
  # pretty labels for each block (one per slice)
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
    
    # block grouping
    column_split = block_factor,
    column_gap   = grid::unit(column_gap_mm, "mm"),
    
    # block labels ABOVE (no colored strip)
    column_title = slice_titles,
    column_title_side = "top",
    column_title_gp = grid::gpar(fontsize = block_label_font, fontface = "bold"),
    column_title_rot = block_label_rot,
    
    # gene names
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    row_names_max_width = grid::unit(row_name_max_mm, "mm"),
    
    # time labels
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
  
  out <- list(
    ht = ht,
    mat = big_mat,
    block_order = block_order,
    block_pretty = block_pretty
  )
  
  if (draw_now) {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = grid::unit(c(top_padding_mm, 2, 2, 6), "mm")
    )
    # overall title (top of device), separate from per-slice titles
    grid::grid.text(
      title,
      x = grid::unit(0.5, "npc"),
      y = grid::unit(0.995, "npc"),
      just = "top",
      gp = grid::gpar(fontsize = 14, fontface = "bold")
    )
  }
  
  invisible(out)
}

# ============================================================
# Build callable objects

# (A) Single example
h2o2_ht <- make_gasch_osr_heatmap(
  condition_keep = "constant 0.32 mm h2o2",
  gasch_osr_long = gasch_osr_long,
  osr_backbone = osr_backbone,
  osr_category_colors = osr_category_colors,
  draw_now = TRUE
)

# (B) Per-condition heatmap objects
gasch_conditions_keep <- setdiff(gasch_conditions, "1mm menadione")

gasch_ht_list <- setNames(vector("list", length(gasch_conditions_keep)), gasch_conditions_keep)
for (cond in gasch_conditions_keep) {
  gasch_ht_list[[cond]] <- make_gasch_osr_heatmap(
    condition_keep = cond,
    gasch_osr_long = gasch_osr_long,
    osr_backbone = osr_backbone,
    osr_category_colors = osr_category_colors,
    draw_now = FALSE
  )
}

# (C) Heat-shock temp series combo
heatshock_temp_series_conditions <- c(
  "heat shock 17 to 37,",
  "heat shock 21 to 37,",
  "heat shock 25 to 37,",
  "heat shock 29 to 37,",
  "heat shock 33 to 37,",
  "heat shock 37 to 37,"
)

heat_shock_temp_series_ht <- make_gasch_osr_heatmap_multi_condition(
  conditions_keep = heatshock_temp_series_conditions,
  gasch_osr_long = gasch_osr_long,
  osr_backbone = osr_backbone,
  osr_category_colors = osr_category_colors,
  draw_now = TRUE
)

# single-condition heatmaps
ypd <- ComplexHeatmap::draw(gasch_ht_list[["ypd"]]$ht)
menadione_1mm <- ComplexHeatmap::draw(gasch_ht_list[["1 mm menadione"]]$ht)
diamide_1p5mm <- ComplexHeatmap::draw(gasch_ht_list[["1.5 mm diamide"]]$ht)
sorbitol_1m <- ComplexHeatmap::draw(gasch_ht_list[["1m sorbitol -"]]$ht)
aa_starv <- ComplexHeatmap::draw(gasch_ht_list[["aa starv"]]$ht)
dtt <- ComplexHeatmap::draw(gasch_ht_list[["dtt"]]$ht)
heat_shock <- ComplexHeatmap::draw(gasch_ht_list[["heat shock"]]$ht)
nitrogen <- ComplexHeatmap::draw(gasch_ht_list[["nitrogen depletion"]]$ht)

# ============================================================
# OVERVIEW HEATMAP: combine all selected conditions into ONE
block_order <- c(
  "heat shock",
  "heat shock temp series",
  "constant 0.32 mm h2o2",
  "1 mm menadione",
  "dtt",
  "1.5 mm diamide",
  "1m sorbitol -",
  "aa starv",
  "nitrogen depletion",
  "ypd"
)

pretty_block <- c(
  "heat shock"               = "37\u00b0C heat shock",
  "heat shock temp series"   = "Variable temperature shocks",
  "constant 0.32 mm h2o2"    = "Hydrogen peroxide",
  "1 mm menadione"           = "Menadione",
  "dtt"                      = "DTT",
  "1.5 mm diamide"           = "Diamide",
  "1m sorbitol -"            = "Sorbitol osmotic shock",
  "aa starv"                 = "Amino acid starvation",
  "nitrogen depletion"       = "Nitrogen depletion",
  "ypd"                      = "Stationary phase (YPD)"
)

plots_dir <- here("08.plots", "Scer")

fig_legend <- paste(
  " ",
  sep = ""
)

overview_ht <- make_gasch_osr_overview_heatmap(
  gasch_ht_list = gasch_ht_list,
  heat_shock_temp_series_ht = heat_shock_temp_series_ht,
  block_order = block_order,
  pretty_block = pretty_block,
  osr_backbone = osr_backbone,
  osr_category_colors = osr_category_colors,
  block_label_rot = 45,
  block_label_font = 10,
  top_padding_mm = 7,
  draw_now = FALSE
)
draw_overview_with_legend <- function(ht_obj, legend_text,
                                      legend_font = 10,
                                      legend_height_mm = 32,
                                      legend_side_margin_mm = 10) {
  
  ht_grob <- grid::grid.grabExpr(
    ComplexHeatmap::draw(
      ht_obj$ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = grid::unit(c(2, 2, 2, 6), "mm")
    )
  )
  
  legend_grob <- grid::grobTree(
    vp = grid::viewport(
      x = grid::unit(legend_side_margin_mm, "mm"),
      y = grid::unit(1, "npc"),
      just = c("left", "top"),
      width  = grid::unit(1, "npc") - grid::unit(2 * legend_side_margin_mm, "mm"),
      height = grid::unit(1, "npc")
    ),
    grid::textGrob(
      legend_text,
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = legend_font, lineheight = 1.05)
    )
  )
  
  grid::grid.draw(
    gridExtra::arrangeGrob(
      ht_grob,
      legend_grob,
      ncol = 1,
      heights = grid::unit.c(
        grid::unit(1, "null"),
        grid::unit(legend_height_mm, "mm")
      )
    )
  )
}
out_pdf <- file.path(plots_dir, "Gasch2000_OSR_overview_heatmap.pdf")
pdf(out_pdf, width = 12, height = 6)
draw_overview_with_legend(overview_ht, fig_legend)
dev.off()
