---
title: "Gasch et al. 2000 S. cer ESR genes  "
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
# Inputs:
# 1) Whole Gasch dataset: 20260109_Gasch2000_ProcessedLog2Data.xlsx (sheet2)
# 2) ESR gene set:        20250109_ESR_clusters_UPDATED_2017.xlsx (sheet2)
#
# Steps:
# - Read Gasch full dataset (UID/NAME/GWEIGHT + conditions)
# - Read ESR set (UID/NAME) and keep NAME so we can parse functional_category + order
# - Subset Gasch to ESR genes by UID
# - Use OSR pipeline’s condition/time parsing to generate tidy long table
# - Build per-condition ESR heatmaps (Gasch conditions)

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
  library(here)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(gridExtra)
  
})


`%||%` <- function(x, y) if (!is.null(x)) x else y

# ============================================================
# Paths

scer_dir  <- here("06.results", "Scer")
plots_dir <- here("08.plots", "Scer")

gasch_path <- file.path(scer_dir, "20260109_Gasch2000_ProcessedLog2Data.xlsx")
if (!file.exists(gasch_path)) stop("Gasch file not found: ", gasch_path)

esr_path  <- file.path(scer_dir, "20250109_ESR_clusters_UPDATED_2017.xlsx")
esr_sheet <- 2
if (!file.exists(esr_path)) stop("ESR file not found: ", esr_path)

if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# ============================================================
# Display params

lfc_limits <- c(-4, 4)
na_col     <- "grey85"

single_cell_w_mm <- 5
single_cell_h_mm <- 0.20
single_col_font  <- 8

overview_cell_w_mm <- 3.0
overview_cell_h_mm <- 0.20
overview_col_font  <- 7
overview_block_rot <- 45

# ============================================================
# Helpers: Gasch header cleanup + condition/time parsing

.clean_col_text <- function(x) {
  x %>%
    str_replace_all("\\s+", " ") %>%
    str_trim() %>%
    str_replace_all("030inutes", "030 minutes") %>%
    str_replace_all("\\(\\s*", "(") %>%
    str_replace_all("\\s*\\)", ")")
}

.standardize_time_unit <- function(u) {
  u <- str_to_lower(str_trim(u))
  case_when(
    u %in% c("min", "mins", "minute", "minutes", "min.") ~ "min",
    u %in% c("h", "hr", "hrs", "hour", "hours") ~ "h",
    u %in% c("d", "day", "days") ~ "d",
    TRUE ~ NA_character_
  )
}

.time_to_minutes <- function(v, u) {
  u <- .standardize_time_unit(u)
  case_when(
    is.na(v) | is.na(u) ~ NA_real_,
    u == "min" ~ as.numeric(v),
    u == "h"   ~ as.numeric(v) * 60,
    u == "d"   ~ as.numeric(v) * 1440,
    TRUE ~ NA_real_
  )
}

.clean_condition <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("redo|rescan", "") %>%
    str_replace_all("\\*.*?\\)", ")") %>%
    str_replace_all("\\(.*?problem.*?\\)", "") %>%
    str_replace_all("reference pool", "") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

# ============================================================
# STANDARDIZE CONDITION NAMES (CRITICAL)
# Ensures gasch_ht_list has EXACT keys used downstream.

.standardize_gasch_condition <- function(cond) {
  c0 <- str_to_lower(str_trim(cond))
  
  # normalize menadione variants
  c0 <- case_when(
    str_detect(c0, "\\bmenadione\\b") ~ "1 mm menadione",
    TRUE ~ c0
  )
  
  # normalize peroxide variants (keep your canonical label)
  c0 <- case_when(
    str_detect(c0, "h2o2") ~ "constant 0.32 mm h2o2",
    TRUE ~ c0
  )
  
  # normalize ypd variants
  c0 <- case_when(
    str_detect(c0, "^ypd\\b") ~ "ypd",
    TRUE ~ c0
  )
  
  # leave others as-is
  c0
}

# ============================================================
# ESR parsing + functional category collapse

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
# Read ESR (UID + NAME)

esr_set_raw <- suppressMessages(
  read_excel(esr_path, sheet = esr_sheet, col_types = "text", .name_repair = "minimal")
)
names(esr_set_raw) <- make.names(names(esr_set_raw), unique = TRUE)

if (!all(c("UID", "NAME") %in% names(esr_set_raw))) {
  stop("Expected columns UID and NAME in ESR file sheet ", esr_sheet,
       ". Found: ", paste(names(esr_set_raw), collapse = ", "))
}

esr_set_raw <- esr_set_raw %>% select(UID, NAME)
esr_name_parsed <- parse_esr_NAME(esr_set_raw$NAME)

esr_backbone <- esr_set_raw %>%
  transmute(
    gene_id  = str_to_upper(str_trim(as.character(UID))),
    name_raw = str_trim(as.character(NAME))
  ) %>%
  bind_cols(esr_name_parsed %>% select(gene_name, desc_raw)) %>%
  mutate(
    scer_gene_name = gene_name,
    functional_category = collapse_esr_function(desc_raw)
  ) %>%
  left_join(esr_category_order_tbl, by = "functional_category") %>%
  mutate(
    scer_gene_name = if_else(is.na(scer_gene_name) | scer_gene_name == "", gene_id, scer_gene_name),
    scer_gene_name = str_to_upper(str_trim(scer_gene_name)),
    functional_category = factor(functional_category, levels = esr_category_order_tbl$functional_category),
    category_order = if_else(is.na(category_order), 99, category_order)
  ) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(category_order, scer_gene_name)

write.csv(
  esr_backbone %>% select(gene_id, scer_gene_name, functional_category, category_order, name_raw, desc_raw),
  file = file.path(scer_dir, "Gasch2000_ESR_backbone_parsed.csv"),
  row.names = FALSE
)

# ============================================================
# Read Gasch sheet2 with SAFE names but PARSE from RAW headers

gasch_raw <- suppressMessages(
  read_excel(gasch_path, sheet = 2, col_types = "text", .name_repair = "minimal")
)

gasch_colnames_raw <- names(gasch_raw)
names(gasch_raw) <- make.names(names(gasch_raw), unique = TRUE)

colmap <- tibble(
  col_safe = names(gasch_raw),
  col_raw  = gasch_colnames_raw
) %>%
  mutate(
    col_raw = if_else(is.na(col_raw) | str_trim(col_raw) == "", col_safe, col_raw),
    col_raw = as.character(col_raw)
  )

if (!all(c("UID", "NAME", "GWEIGHT") %in% names(gasch_raw))) {
  stop("Expected UID/NAME/GWEIGHT in Gasch sheet2 AFTER make.names(). Found: ",
       paste(names(gasch_raw), collapse = ", "))
}

gasch_df <- gasch_raw %>%
  rename(
    gene_id   = UID,
    gene_name = NAME,
    gweight   = GWEIGHT
  ) %>%
  mutate(
    gene_id   = str_to_upper(str_trim(as.character(gene_id))),
    gene_name = str_trim(as.character(gene_name))
  ) %>%
  mutate(across(-c(gene_id, gene_name, gweight), ~ {
    x <- str_trim(as.character(.x))
    x <- na_if(x, "")
    suppressWarnings(as.numeric(x))
  }))

gasch_long <- gasch_df %>%
  pivot_longer(
    cols = -c(gene_id, gene_name, gweight),
    names_to = "col_safe",
    values_to = "log2FC"
  ) %>%
  left_join(colmap, by = "col_safe") %>%
  mutate(condition_time_raw = .clean_col_text(col_raw)) %>%
  select(-col_safe, -col_raw)

gasch_long <- gasch_long %>%
  mutate(
    time_in_paren = str_match(
      condition_time_raw,
      "\\((\\s*[0-9]+\\.?[0-9]*)\\s*(min|mins|minutes|h|hr|hrs|hours|d|day|days)\\s*\\)"
    ),
    time_value_paren = suppressWarnings(as.numeric(time_in_paren[, 2])),
    time_unit_paren  = time_in_paren[, 3],
    
    time_trail = str_match(
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
      !is.na(time_value_paren) ~ str_remove(
        condition_time_raw,
        "\\s*\\(\\s*[0-9]+\\.?[0-9]*\\s*(min|mins|minutes|h|hr|hrs|hours|d|day|days)\\s*\\)\\s*.*$"
      ),
      !is.na(time_value_trail) ~ str_remove(
        condition_time_raw,
        "\\s*[0-9]+\\.?[0-9]*\\s*(minutes|minute|mins|min\\.?|hours|hour|hrs|hr|h|days|day|d)\\.?$"
      ),
      TRUE ~ condition_time_raw
    ),
    condition = .clean_condition(condition_raw),
    # APPLY STANDARDIZATION HERE
    condition = .standardize_gasch_condition(condition),
    
    time_label = case_when(
      time_unit == "min" ~ paste0(sprintf("%g", time_value), "m"),
      time_unit == "h"   ~ paste0(sprintf("%g", time_value), "h"),
      time_unit == "d"   ~ paste0(sprintf("%g", time_value), "d"),
      TRUE ~ NA_character_
    )
  ) %>%
  select(gene_id, gene_name, gweight, condition, time_value, time_unit, time_min, time_label, log2FC) %>%
  filter(!is.na(time_min), !is.na(time_label), time_label != "")

gasch_conditions <- gasch_long %>%
  distinct(condition) %>%
  arrange(condition) %>%
  pull(condition)

message("Parsed conditions (standardized) (n=", length(gasch_conditions), "):")
print(gasch_conditions)

write.csv(
  tibble(condition = gasch_conditions),
  file = file.path(scer_dir, "Gasch2000_parsed_conditions_STANDARDIZED.csv"),
  row.names = FALSE
)

# ============================================================
# Subset ESR genes from Gasch (join on UID)

gasch_esr_long <- gasch_long %>%
  inner_join(
    esr_backbone %>% select(gene_id, scer_gene_name, functional_category, category_order, name_raw),
    by = "gene_id"
  ) %>%
  mutate(functional_category = factor(functional_category, levels = esr_category_order_tbl$functional_category))

# ============================================================
# MASTER backbone row order used everywhere

esr_backbone_master <- esr_backbone %>%
  arrange(category_order, scer_gene_name) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(scer_gene_name = make.unique(scer_gene_name))

master_row_order <- esr_backbone_master$scer_gene_name

# ============================================================
# Build matrix + annotation per condition (duplicate-safe)

make_gasch_esr_inputs <- function(condition_keep, gasch_esr_long, backbone_master) {
  
  df <- gasch_esr_long %>% filter(condition == condition_keep)
  if (nrow(df) == 0) stop("No rows for condition: ", condition_keep)
  
  time_levels <- df %>%
    distinct(time_label, time_min) %>%
    arrange(time_min, time_label) %>%
    pull(time_label)
  
  sum_df <- df %>%
    select(gene_id, time_label, log2FC) %>%
    mutate(time_label = factor(time_label, levels = time_levels)) %>%
    group_by(gene_id, time_label) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop")
  
  joined <- backbone_master %>% left_join(sum_df, by = "gene_id")
  
  mat <- joined %>%
    select(scer_gene_name, time_label, log2FC) %>%
    mutate(time_label = as.character(time_label)) %>%
    filter(!is.na(time_label)) %>%
    pivot_wider(
      names_from  = time_label,
      values_from = log2FC,
      values_fn   = mean
    ) %>%
    column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  missing_cols <- setdiff(time_levels, colnames(mat))
  if (length(missing_cols) > 0) {
    for (cc in missing_cols) mat <- cbind(mat, setNames(rep(NA_real_, nrow(mat)), cc))
  }
  mat <- mat[, time_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  
  # enforce master row order SAFELY
  common_rows <- intersect(master_row_order, rownames(mat))
  mat <- mat[common_rows, , drop = FALSE]
  
  row_annot <- backbone_master %>%
    filter(scer_gene_name %in% common_rows) %>%
    arrange(match(scer_gene_name, common_rows)) %>%
    select(scer_gene_name, functional_category, category_order)
  
  list(mat = mat, row_annot = row_annot, time_levels = time_levels)
}

# ============================================================
# Plot per-condition ESR heatmap

make_gasch_esr_heatmap <- function(condition_keep,
                                   gasch_esr_long,
                                   backbone_master,
                                   esr_category_colors,
                                   lfc_limits = c(-4, 4),
                                   na_col = "grey85",
                                   col_font = 8,
                                   cell_w_mm = 5,
                                   cell_h_mm = 0.20,
                                   draw_now = TRUE) {
  
  inputs <- make_gasch_esr_inputs(condition_keep, gasch_esr_long, backbone_master)
  mat <- inputs$mat
  row_annot <- inputs$row_annot
  
  row_ha <- rowAnnotation(
    Function = row_annot$functional_category,
    col = list(Function = esr_category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
  )
  
  ht <- Heatmap(
    mat,
    name = "log2FC",
    na_col = na_col,
    col = colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    left_annotation = row_ha,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = col_font, fontface = "bold"),
    width  = unit(ncol(mat) * cell_w_mm, "mm"),
    height = unit(nrow(mat) * cell_h_mm, "mm"),
    column_title = str_to_title(condition_keep),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat)
  if (draw_now) draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  invisible(out)
}

# ============================================================
# Build per-condition heatmap objects (ALL standardized conditions)

gasch_ht_list <- setNames(vector("list", length(gasch_conditions)), gasch_conditions)

for (cond in gasch_conditions) {
  gasch_ht_list[[cond]] <- make_gasch_esr_heatmap(
    condition_keep   = cond,
    gasch_esr_long   = gasch_esr_long,
    backbone_master  = esr_backbone_master,
    esr_category_colors = esr_category_colors,
    lfc_limits = lfc_limits,
    na_col = na_col,
    col_font = single_col_font,
    cell_w_mm = single_cell_w_mm,
    cell_h_mm = single_cell_h_mm,
    draw_now = FALSE
  )
}

message("Built heatmaps for ", length(gasch_ht_list), " standardized conditions.")

# ============================================================
# (C) Variable-temperature heat-shock series block (multi-condition)
# IMPORTANT: we map these raw names into a SINGLE synthetic block.

heatshock_temp_series_conditions_raw <- c(
  "heat shock 17 to 37,",
  "heat shock 21 to 37,",
  "heat shock 25 to 37,",
  "heat shock 29 to 37,",
  "heat shock 33 to 37,",
  "heat shock 37 to 37,"
)

# these should still exist (they are not standardized away)
missing_hs <- setdiff(heatshock_temp_series_conditions_raw, unique(gasch_long$condition_raw %||% character(0)))
# we didn't keep condition_raw, so just check against pre-standardized set:
# safest is check against the *pre-standard* list from file written earlier if needed.
# For now, check in the (non-standardized) printed list you saw; if you changed the excel, update.

make_gasch_esr_heatmap_multi_condition <- function(conditions_keep,
                                                   gasch_esr_long,
                                                   backbone_master,
                                                   esr_category_colors,
                                                   lfc_limits = c(-4, 4),
                                                   na_col = "grey85",
                                                   col_font = 7,
                                                   cell_w_mm = 3.0,
                                                   cell_h_mm = 0.20,
                                                   title = "Variable temperature shocks",
                                                   draw_now = TRUE) {
  
  df <- gasch_esr_long %>% filter(condition %in% conditions_keep)
  if (nrow(df) == 0) {
    stop("No rows found for multi-condition block. Conditions requested:\n- ",
         paste(conditions_keep, collapse = "\n- "))
  }
  
  df <- df %>% mutate(condition = factor(condition, levels = conditions_keep))
  
  time_tbl <- df %>%
    distinct(condition, time_label, time_min) %>%
    arrange(condition, time_min, time_label)
  
  sum_df <- df %>%
    select(gene_id, condition, time_label, log2FC) %>%
    group_by(gene_id, condition, time_label) %>%
    summarise(log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
    mutate(col_key = paste0(as.character(condition), " | ", as.character(time_label)))
  
  col_levels <- time_tbl %>%
    mutate(col_key = paste0(as.character(condition), " | ", as.character(time_label))) %>%
    pull(col_key) %>%
    unique()
  
  joined <- backbone_master %>% left_join(sum_df, by = "gene_id")
  
  mat <- joined %>%
    select(scer_gene_name, col_key, log2FC) %>%
    filter(!is.na(col_key)) %>%
    pivot_wider(
      names_from = col_key,
      values_from = log2FC,
      values_fn = mean
    ) %>%
    column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  missing_cols <- setdiff(col_levels, colnames(mat))
  if (length(missing_cols) > 0) {
    for (cc in missing_cols) mat <- cbind(mat, setNames(rep(NA_real_, nrow(mat)), cc))
  }
  mat <- mat[, col_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  
  common_rows <- intersect(master_row_order, rownames(mat))
  mat <- mat[common_rows, , drop = FALSE]
  
  row_annot <- backbone_master %>%
    filter(scer_gene_name %in% common_rows) %>%
    arrange(match(scer_gene_name, common_rows))
  
  row_ha <- rowAnnotation(
    Function = row_annot$functional_category,
    col = list(Function = esr_category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
  )
  
  ht <- Heatmap(
    mat,
    name = "log2FC",
    na_col = na_col,
    col = colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    left_annotation = row_ha,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = col_font, fontface = "bold", rot = 90),
    width  = unit(ncol(mat) * cell_w_mm, "mm"),
    height = unit(nrow(mat) * cell_h_mm, "mm"),
    column_title = title,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat)
  if (draw_now) draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  invisible(out)
}

# NOTE: keep these as their raw condition strings (they are NOT standardized)
heat_shock_temp_series_ht <- make_gasch_esr_heatmap_multi_condition(
  conditions_keep = heatshock_temp_series_conditions_raw,
  gasch_esr_long = gasch_esr_long,
  backbone_master = esr_backbone_master,
  esr_category_colors = esr_category_colors,
  lfc_limits = lfc_limits,
  na_col = na_col,
  draw_now = FALSE,
  title = "Variable temperature shocks"
)

gasch_ht_list[["heat shock temp series"]] <- heat_shock_temp_series_ht

# ============================================================
# Overview heatmap builder (SAFE row indexing + pads missing rows)

make_gasch_esr_overview_heatmap <- function(gasch_ht_list,
                                            block_order,
                                            pretty_block,
                                            backbone_master,
                                            esr_category_colors,
                                            lfc_limits = c(-4, 4),
                                            na_col = "grey85",
                                            col_font = 7,
                                            cell_w_mm = 3.0,
                                            cell_h_mm = 0.20,
                                            block_label_rot = 45,
                                            block_label_font = 10,
                                            column_gap_mm = 1,
                                            top_padding_mm = 10,
                                            title = "",
                                            draw_now = TRUE) {
  
  backbone2 <- backbone_master %>%
    mutate(scer_gene_name = make.unique(scer_gene_name)) %>%
    distinct(gene_id, .keep_all = TRUE)
  
  row_order <- backbone2$scer_gene_name
  
  block_mats <- list()
  block_labels_for_cols <- character(0)
  col_labels_for_cols <- character(0)
  
  for (bk in block_order) {
    if (!bk %in% names(gasch_ht_list)) stop("Missing block in gasch_ht_list: ", bk)
    mat <- gasch_ht_list[[bk]]$mat
    
    # align and PAD to full row_order
    common_rows <- intersect(row_order, rownames(mat))
    mat2 <- mat[common_rows, , drop = FALSE]
    
    missing_rows <- setdiff(row_order, rownames(mat2))
    if (length(missing_rows) > 0) {
      pad <- matrix(NA_real_, nrow = length(missing_rows), ncol = ncol(mat2),
                    dimnames = list(missing_rows, colnames(mat2)))
      mat2 <- rbind(mat2, pad)
    }
    mat2 <- mat2[row_order, , drop = FALSE]
    
    block_mats[[bk]] <- mat2
    block_labels_for_cols <- c(block_labels_for_cols, rep(bk, ncol(mat2)))
    col_labels_for_cols   <- c(col_labels_for_cols, colnames(mat2))
  }
  
  big_mat <- do.call(cbind, block_mats)
  storage.mode(big_mat) <- "numeric"
  
  row_annot <- backbone2 %>%
    select(scer_gene_name, functional_category, category_order) %>%
    arrange(match(scer_gene_name, rownames(big_mat)))
  
  row_ha <- rowAnnotation(
    Function = row_annot$functional_category,
    col = list(Function = esr_category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
  )
  
  block_factor <- factor(block_labels_for_cols, levels = block_order)
  
  block_levels <- levels(block_factor)
  block_pretty <- pretty_block[block_levels]
  block_pretty[is.na(block_pretty) | block_pretty == ""] <- block_levels
  slice_titles <- unname(block_pretty)
  
  ht <- Heatmap(
    big_mat,
    name = "log2FC",
    na_col = na_col,
    col = colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    left_annotation = row_ha,
    column_split = block_factor,
    column_gap   = unit(column_gap_mm, "mm"),
    column_title = slice_titles,
    column_title_side = "top",
    column_title_gp = gpar(fontsize = block_label_font, fontface = "bold"),
    column_title_rot = block_label_rot,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_labels = col_labels_for_cols,
    column_names_gp = gpar(fontsize = col_font, fontface = "bold", rot = 90),
    width  = unit(ncol(big_mat) * cell_w_mm, "mm"),
    height = unit(nrow(big_mat) * cell_h_mm, "mm"),
    border = FALSE
  )
  
  out <- list(ht = ht, mat = big_mat)
  
  if (draw_now) {
    draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = unit(c(top_padding_mm, 2, 2, 6), "mm")
    )
    if (!is.null(title) && title != "") {
      grid.text(title, x = unit(0.5, "npc"), y = unit(0.995, "npc"),
                just = "top", gp = gpar(fontsize = 14, fontface = "bold"))
    }
  }
  
  invisible(out)
}

# ============================================================
# SINGLE-CONDITION HEATMAPS (your standard names)
# IMPORTANT: because we standardized, keys below will exist.

ypd <- ComplexHeatmap::draw(gasch_ht_list[["ypd"]]$ht)
h2o2 <- ComplexHeatmap::draw(gasch_ht_list[["constant 0.32 mm h2o2"]]$ht)
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

overview <- make_gasch_esr_overview_heatmap(
  gasch_ht_list = gasch_ht_list,
  block_order = block_order,
  pretty_block = pretty_block,
  backbone_master = esr_backbone_master,
  esr_category_colors = esr_category_colors,
  lfc_limits = lfc_limits,
  na_col = na_col,
  col_font = overview_col_font,
  cell_w_mm = overview_cell_w_mm,
  cell_h_mm = overview_cell_h_mm,
  block_label_rot = overview_block_rot,
  column_gap_mm = 1,
  top_padding_mm = 12,
  title = " ",
  draw_now = TRUE
)

pdf(file.path(plots_dir, "Gasch2000_ESR_overview_heatmap.pdf"),
    width = 16, height = 12, useDingbats = FALSE)
ComplexHeatmap::draw(overview$ht, heatmap_legend_side = "right", annotation_legend_side = "right",
                     merge_legends = TRUE, padding = unit(c(12,2,2,6), "mm"))
grid::grid.text(" ",
                x = unit(0.5, "npc"), y = unit(0.995, "npc"),
                just = "top", gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
