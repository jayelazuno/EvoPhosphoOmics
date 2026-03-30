---
title: "Characterizing the evolutionary divergence of PHO regulon structure and 
dynamics across four yeast species"
author: "Joshua Ayelazuno"
date: "2025-11-11"
output: html_document
---
# -------------------- Packages ------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("patchwork", quietly = TRUE)) BiocManager::install("patchwork")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(readxl)
  library(here)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)      # gpar()
  library(patchwork)
})

# ============================================================
# PHO regulon cross-species time-course heatmaps (overview)
# Goal: visualize conserved vs species-specific PHO regulon dynamics across
# 4 budding yeasts by projecting each species’ time-course log2FC onto the
# same S. cerevisiae PHO regulon backbone, then concatenating heatmaps into
# a single overview with species labels as the top (slice) titles.
#
# Design/implementation:
# - Use curated S. cerevisiae PHO regulon (gene_id + display gene_name + function).
# - Map orthologs via OrthoFinder orthogroups, then apply species-specific ID maps
#   (Cgla QNG->GWK; Klac XP->KLLA).
# - Standardize each species’ time axis to ordered labels (T01..T10, species-specific),
#   build complete matrices (all PHO genes x all timepoints), and mark missing orthologs
#   as NA (drawn in grey).
# - Build 4 species heatmaps and ONE combined overview heatmap where the top labels
#   are species names; spacing between species blocks is controlled by column_gap_mm.
# ============================================================

`%||%` <- function(x, y) if (!is.null(x)) x else y

.tp_sort_key <- function(tp) {
  tp_chr <- as.character(tp)
  num <- suppressWarnings(as.numeric(stringr::str_extract(tp_chr, "\\d+\\.?\\d*")))
  tibble(tp_chr = tp_chr, num = num) %>%
    mutate(ord_num = ifelse(is.na(num), Inf, num)) %>%
    arrange(ord_num, tp_chr) %>%
    pull(tp_chr)
}

make_tp_map_10 <- function(lfc_tidy_df, species_code) {
  stopifnot(all(c("gene_id", "tp", "log2FC") %in% names(lfc_tidy_df)))
  
  tp_chr <- unique(as.character(lfc_tidy_df$tp))
  
  # known missing timepoints for Calb/Klac
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
  time_levels <- tp_map$time_std %>% as.character()
  
  lfc_tidy_df %>%
    mutate(
      gene_id = stringr::str_trim(as.character(gene_id)),
      tp      = stringr::str_trim(as.character(tp))
    ) %>%
    left_join(tp_map %>% mutate(tp = stringr::str_trim(as.character(tp))), by = "tp") %>%
    filter(!is.na(time_std)) %>%
    mutate(time_std = factor(as.character(time_std), levels = time_levels))
}

standardize_function <- function(x) {
  x <- stringr::str_trim(as.character(x))
  x <- dplyr::if_else(is.na(x) | x == "", "Other", x)
  stringr::str_to_title(x)
}

make_na_legend <- function(label = "No ortholog (NA)", col = "grey85") {
  ComplexHeatmap::Legend(
    labels = label,
    legend_gp = grid::gpar(fill = col, col = NA),
    title = NULL
  )
}

# Complete matrix over backbone display names; NA gene_id rows remain NA (grey) if requested.
make_complete_hm_matrix <- function(gene_df, lfc_std,
                                    time_levels,
                                    keep_na_if_gene_id_missing = FALSE) {
  
  gene_df <- gene_df %>%
    mutate(
      gene_id        = dplyr::na_if(stringr::str_trim(as.character(gene_id)), ""),
      scer_gene_name = stringr::str_trim(as.character(scer_gene_name)),
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
      gene_id  = dplyr::na_if(stringr::str_trim(as.character(gene_id)), ""),
      time_std = factor(as.character(time_std), levels = time_levels),
      log2FC   = suppressWarnings(as.numeric(log2FC))
    ) %>%
    select(gene_id, time_std, log2FC)
  
  id_to_name <- gene_df %>% select(gene_id, scer_gene_name) %>% distinct()
  
  complete_long <- complete_grid %>%
    left_join(id_to_name, by = "scer_gene_name") %>%
    left_join(observed, by = c("gene_id", "time_std"))
  
  if (keep_na_if_gene_id_missing) {
    # no ortholog => keep NA for all times (draw as na_col)
    complete_long <- complete_long %>%
      mutate(log2FC = if_else(is.na(gene_id), NA_real_, tidyr::replace_na(log2FC, 0)))
  } else {
    complete_long <- complete_long %>%
      mutate(log2FC = tidyr::replace_na(log2FC, 0))
  }
  
  mat <- complete_long %>%
    select(scer_gene_name, time_std, log2FC) %>%
    pivot_wider(names_from = time_std, values_from = log2FC) %>%
    tibble::column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  storage.mode(mat) <- "numeric"
  mat <- mat[, time_levels, drop = FALSE]
  
  row_annot <- gene_df %>%
    filter(scer_gene_name %in% rownames(mat)) %>%
    arrange(match(scer_gene_name, rownames(mat))) %>%
    select(scer_gene_name, functional_category, category_order)
  
  list(mat = mat, row_annot = row_annot)
}

plot_pho_heatmap <- function(mat, row_annot, category_colors,
                             lfc_limits = c(-6, 6),
                             heatmap_name = "log2FC",
                             na_col = "grey85",
                             row_font = 10,
                             col_font = 9,
                             cell_w_mm = 6,
                             cell_h_mm = 4) {
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = row_annot$functional_category,
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
  
  ComplexHeatmap::Heatmap(
    mat,
    name = heatmap_name,
    na_col = na_col,
    col = circlize::colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_side = "left",
    column_names_side = "bottom",
    row_names_gp = grid::gpar(fontsize = row_font, fontface = "bold"),
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold", rot = 90),
    left_annotation = row_ha,
    width  = grid::unit(ncol(mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_h_mm, "mm"),
    heatmap_legend_param = list(
      title = heatmap_name,
      at = c(lfc_limits[1], 0, lfc_limits[2]),
      title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 9,  fontface = "bold")
    ),
    border = FALSE,
    column_title = NULL,
    row_title = NULL
  )
}

make_pho_species_overview_heatmap <- function(mat_list,
                                              row_annot,
                                              category_colors,
                                              block_order,
                                              pretty_block,
                                              lfc_limits = c(-6, 6),
                                              na_col = "grey85",
                                              row_font = 10,
                                              col_font = 8,
                                              cell_w_mm = 3.5,
                                              cell_h_mm = 3.8,
                                              column_gap_mm = 1,
                                              block_label_rot = 0,
                                              block_label_font = 14,
                                              row_name_max_mm = 70,
                                              draw_now = TRUE) {
  
  row_order <- rownames(mat_list[[block_order[1]]])
  
  block_mats <- list()
  block_labels_for_cols <- character(0)
  col_labels_for_cols <- character(0)
  
  for (bk in block_order) {
    mat <- mat_list[[bk]][row_order, , drop = FALSE]
    block_mats[[bk]] <- mat
    block_labels_for_cols <- c(block_labels_for_cols, rep(bk, ncol(mat)))
    col_labels_for_cols   <- c(col_labels_for_cols, colnames(mat))
  }
  
  big_mat <- do.call(cbind, block_mats)
  storage.mode(big_mat) <- "numeric"
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = row_annot$functional_category,
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
  
  block_factor <- factor(block_labels_for_cols, levels = block_order)
  block_pretty <- pretty_block[levels(block_factor)]
  block_pretty[is.na(block_pretty) | block_pretty == ""] <- levels(block_factor)
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
    column_title_gp = grid::gpar(fontsize = block_label_font, fontface = "bold.italic"),
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
  
  if (draw_now) {
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = grid::unit(c(10, 2, 2, 6), "mm")
    )
  }
  
  invisible(out)
}

# -------------------- Paths / Inputs ------------------------
plots_dir <- here("08.plots", "Scer")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

category_colors <- c(
  "Transporters"      = "#E41A1C",
  "Poly P Metabolism" = "#984EA3",
  "Phosphatases"      = "#377EB8",
  "Regulatory"        = "#4DAF4A",
  "Other"             = "#999999"
)

# -------------------- Backbone: PHO regulon (Scer) ----------
scer_dir <- here("06.results", "Scer")

pho_regulon <- readxl::read_excel(
  path  = file.path(scer_dir, "20251229_PHO_Regulon.xlsx"),
  sheet = 1,
  col_types = "text"
) %>%
  transmute(
    gene_id = stringr::str_trim(as.character(gene_id)),
    gene_name = stringr::str_trim(as.character(gene_name)),
    functional_category = standardize_function(functional_category),
    category_order = suppressWarnings(as.numeric(category_order))
  ) %>%
  filter(!is.na(gene_id) & gene_id != "") %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  mutate(
    scer_gene_name = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name)
  ) %>%
  arrange(category_order, scer_gene_name)

pho_backbone <- pho_regulon %>%
  transmute(gene_id, scer_gene_name, functional_category, category_order)

# -------------------- Read LFC tidy (all species) ------------
read_lfc_species <- function(sp_code) {
  f <- list.files(here("06.results", sp_code, "paired", "combined"),
                  pattern = "LFC_tidy", full.names = TRUE)
  if (length(f) == 0) stop("No LFC_tidy file found for: ", sp_code)
  readr::read_csv(f[1], show_col_types = FALSE) %>%
    mutate(
      gene_id = stringr::str_trim(as.character(gene_id)),
      tp      = stringr::str_trim(as.character(tp))
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

# -------------------- Orthogroups: map Scer -> others --------
orthofinder_dir <- here("05.metadata", "Orthogroups")
orthogroups <- readr::read_tsv(file.path(orthofinder_dir, "Orthogroups.tsv"),
                               show_col_types = FALSE)

required_cols <- c("Orthogroup", "S_cerevisiae", "C_glabrata", "C_albicans", "K_lactis")
if (!all(required_cols %in% names(orthogroups))) {
  stop("Orthogroups.tsv missing required columns: ",
       paste(setdiff(required_cols, names(orthogroups)), collapse = ", "))
}

scer_ids_regex <- paste(na.omit(pho_backbone$gene_id), collapse = "|")

pho_orthologs <- pho_backbone %>%
  left_join(
    orthogroups %>%
      filter(stringr::str_detect(S_cerevisiae, scer_ids_regex)) %>%
      mutate(
        scer_gene_id = purrr::map_chr(S_cerevisiae, function(genes) {
          matches <- stringr::str_extract_all(genes, scer_ids_regex)[[1]]
          if (length(matches) > 0) matches[1] else NA_character_
        })
      ) %>%
      select(
        scer_gene_id,
        Orthogroup,
        Scer_genes = S_cerevisiae,
        Cgla_genes = C_glabrata,
        Calb_genes = C_albicans,
        Klac_genes = K_lactis
      ),
    by = c("gene_id" = "scer_gene_id")
  ) %>%
  rename(scer_gene_id = gene_id) %>%
  select(Orthogroup, scer_gene_name, scer_gene_id, functional_category, category_order,
         Scer_genes, Cgla_genes, Calb_genes, Klac_genes)

pho_orthologs_long <- pho_orthologs %>%
  pivot_longer(cols = ends_with("_genes"), names_to = "species", values_to = "gene_ids") %>%
  mutate(
    species = stringr::str_remove(species, "_genes"),
    species = case_when(
      species == "Scer" ~ "S_cerevisiae",
      species == "Cgla" ~ "C_glabrata",
      species == "Calb" ~ "C_albicans",
      species == "Klac" ~ "K_lactis",
      TRUE ~ species
    )
  ) %>%
  tidyr::separate_rows(gene_ids, sep = ",\\s*") %>%
  mutate(
    gene_ids = if_else(!is.na(gene_ids) & gene_ids != "",
                       stringr::str_extract(gene_ids, "^[^\\s|]+"),
                       NA_character_)
  )

# species-specific ID normalization
xp_klla_map <- readr::read_csv(
  here("06.results", "Klac", "xp_to_klla0_mapping.csv"),
  show_col_types = FALSE
)

qng_cagl_map <- readr::read_csv(
  here("06.results", "Cgla", "qng_gwk_cagl_threeway_map_complete.csv"),
  show_col_types = FALSE
)

pho_orthologs_long <- pho_orthologs_long %>%
  mutate(original_id = gene_ids) %>%
  left_join(xp_klla_map %>% select(xp_id, locus_tag),
            by = c("gene_ids" = "xp_id")) %>%
  left_join(qng_cagl_map %>% select(qng_id, gwk60_id, cagl_id),
            by = c("gene_ids" = "qng_id")) %>%
  mutate(
    gene_ids = case_when(
      species == "K_lactis"   & stringr::str_detect(original_id, "^XP_") & !is.na(locus_tag) ~ locus_tag,
      species == "C_glabrata" & stringr::str_detect(original_id, "^QNG") & !is.na(gwk60_id)  ~ gwk60_id,
      TRUE ~ gene_ids
    )
  ) %>%
  select(-locus_tag, -gwk60_id)

# ============================================================
# WRITE mapped PHO regulon CSVs 
# ============================================================
write_species_mapped_csv <- function(species_full, species_code, pho_backbone, pho_orthologs_long) {
  # species_full: "C_glabrata" etc; species_code: "Cgla" etc
  
  out_path <- here("06.results", species_code, paste0("PHO_regulon_", species_code, "_mapped.csv"))
  
  df <- pho_backbone %>%
    select(scer_gene_name, functional_category, category_order) %>%
    left_join(
      pho_orthologs_long %>%
        filter(species == species_full) %>%
        group_by(scer_gene_name) %>%
        summarise(gene_id = dplyr::first(gene_ids), .groups = "drop") %>%
        mutate(gene_id = stringr::str_trim(as.character(gene_id))),
      by = "scer_gene_name"
    ) %>%
    mutate(
      functional_category = standardize_function(functional_category),
      category_order = suppressWarnings(as.numeric(category_order))
    ) %>%
    arrange(category_order, scer_gene_name)
  
  readr::write_csv(df, out_path)
  message("Wrote: ", out_path)
  invisible(out_path)
}

# Scer mapped is basically the backbone (but we still write it for symmetry)
out_scer_mapped <- here("06.results", "Scer", "PHO_regulon_Scer_mapped.csv")
write_csv(
  pho_backbone %>% arrange(category_order, scer_gene_name),
  out_scer_mapped
)


write_species_mapped_csv("C_glabrata", "Cgla", pho_backbone, pho_orthologs_long)
write_species_mapped_csv("C_albicans", "Calb", pho_backbone, pho_orthologs_long)
write_species_mapped_csv("K_lactis",   "Klac", pho_backbone, pho_orthologs_long)

# ============================================================
# READ mapped CSVs back in + build heatmaps
# ============================================================
read_species_mapped_csv <- function(species_code) {
  p <- here("06.results", species_code, paste0("PHO_regulon_", species_code, "_mapped1.csv"))
  if (!file.exists(p)) stop("Missing mapped file: ", p)
  
  readr::read_csv(p, show_col_types = FALSE) %>%
    transmute(
      gene_id = dplyr::na_if(stringr::str_trim(as.character(gene_id)), ""),
      scer_gene_name = stringr::str_trim(as.character(scer_gene_name)),
      functional_category = standardize_function(functional_category),
      category_order = suppressWarnings(as.numeric(category_order))
    ) %>%
    arrange(category_order, scer_gene_name)
}

species_order <- c("Scer", "Cgla", "Calb", "Klac")
mapped_tables <- list(
  Scer = read_species_mapped_csv("Scer"),
  Cgla = read_species_mapped_csv("Cgla"),
  Calb = read_species_mapped_csv("Calb"),
  Klac = read_species_mapped_csv("Klac")
)

# Build per-species matrices from mapped tables
build_species_matrix <- function(species_code, mapped_gene_df, lfc_df, tp_map) {
  lfc_std <- apply_tp_map(lfc_df, tp_map)
  time_levels <- tp_map$time_std %>% as.character()
  
  inputs <- make_complete_hm_matrix(
    gene_df = mapped_gene_df,
    lfc_std = lfc_std %>% select(gene_id, time_std, log2FC),
    time_levels = time_levels,
    keep_na_if_gene_id_missing = (species_code != "Scer")
  )
  
  # enforce shared row order based on Scer backbone display order
  row_order <- mapped_tables$Scer$scer_gene_name
  inputs$mat <- inputs$mat[row_order, , drop = FALSE]
  
  # row annotation should come from the Scer backbone ordering (shared)
  row_annot <- mapped_tables$Scer %>%
    select(scer_gene_name, functional_category, category_order) %>%
    mutate(functional_category = standardize_function(functional_category)) %>%
    arrange(match(scer_gene_name, row_order))
  
  list(mat = inputs$mat, row_annot = row_annot)
}

species_mats <- list()
for (sp in species_order) {
  species_mats[[sp]] <- build_species_matrix(
    species_code = sp,
    mapped_gene_df = mapped_tables[[sp]],
    lfc_df = lfc_by_species[[sp]],
    tp_map = tp_maps[[sp]]
  )
}

# Individual heatmaps (optional interactive draw)
ht_single <- list(
  Scer = plot_pho_heatmap(species_mats$Scer$mat, species_mats$Scer$row_annot, category_colors, na_col = "grey85"),
  Cgla = plot_pho_heatmap(species_mats$Cgla$mat, species_mats$Cgla$row_annot, category_colors, na_col = "grey85"),
  Calb = plot_pho_heatmap(species_mats$Calb$mat, species_mats$Calb$row_annot, category_colors, na_col = "grey85"),
  Klac = plot_pho_heatmap(species_mats$Klac$mat, species_mats$Klac$row_annot, category_colors, na_col = "grey85")
)

# Build overview: use matrices (not Heatmap objects) to control concat + labeling
mat_list <- list(
  Scer = species_mats$Scer$mat,
  Cgla = species_mats$Cgla$mat,
  Calb = species_mats$Calb$mat,
  Klac = species_mats$Klac$mat
)

species_pretty <- c(
  Scer = "S. cerevisiae",
  Cgla = "C. glabrata",
  Calb = "C. albicans",
  Klac = "K. lactis"
)

# --- CHANGE THESE TWO ARGS in the overview call ---
pho_overview <- make_pho_species_overview_heatmap(
  mat_list = mat_list,
  row_annot = species_mats$Scer$row_annot,
  category_colors = category_colors,
  block_order = species_order,
  pretty_block = species_pretty,
  na_col = "grey85",
  column_gap_mm = 1,
  block_label_rot = 45,      # <-- 45 degrees
  block_label_font = 14,
  draw_now = FALSE
)
#  Save plots 
na_leg <- make_na_legend("No ortholog (NA)", col = "grey85")

out_pdf <- file.path(plots_dir, "PHO_regulon_4species_overview_heatmap.pdf")
pdf(out_pdf, width = 8, height = 6)
ComplexHeatmap::draw(
  pho_overview$ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right",
  heatmap_legend_list = list(na_leg),
  merge_legends = TRUE,
  padding = grid::unit(c(10, 2, 2, 6), "mm")
)
dev.off()



# ============================================================
# Panel B: Induction time is defined as the time point at which a PHO gene reaches its maximum log2FC.
#This compares temporal execution order of the PHO regulon independent of response magnitude.
# Induction time = PEAK log2FC
# ============================================================


cgla_inputs <- make_complete_hm_matrix(
  gene_df = read_csv(
    here("06.results", "Cgla", "PHO_regulon_Cgla_mapped1.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      gene_id = str_trim(as.character(gene_id)),
      scer_gene_name = paste0("Cg-", scer_gene_name)
    ),
  lfc_std = apply_tp_map(lfc_by_species$Cgla, tp_maps$Cgla),
  time_levels = tp_maps$Cgla$time_std %>% as.character(),
  keep_na_if_gene_id_missing = TRUE
)

calb_inputs <- make_complete_hm_matrix(
  gene_df = read_csv(
    here("06.results", "Calb", "PHO_regulon_Calb_mapped1.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      gene_id = str_trim(as.character(gene_id)),
      scer_gene_name = paste0("Ca-", scer_gene_name)
    ),
  lfc_std = apply_tp_map(lfc_by_species$Calb, tp_maps$Calb),
  time_levels = tp_maps$Calb$time_std %>% as.character(),
  keep_na_if_gene_id_missing = TRUE
)

klac_inputs <- make_complete_hm_matrix(
  gene_df = read_csv(
    here("06.results", "Klac", "PHO_regulon_Klac_mapped1.csv"),
    show_col_types = FALSE
  ) %>%
    mutate(
      gene_id = str_trim(as.character(gene_id)),
      scer_gene_name = paste0("Kl-", scer_gene_name)
    ),
  lfc_std = apply_tp_map(lfc_by_species$Klac, tp_maps$Klac),
  time_levels = tp_maps$Klac$time_std %>% as.character(),
  keep_na_if_gene_id_missing = TRUE
)


# -------------------- helpers -------------------------------

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

# -------------------- Panel A PHO matrices ------------------
# MUST be the same matrices used for Panel A heatmaps

mats_pho <- list(
  Scer = scer_inputs$mat,
  Cgla = cgla_inputs$mat,
  Calb = calb_inputs$mat,
  Klac = klac_inputs$mat
)

# sanity check
sapply(mats_pho, nrow)

# -------------------- compute peak-time ranks ----------------
rank_tbl <- imap_dfr(mats_pho, function(mat, sp) {
  compute_peak_time(mat) %>%
    mutate(species = sp) %>%
    rank_peak_time()
})

# Inspect earliest PHO genes per species (optional)
rank_tbl %>%
  filter(!is.na(rank)) %>%
  arrange(species, rank) %>%
  group_by(species) %>%
  slice_head(n = 8) %>%
  ungroup() %>%
  print(n = 40)

# -------------------- pairwise correlations ------------------
species <- names(mats_pho)

pair_tbl <- map_dfr(
  combn(species, 2, simplify = FALSE),
  \(p) pairwise_rank_cor(rank_tbl, p[1], p[2])
)

pair_tbl %>% arrange(desc(spearman_rho)) %>% print(n = Inf)

# -------------------- correlation heatmap --------------------
cor_mat <- pair_tbl %>%
  select(sp1, sp2, spearman_rho) %>%
  bind_rows(pair_tbl %>% transmute(sp1 = sp2, sp2 = sp1, spearman_rho)) %>%
  bind_rows(tibble(sp1 = species, sp2 = species, spearman_rho = 1)) %>%
  distinct() %>%
  mutate(
    sp1 = factor(sp1, levels = species),
    sp2 = factor(sp2, levels = species)
  )

ggplot(cor_mat, aes(sp1, sp2, fill = spearman_rho)) +
  geom_tile(color = "white") +
  coord_equal() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman \u03c1 (PHO genes)"
  ) +
  labs(
    x = NULL, y = NULL
  ) +
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
#Rank–rank scatter plots 
#computes peak-time for each gene (max log2FC across T01–T10)
#converts peak-times → ranks (earliest peak = 1)
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

# ---- build Scer vs other scatter df from rank_tbl ----
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

# ---- plot function ----
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

# ----plots (Scer vs each species) ----
pC_cgla <- plot_rankrank(rank_tbl, "Cgla")
pC_calb <- plot_rankrank(rank_tbl, "Calb")
pC_klac <- plot_rankrank(rank_tbl, "Klac")

pC_cgla
pC_calb
pC_klac

# arrange in one row (if you want a single multi-panel figure)
# install.packages("patchwork") # if needed
(pC_cgla | pC_calb | pC_klac)

# ============================================================
#  Directionality of rewiring (Δrank vs Scer)
# Δrank = rank_other - rank_Scer
# Positive = later than Scer, Negative = earlier than Scer
# ============================================================

# ---- make Δrank table for one species vs Scer ----
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

# ---- plot Δrank (bar-style, ordered by shift) ----
plot_delta_rank <- function(delta_df, other_species) {
  
  # order genes by delta (or absolute shift if you prefer)
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

# ---- build Δrank dataframes ----
d_cgla <- make_delta_rank(rank_tbl, "Cgla")
d_calb <- make_delta_rank(rank_tbl, "Calb")
d_klac <- make_delta_rank(rank_tbl, "Klac")

# Inspect tables (optional)
d_cgla %>% arrange(desc(abs_shift)) %>% print(n = 40)
d_calb %>% arrange(desc(abs_shift)) %>% print(n = 40)
d_klac %>% arrange(desc(abs_shift)) %>% print(n = 40)

# ----call individually ----
pD_cgla <- plot_delta_rank(d_cgla, "Cgla")
pD_calb <- plot_delta_rank(d_calb, "Calb")
pD_klac <- plot_delta_rank(d_klac, "Klac")

pD_cgla
pD_calb
pD_klac

(pD_cgla | pD_calb | pD_klac)




