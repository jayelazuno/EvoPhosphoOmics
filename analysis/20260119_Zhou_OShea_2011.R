# ============================================================
# Zhou & O'Shea 2011 (GSE23580) — OSR + ESR overview heatmaps
# Goal (Figure 3/4 reanalysis):
#   - Subset Pi_starvation contrasts only (No_Pi vs High_Pi)
#   - Compute Δexpression = (No_Pi − High_Pi) for each genotype
#   - Subset to curated OSR and curated ESR backbones
#   - Make overview heatmaps (replicates shown; columns split by genotype)
# Outputs:
#   06.results/Scer/GSE23580_* tables
#   08.plots/Scer/Zhou2011_GSE23580_OSR_overview.pdf
#   08.plots/Scer/Zhou2011_GSE23580_ESR_overview.pdf
# ============================================================

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

`%||%` <- function(x, y) if (!is.null(x)) x else y


# Paths
scer_dir  <- here("06.results", "Scer")
plots_dir <- here("08.plots", "Scer")

# Packages for GEO
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("GEOquery", "Biobase")) {
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
  if (!requireNamespace("Cairo", quietly = TRUE)) install.packages("Cairo")
  
}
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(Cairo)
})

# ----------------------------
# Load GEO (prefer cached objects if present)
gse_id <- "GSE23580"
eset_rds <- file.path(scer_dir, "GSE23580_eset.rds")

if (file.exists(eset_rds)) {
  eset <- readRDS(eset_rds)
} else {
  gse_obj <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, destdir = scer_dir)
  eset <- if (is.list(gse_obj)) gse_obj[[1]] else gse_obj
  saveRDS(eset, eset_rds)
}
stopifnot(inherits(eset, "ExpressionSet"))

expr_mat <- Biobase::exprs(eset)
pheno    <- Biobase::pData(eset)
feature  <- Biobase::fData(eset)

saveRDS(expr_mat, file.path(scer_dir, "GSE23580_expr_mat.rds"))
saveRDS(pheno,    file.path(scer_dir, "GSE23580_pheno.rds"))
saveRDS(feature,  file.path(scer_dir, "GSE23580_feature.rds"))

write.csv(expr_mat, file.path(scer_dir, "GSE23580_expr_mat.csv"))
write.csv(pheno,    file.path(scer_dir, "GSE23580_pheno.csv"), row.names = FALSE)
write.csv(feature,  file.path(scer_dir, "GSE23580_feature.csv"), row.names = FALSE)

# ============================================================
# Build design table from pheno (same logic you used)

extract_after <- function(x, key) {
  str_trim(str_remove(as.character(x), paste0("^", key, ":")))
}

design_tbl <- pheno %>%
  transmute(
    gsm   = geo_accession,
    title = title,
    
    genotype_ch1 = extract_after(`characteristics_ch1.1`, "genotype/variation"),
    pi_ch1 = case_when(
      str_detect(`characteristics_ch1.2`, "No Pi")   ~ "No_Pi",
      str_detect(`characteristics_ch1.2`, "High Pi") ~ "High_Pi",
      TRUE ~ NA_character_
    ),
    
    genotype_ch2 = extract_after(`characteristics_ch2.1`, "genotype/variation"),
    pi_ch2 = case_when(
      str_detect(`characteristics_ch2.2`, "No Pi")   ~ "No_Pi",
      str_detect(`characteristics_ch2.2`, "High Pi") ~ "High_Pi",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    # parse replicate label from title when present
    replicate = str_match(title, "Replicate\\s*([0-9]+)")[, 2] %||% NA_character_,
    replicate = if_else(is.na(replicate) | replicate == "", "NA", replicate),
    replicate = paste0("R", replicate),
    
    # classify comparison
    comparison_class = case_when(
      pi_ch1 != pi_ch2 & genotype_ch1 == genotype_ch2 ~ "Pi_starvation",
      pi_ch1 == pi_ch2 & genotype_ch1 != genotype_ch2 ~ "Genotype",
      TRUE ~ "Mixed"
    ),
    
    # for Pi_starvation, compute sign to enforce Δ = (No_Pi − High_Pi)
    delta_sign = case_when(
      comparison_class == "Pi_starvation" & pi_ch1 == "No_Pi"   & pi_ch2 == "High_Pi" ~  1,
      comparison_class == "Pi_starvation" & pi_ch1 == "High_Pi" & pi_ch2 == "No_Pi"   ~ -1,
      TRUE ~ NA_real_
    ),
    
    genotype = if_else(comparison_class == "Pi_starvation", genotype_ch1, NA_character_),
    genotype = str_squish(genotype)
  )

write.csv(design_tbl, file.path(scer_dir, "GSE23580_design_table.csv"), row.names = FALSE)

# ============================================================
# Subset Pi_starvation only

pi_tbl <- design_tbl %>%
  filter(comparison_class == "Pi_starvation") %>%
  filter(!is.na(delta_sign), !is.na(genotype)) %>%
  mutate(
    genotype = str_replace_all(genotype, "\\s+", " "),
    genotype = str_replace_all(genotype, "Δ", "Δ"),  # keep symbol as-is
    col_group = genotype,
    col_key   = paste0(genotype, " | ", replicate, " | ", gsm)
  )

# ensure GSMs exist as columns
pi_tbl <- pi_tbl %>% filter(gsm %in% colnames(expr_mat))
stopifnot(nrow(pi_tbl) > 0)

# ============================================================
# Feature -> gene_id mapping (robust picker)

.pick_gene_id_column <- function(feature_df) {
  nm <- names(feature_df)
  
  # try common candidates first
  preferred <- c("ORF", "orf", "ORF_NAME", "ORF.name", "GENE", "Gene", "gene", "gene_symbol",
                 "GENE_SYMBOL", "Gene.symbol", "GENE_SYMBOLS")
  hit <- preferred[preferred %in% nm]
  if (length(hit) > 0) return(hit[1])
  
  # otherwise: pick a column that contains many yeast ORF-looking entries
  score_col <- function(v) {
    v <- as.character(v)
    mean(str_detect(v, "^Y[A-P][LR][0-9]{3}[CW]$"), na.rm = TRUE)
  }
  scores <- vapply(feature_df, score_col, numeric(1))
  best <- names(scores)[which.max(scores)]
  if (is.na(best) || scores[best] < 0.05) {
    stop("Could not confidently identify a gene_id column in feature data.\n",
         "Top candidate was: ", best, " with ORF-pattern fraction=", round(scores[best], 3), "\n",
         "Available columns:\n- ", paste(nm, collapse = "\n- "))
  }
  best
}


# feature is Biobase::fData(eset)
feature_df <- as.data.frame(feature) %>%
  rownames_to_column("probe_id") %>%          # <-- THIS is the critical fix
  as_tibble(.name_repair = "minimal")

feature_map <- feature_df %>%
  mutate(
    .all_text = apply(across(everything(), as.character), 1, paste, collapse = " ")
  ) %>%
  transmute(
    probe_id,
    gene_id = str_extract(str_to_upper(.all_text), "Y[A-P][LR][0-9]{3}[CW]")
  ) %>%
  filter(!is.na(gene_id)) %>%
  distinct(probe_id, .keep_all = TRUE)

# cat("feature_map rows:", nrow(feature_map), "\n")
# cat("example probe_ids:", paste(head(feature_map$probe_id, 5), collapse = ", "), "\n")
# cat("example gene_ids:", paste(head(feature_map$gene_id, 5), collapse = ", "), "\n")

pi_expr <- as.data.frame(expr_mat) %>%
  rownames_to_column("probe_id") %>%
  pivot_longer(-probe_id, names_to = "gsm", values_to = "expr") %>%
  inner_join(pi_tbl, by = "gsm") %>%
  mutate(
    genotype = genotype_ch1,
    
    # much clearer replicate labels
    replicate = dplyr::case_when(
      stringr::str_detect(title, "Mutant Cycle Replicate\\s*1") ~ "MC1",
      stringr::str_detect(title, "Mutant Cycle Replicate\\s*2") ~ "MC2",
      stringr::str_detect(title, "Mutant Cycle Replicate\\s*3") ~ "MC3",
      stringr::str_detect(title, "Wild type no vs high Pi conditions Replicate\\s*1") ~ "WT4_1",
      stringr::str_detect(title, "Wild type no vs high Pi conditions Replicate\\s*2") ~ "WT4_2",
      stringr::str_detect(title, "Wild type no vs high Pi conditions Replicate\\s*3") ~ "WT4_3",
      stringr::str_detect(title, "Wild type no vs high Pi conditions Replicate\\s*4") ~ "WT4_4",
      TRUE ~ "UNK"
    ),
    
    col_key = paste0(genotype, " | ", replicate, " | ", gsm),
    delta = expr
  )


pi_gene <- pi_expr %>%
  inner_join(feature_map, by = "probe_id") %>%
  group_by(gene_id, genotype, replicate, col_key, gsm) %>%
  summarise(delta = mean(delta, na.rm = TRUE), .groups = "drop")

# cat("pi_gene rows:", nrow(pi_gene), "\n")
# cat("pi_gene ORFs:", dplyr::n_distinct(pi_gene$gene_id), "\n")
# cat("OSR overlap:", length(intersect(unique(pi_gene$gene_id), unique(osr_backbone$gene_id))), "\n")

write.csv(pi_gene, file.path(scer_dir, "GSE23580_Pi_starvation_gene_delta_long.csv"), row.names = FALSE)

# ============================================================
# Load curated backbones: OSR + ESR

#OSR backbone (your curated file) 
osr_backbone_path <- file.path(scer_dir, "OSR_Scer_mapped.csv")
if (!file.exists(osr_backbone_path)) stop("OSR backbone not found: ", osr_backbone_path)

standardize_function <- function(x) {
  x <- str_trim(as.character(x))
  if_else(is.na(x) | x == "", "Unknown Function", x) %>% str_to_title()
}

# OSR colors (from your Gasch OSR script)
osr_category_colors <- c(
  "Antioxidant"           = "#E41A1C",
  "Chaperon"              = "#377EB8",
  "Amino Acid Metabolism" = "#4DAF4A",
  "Carbon Metabolism"     = "#984EA3",
  "Protein Degradation"   = "#FF7F00",
  "Not Classified"        = "#999999",
  "Unknown Function"      = "#A65628"
)

osr_backbone <- readr::read_csv(
  osr_backbone_path,
  col_types = readr::cols(.default = readr::col_character())
) %>%
  transmute(
    gene_id = str_to_upper(str_trim(as.character(gene_id))),
    scer_gene_name = str_to_upper(str_trim(as.character(scer_gene_name))),
    functional_category = standardize_function(functional_category),
    category_order = suppressWarnings(as.numeric(category_order))
  ) %>%
  mutate(
    scer_gene_name = if_else(is.na(scer_gene_name) | scer_gene_name == "", gene_id, scer_gene_name),
    category_order = if_else(is.na(category_order), 999, category_order),
    functional_category = factor(functional_category, levels = names(osr_category_colors))
  ) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(category_order, scer_gene_name) %>%
  mutate(scer_gene_name = make.unique(scer_gene_name))

# ---- ESR backbone (parse from your ESR xlsx like Gasch ESR script) ----
esr_path  <- file.path(scer_dir, "20250109_ESR_clusters_UPDATED_2017.xlsx")
esr_sheet <- 2
if (!file.exists(esr_path)) stop("ESR file not found: ", esr_path)

parse_esr_NAME <- function(name_vec) {
  x <- str_squish(as.character(name_vec))
  
  gene_id <- str_extract(x, "^[A-Z0-9]{2,3}\\d{3}[A-Z]")
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

esr_set_raw <- suppressMessages(
  readxl::read_excel(esr_path, sheet = esr_sheet, col_types = "text", .name_repair = "minimal")
)
names(esr_set_raw) <- make.names(names(esr_set_raw), unique = TRUE)
if (!all(c("UID", "NAME") %in% names(esr_set_raw))) {
  stop("Expected columns UID and NAME in ESR file sheet ", esr_sheet,
       ". Found: ", paste(names(esr_set_raw), collapse = ", "))
}

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
    category_order = if_else(is.na(category_order), 999, category_order),
    functional_category = factor(functional_category, levels = esr_category_order_tbl$functional_category)
  ) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  arrange(category_order, scer_gene_name) %>%
  mutate(scer_gene_name = make.unique(scer_gene_name))

write.csv(esr_backbone, file.path(scer_dir, "Zhou2011_ESR_backbone_parsed.csv"), row.names = FALSE)

# ============================================================

fraction_induced_by_category <- function(backbone, delta_long, thresh = 1) {
  df <- delta_long %>%
    inner_join(backbone %>% select(gene_id, functional_category), by = "gene_id") %>%
    group_by(gene_id, functional_category) %>%
    summarise(delta = mean(delta, na.rm = TRUE), .groups = "drop") %>%
    mutate(induced = delta >= thresh)
  
  df %>%
    group_by(functional_category) %>%
    summarise(
      n_genes = n(),
      frac_induced = mean(induced, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(frac_induced), desc(n_genes))
}

osr_frac <- fraction_induced_by_category(osr_backbone, pi_gene, thresh = 1)
esr_frac <- fraction_induced_by_category(esr_backbone, pi_gene, thresh = 1)

write.csv(osr_frac, file.path(scer_dir, "Zhou2011_GSE23580_OSR_fraction_induced_by_category.csv"), row.names = FALSE)
write.csv(esr_frac, file.path(scer_dir, "Zhou2011_GSE23580_ESR_fraction_induced_by_category.csv"), row.names = FALSE)

#message("DONE: wrote OSR/ESR overview PDFs + fraction-induced tables to:\n  ", plots_dir)



.make_overview_from_delta <- function(delta_long,
                                      backbone,
                                      category_colors,
                                      out_prefix,
                                      lfc_limits = c(-4, 4),
                                      na_col = "grey85",
                                      show_row_names = TRUE,
                                      row_font = 14,
                                      col_font = 8,
                                      cell_w_mm = 3,
                                      cell_h_mm = 0.2,
                                      pdf_dir = NULL,
                                      pdf_w = 7,
                                      pdf_h = 10,
                                      draw_now = TRUE) {
  
  stopifnot(all(c("gene_id","scer_gene_name","functional_category") %in% names(backbone)))
  stopifnot(all(c("gene_id","genotype","replicate","gsm","col_key","delta") %in% names(delta_long)))
  
  # ----------------------------
  geno_levels <- c("Wild type", "pho2Δ", "pho4Δ", "pho2Δ pho4Δ")
  
  geno_simple <- function(x) {
    dplyr::case_when(
      x == "Wild type"     ~ "WT",
      x == "pho2Δ"         ~ "pho2Δ",
      x == "pho4Δ"         ~ "pho4Δ",
      x == "pho2Δ pho4Δ"   ~ "pho2Δ pho4Δ",
      TRUE ~ as.character(x)
    )
  }
  
  block_label <- function(genotype) {
    dplyr::case_when(
      genotype == "Wild type"     ~ "WT / Pho4⁺",
      genotype == "pho2Δ"         ~ "pho2Δ / Pho4⁺",
      genotype == "pho4Δ"         ~ "pho4Δ / Pho4⁻",
      genotype == "pho2Δ pho4Δ"   ~ "pho2Δ pho4Δ / Pho4⁻ Pho2⁻",
      TRUE ~ as.character(genotype)
    )
  }
  
  batch_label <- function(repl) {
    repl <- as.character(repl)
    dplyr::case_when(
      stringr::str_detect(repl, "^MC")   ~ "MC",
      stringr::str_detect(repl, "^WT4")  ~ "WT4",
      TRUE ~ repl
    )
  }
  
  # ----------------------------
  df <- delta_long %>%
    dplyr::inner_join(backbone %>% dplyr::select(gene_id), by = "gene_id") %>%
    dplyr::mutate(
      genotype    = factor(genotype, levels = geno_levels),
      geno_simple = geno_simple(as.character(genotype)),
      block_lab   = block_label(as.character(genotype)),
      batch       = batch_label(replicate)
    )
  
  stopifnot(nrow(df) > 0)
  
  # ----------------------------
  col_tbl <- df %>%
    dplyr::distinct(genotype, replicate, batch, gsm, col_key, geno_simple, block_lab) %>%
    dplyr::arrange(genotype, batch, replicate, gsm)
  
  col_levels <- col_tbl$col_key
  col_labels <- paste0(col_tbl$geno_simple, "-", col_tbl$batch)
  
  # ----------------------------
  sum_df <- df %>%
    dplyr::group_by(gene_id, col_key) %>%
    dplyr::summarise(delta = mean(delta, na.rm = TRUE), .groups = "drop")
  
  joined <- backbone %>% dplyr::left_join(sum_df, by = "gene_id")
  
  mat <- joined %>%
    dplyr::select(scer_gene_name, col_key, delta) %>%
    dplyr::filter(!is.na(col_key)) %>%
    tidyr::pivot_wider(names_from = col_key, values_from = delta, values_fn = mean) %>%
    tibble::column_to_rownames("scer_gene_name") %>%
    as.matrix()
  
  # ensure columns exist in correct order
  missing_cols <- setdiff(col_levels, colnames(mat))
  if (length(missing_cols) > 0) {
    for (cc in missing_cols) mat <- cbind(mat, setNames(rep(NA_real_, nrow(mat)), cc))
  }
  mat <- mat[, col_levels, drop = FALSE]
  storage.mode(mat) <- "numeric"
  
  # ----------------------------
  row_annot <- backbone %>%
    dplyr::select(scer_gene_name, functional_category, dplyr::any_of("category_order")) %>%
    dplyr::mutate(scer_gene_name = as.character(scer_gene_name)) %>%
    dplyr::filter(scer_gene_name %in% rownames(mat)) %>%
    dplyr::arrange(match(scer_gene_name, rownames(mat)))
  
  stopifnot(nrow(row_annot) == nrow(mat))
  
  row_ha <- ComplexHeatmap::rowAnnotation(
    Function = row_annot$functional_category,
    col = list(Function = category_colors),
    show_annotation_name = TRUE,
    annotation_name_gp = grid::gpar(fontsize = 10, fontface = "bold")
  )
  
  # ----------------------------
  # TOP BLOCK LABELS: ONE label per block (slice), centered, rotated, no clipping
  block_levels <- c(
    "WT / Pho4⁺",
    "pho2Δ / Pho4⁺",
    "pho4Δ / Pho4⁻",
    "pho2Δ pho4Δ / Pho4⁻ Pho2⁻"
  )
  
  split_vec <- factor(col_tbl$block_lab, levels = block_levels)
  
  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    Block = ComplexHeatmap::anno_block(
      gp = grid::gpar(col = NA),
      labels = NULL,
      panel_fun = function(index, nm) {
        grid::grid.text(
          nm,
          x = grid::unit(0.5, "npc"),
          y = grid::unit(0.1, "npc"),   # ↓ closer to heatmap (was ~0.9)
          just = "left",
          rot = 45,
          gp = grid::gpar(fontsize = 12, fontface = "bold")  # bold
        )
      }
    ),
    annotation_height = grid::unit(16, "mm")  # ↓ tighter (was ~22)
  )
  
  # ----------------------------
  ht <- ComplexHeatmap::Heatmap(
    mat,
    name = "Δ (NoPi − HighPi)",
    na_col = na_col,
    col = circlize::colorRamp2(
      c(lfc_limits[1], lfc_limits[1]/2, 0, lfc_limits[2]/2, lfc_limits[2]),
      c("#2166AC", "#4393C3", "white", "#D6604D", "#B2182B")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    
    left_annotation = row_ha,
    
    column_split   = split_vec,
    column_title   = NULL,    # remove repeated slice titles
    top_annotation = top_ha,  # keep ONE label per block
    
    show_row_names = show_row_names,
    row_names_side = "left",
    row_names_gp   = grid::gpar(fontsize = row_font, fontface = "bold"),
    
    show_column_names = TRUE,
    column_labels = col_labels,
    column_names_gp = grid::gpar(fontsize = col_font, fontface = "bold", rot = 90),
    
    width  = grid::unit(ncol(mat) * cell_w_mm, "mm"),
    height = grid::unit(nrow(mat) * cell_h_mm, "mm"),
    
    border = FALSE
  )
  
  out <- list(ht = ht, mat = mat, col_tbl = col_tbl)
  
  # ----------------------------
  if (!is.null(pdf_dir)) {
    if (!dir.exists(pdf_dir)) dir.create(pdf_dir, recursive = TRUE)
    pdf_file <- file.path(pdf_dir, paste0(out_prefix, ".pdf"))
    
    opened <- FALSE
    if (!opened && requireNamespace("Cairo", quietly = TRUE)) {
      try({ Cairo::CairoPDF(file = pdf_file, width = pdf_w, height = pdf_h); opened <- TRUE }, silent = TRUE)
    }
    if (!opened && capabilities("cairo")) {
      try({ grDevices::cairo_pdf(filename = pdf_file, width = pdf_w, height = pdf_h); opened <- TRUE }, silent = TRUE)
    }
    if (!opened) {
      grDevices::pdf(pdf_file, width = pdf_w, height = pdf_h, useDingbats = FALSE)
      opened <- TRUE
    }
    
    ComplexHeatmap::draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      merge_legends = TRUE,
      padding = grid::unit(c(10, 2, 2, 6), "mm") # enough headroom so text never clips
    )
    grDevices::dev.off()
    message("Wrote: ", pdf_file)
  }
  
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

osr_overview <- .make_overview_from_delta(
  delta_long      = pi_gene,
  backbone        = osr_backbone,
  category_colors = osr_category_colors,
  out_prefix      = "Zhou2011_GSE23580_OSR_overview",
  lfc_limits      = c(-4, 4),
  show_row_names  = TRUE,
  row_font        = 8,
  col_font        = 8,
  cell_w_mm       = 3,
  cell_h_mm       = 3,
  pdf_dir         = plots_dir,
  pdf_w           = 5.5,
  pdf_h           = 10.7,
  draw_now        = FALSE
)

esr_overview <- .make_overview_from_delta(
  delta_long      = pi_gene,
  backbone        = esr_backbone,
  category_colors = esr_category_colors,
  out_prefix      = "Zhou2011_GSE23580_ESR_overview",
  lfc_limits      = c(-4, 4),
  show_row_names  = FALSE,
  row_font        = 4,
  col_font        = 8,
  cell_w_mm       = 3,
  cell_h_mm       = 0.2,
  pdf_dir         = plots_dir,
  pdf_w           = 5.5,
  pdf_h           = 10.7,
  draw_now        = FALSE
)
