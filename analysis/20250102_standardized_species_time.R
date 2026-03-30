---
title: "Exploratory plots from the RNA-seq data"
author: "Joshua Ayelazuno"
date: "2025-11-11"
output: html_document
---
  
  ## ================= Package setup (BiocManager style) =======================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

req <- c(
  # tidyverse core
  "readr","dplyr","tidyr","stringr","tibble","purrr","here","ggplot2",
  # plotting / export
  "UpSetR","ComplexUpset","officer","rvg"
)

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
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(purrr); library(tibble); library(here); library(ggplot2)
  library(UpSetR)
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    library(ComplexUpset)
  }
  library(officer); library(rvg)
})

## =======================
## 0) Load + prep helpers
## =======================

# ---- NEW: rules to make timepoints comparable across species ----
# Calb/Klac have extra timepoints (t0150, t0210) that I dropped for comparable plotting.
tp_drop_rules_default <- list(
  Calb = c("t0150","t0210"),
  Klac = c("t0150","t0210")
)

# ---- NEW: standardize time labels to T01..Tn (per species AFTER dropping extras) ----
make_T_labels <- function(tp_vec) {
  tp_vec <- tp_vec[order(as.integer(sub("^t", "", tp_vec)))]
  setNames(sprintf("T%02d", seq_along(tp_vec)), tp_vec)  # names are raw tXXXX
}

apply_time_filter_and_labels <- function(mat, species, drop_rules = tp_drop_rules_default) {
  drop_tp <- drop_rules[[species]]
  if (is.null(drop_tp)) drop_tp <- character(0)
  
  keep <- setdiff(colnames(mat), drop_tp)
  keep <- keep[order(as.integer(sub("^t", "", keep)))]
  mat2 <- mat[, keep, drop = FALSE]
  
  # mapping: raw -> standardized
  tp_map <- make_T_labels(colnames(mat2))
  
  list(mat = mat2, tp_map = tp_map)
}

load_lfc_mats <- function(
    fn_Scer = here("06.results","Scer","paired","combined","20251019__Scer__paired__LFC_matrix.csv"),
    fn_Cgla = here("06.results","Cgla","paired","combined","20251019__Cgla__paired__LFC_matrix.csv"),
    fn_Calb = here("06.results","Calb","paired","combined","20251019__Calb__paired__LFC_matrix.csv"),
    fn_Klac = here("06.results","Klac","paired","combined","20251019__Klac__paired__LFC_matrix.csv"),
    drop_rules = tp_drop_rules_default
) {
  to_mat <- function(df) {
    df2 <- as.data.frame(df, stringsAsFactors = FALSE)
    stopifnot("gene_id" %in% names(df2))
    rownames(df2) <- df2$gene_id
    
    tp_cols <- grep("^t\\d{4}$", names(df2), value = TRUE)
    if (length(tp_cols) == 0) stop("No tXXXX columns found in the input table.")
    
    tp_cols <- tp_cols[order(as.integer(sub("^t", "", tp_cols)))]
    mat <- as.matrix(df2[, tp_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    mat
  }
  
  Scer <- readr::read_csv(fn_Scer, show_col_types = FALSE)
  Cgla <- readr::read_csv(fn_Cgla, show_col_types = FALSE)
  Calb <- readr::read_csv(fn_Calb, show_col_types = FALSE)
  Klac <- readr::read_csv(fn_Klac, show_col_types = FALSE)
  
  mats_raw <- list(
    Scer = to_mat(Scer),
    Cgla = to_mat(Cgla),
    Calb = to_mat(Calb),
    Klac = to_mat(Klac)
  )
  
  # Apply species-specific filtering (drop extra t0150/t0210 for Calb/Klac)
  # and compute standardized T01.. labels per species.
  mats_std <- list()
  tp_maps  <- list()
  for (sp in names(mats_raw)) {
    tmp <- apply_time_filter_and_labels(mats_raw[[sp]], sp, drop_rules = drop_rules)
    mats_std[[sp]] <- tmp$mat
    tp_maps[[sp]]  <- tmp$tp_map
  }
  
  list(mats_raw = mats_raw, mats = mats_std, tp_maps = tp_maps)
}

count_deg_by_timepoint <- function(mat, thr = 1) {
  tibble::tibble(
    timepoint = colnames(mat),
    up   = colSums(mat >=  thr, na.rm = TRUE),
    down = colSums(mat <= -thr, na.rm = TRUE)
  )
}

# ---- UPDATED: prep counts uses standardized time labels for plotting ----
prep_deg_counts <- function(mats, lfc_thr = 1, tp_maps = NULL) {
  deg_counts_raw <- lapply(mats, count_deg_by_timepoint, thr = lfc_thr)
  
  # attach standardized label column (T01..)
  if (!is.null(tp_maps)) {
    deg_counts <- imap(deg_counts_raw, function(df, sp) {
      mp <- tp_maps[[sp]]
      df %>%
        mutate(time_std = unname(mp[timepoint]),
               Tnum = as.integer(sub("^T", "", time_std))) %>%
        arrange(Tnum)
    })
  } else {
    deg_counts <- deg_counts_raw
  }
  
  ylim_max <- max(unlist(lapply(deg_counts, \(d) c(d$up, d$down))), na.rm = TRUE)
  list(deg_counts = deg_counts, ylim_max = ylim_max)
}

## =======================
## 1) Boxplots (per species)
## =======================

# ---- UPDATED: boxplots show T01.. labels on x-axis ----
plot_boxplots_grid <- function(mats, tp_maps = NULL) {
  op <- par(mfrow = c(2,2), mar = c(6,4,3,1) + 0.1)
  on.exit(par(op), add = TRUE)
  
  for (sp in names(mats)) {
    df <- as.data.frame(mats[[sp]])
    
    if (!is.null(tp_maps)) {
      mp <- tp_maps[[sp]]
      colnames(df) <- unname(mp[colnames(df)])  # replace tXXXX with T01..
      # ensure order stays T01..Tn
      df <- df[, sprintf("T%02d", seq_len(ncol(df))), drop = FALSE]
    }
    
    boxplot(df,
            las = 2,
            xlab = "Timepoint (standardized)",
            ylab = "log2FC",
            main = paste0(sp, " LFC distribution"))
    abline(h = 0, col = "grey70", lty = 2)
  }
  invisible(NULL)
}

plot_boxplot_one <- function(mats, tp_maps = NULL, species = c("Scer","Cgla","Calb","Klac")) {
  species <- match.arg(species)
  df <- as.data.frame(mats[[species]])
  
  if (!is.null(tp_maps)) {
    mp <- tp_maps[[species]]
    colnames(df) <- unname(mp[colnames(df)])
    df <- df[, sprintf("T%02d", seq_len(ncol(df))), drop = FALSE]
  }
  
  boxplot(df,
          las = 2,
          xlab = "Timepoint (standardized)",
          ylab = "log2FC",
          main = paste0(species, " LFC distribution"))
  abline(h = 0, col = "grey70", lty = 2)
  invisible(NULL)
}

## =======================
## 2) Up/Down DEG barplots
## =======================

# ---- UPDATED: barplots label x-axis as T01.. and drop Calb/Klac extras already ----
plot_deg_bars_grid_gg <- function(mats, lfc_thr = 1, tp_maps = NULL,
                                  use_std_time = TRUE,
                                  col_up = "firebrick3", col_down = "dodgerblue3") {
  `%||%` <- function(x, y) if (!is.null(x)) x else y
  
  # ---- prep (your helper) ----
  prep <- prep_deg_counts(mats, lfc_thr = lfc_thr, tp_maps = tp_maps)
  deg_counts <- prep$deg_counts
  
  # ---- species display map + fixed plot order ----
  species_map <- c(
    Scer = "S.cerevisiae",
    Cgla = "C.glabrata",
    Calb = "C.albicans",
    Klac = "K.lactis"
  )
  sp_order <- c("Scer", "Cgla", "Calb", "Klac")
  sp_plot  <- intersect(sp_order, names(deg_counts))
  
  # ---- build tidy df: species x time x direction ----
  df_long <- bind_rows(lapply(sp_plot, function(sp) {
    df <- deg_counts[[sp]]
    
    tp <- if (use_std_time && "time_std" %in% names(df)) df$time_std else df$timepoint
    
    tibble(
      species_key = sp,
      species = species_map[[sp]] %||% sp,
      time = tp,
      Up = df$up,
      Down = df$down
    )
  })) %>%
    pivot_longer(c("Up", "Down"), names_to = "Direction", values_to = "Count")
  
  # ---- order timepoints like T01, T02, ... if possible ----
  if (all(grepl("^T\\d+$", df_long$time))) {
    levs <- unique(df_long$time[order(as.integer(sub("^T", "", df_long$time)))])
    df_long <- df_long %>% mutate(time = factor(time, levels = levs))
  }
  
  # ---- order facets ----
  df_long <- df_long %>%
    mutate(species = factor(species, levels = unname(species_map[sp_plot])))
  
  # ---- plot ----
  ggplot(df_long, aes(x = time, y = Count, fill = Direction)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    facet_wrap(~ species, ncol = 2) +
    scale_fill_manual(values = c(Up = col_up, Down = col_down)) +
    labs(
      x = "Timepoint",
      y = paste0("DEG count (|LFC| > ", lfc_thr, ")"),
      fill = NULL
    ) +
    theme_classic(base_size = 14) +
    theme(
      # base text (fallback)
      text = element_text(size = 14, colour = "black"),
      
      # axis titles
      axis.title.x = element_text(
        face = "bold",
        size = 14,
        colour = "black",
        margin = margin(t = 10)
      ),
      axis.title.y = element_text(
        face = "bold",
        size = 14,
        colour = "black",
        margin = margin(r = 10)
      ),
      
      # axis tick labels
      axis.text.x = element_text(
        face = "bold",
        size = 14,
        colour = "black",
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      axis.text.y = element_text(
        face = "bold",
        size = 14,
        colour = "black"
      ),
      
      # facet strips (species names)
      strip.text = element_text(
        face = "bold.italic",
        size = 14,
        colour = "black"
      ),
      
      # legend
      legend.text = element_text(
        face = "bold",
        size = 14,
        colour = "black"
      ),
      legend.title = element_blank(),
      legend.position = "bottom"
    )
}

p <- plot_deg_bars_grid_gg(mats, lfc_thr = 1, tp_maps = tp_maps)
p


## =======================
## 3) PHO84 ortholog overlay
## =======================

# ---- UPDATED: overlay uses ONLY comparable timepoints and labels x-axis as T01.. ----
plot_gene_overlay <- function(mats, tp_maps, ids, n_points = 10,
                              main = "PHO84 ortholog across -Pi time points (T01..T10)") {
  
  cols <- c(Scer="black", Cgla="dodgerblue3", Calb="firebrick3", Klac="darkgreen")
  pch  <- c(Scer=16, Cgla=17, Calb=15, Klac=18)
  
  # build per-species series in standardized T-space
  series <- list()
  for (sp in names(mats)) {
    mat <- mats[[sp]]
    gid <- ids[[sp]]
    if (is.null(gid) || !gid %in% rownames(mat)) next
    
    # already filtered (Calb/Klac extras removed), now ensure time order
    tps <- colnames(mat)
    tps <- tps[order(as.integer(sub("^t", "", tps)))]
    tps <- head(tps, n_points)
    
    y <- suppressWarnings(as.numeric(mat[gid, tps]))
    x <- seq_along(tps)
    
    # label as T01..Tn for THIS species
    Tlabs <- sprintf("T%02d", x)
    
    series[[sp]] <- list(x=x, y=y, labs=Tlabs)
  }
  
  if (length(series) == 0) stop("None of the provided gene IDs were found in the matrices.")
  
  # y-limits
  y_range <- range(unlist(lapply(series, `[[`, "y")), na.rm = TRUE)
  
  # x-axis is standardized T01..Tn (same for everyone)
  max_n <- max(vapply(series, function(s) length(s$x), integer(1)))
  x_all <- seq_len(max_n)
  
  plot(x_all, rep(NA_real_, length(x_all)), type="n",
       xaxt="n",
       xlab = "Timepoint (standardized rank)",
       ylab = "log2FC",
       ylim = y_range,
       main = main)
  
  axis(1, at = x_all, labels = sprintf("T%02d", x_all), las = 2)
  abline(h=0, col="grey60", lty=2)
  
  for (sp in names(series)) {
    lines(series[[sp]]$x, series[[sp]]$y, type="b",
          pch=pch[sp], col=cols[sp])
  }
  
  legend("topleft", legend = names(series),
         col = cols[names(series)], pch = pch[names(series)], bty="n")
  
  invisible(series)
}

## =======================
## 4) UpSet plots (Up/Down)
## =======================

# ---- UPDATED: UpSet uses standardized column names T01.. (after filtering) ----
make_upset_incidence <- function(mat, tp_map, thr = 1, direction = c("up","down")) {
  direction <- match.arg(direction)
  
  # rename columns to standardized labels before building sets
  mat2 <- mat
  colnames(mat2) <- unname(tp_map[colnames(mat2)])
  tp <- colnames(mat2)
  
  sets <- lapply(tp, function(t) {
    if (direction == "up") rownames(mat2)[which(mat2[, t] >= thr)]
    else rownames(mat2)[which(mat2[, t] <= -thr)]
  })
  names(sets) <- tp
  
  all_genes <- unique(unlist(sets))
  if (length(all_genes) == 0) return(NULL)
  
  inc <- as.data.frame(matrix(0L, nrow = length(all_genes), ncol = length(tp)))
  rownames(inc) <- all_genes
  colnames(inc) <- tp
  for (t in tp) inc[sets[[t]], t] <- 1L
  inc$gene_id <- rownames(inc)
  inc
}

plot_upset_one <- function(mat, tp_map, species = "Scer", lfc_thr = 1, direction = c("up","down"),
                           max_sets = 8, max_intersections = 20) {
  direction <- match.arg(direction)
  inc <- make_upset_incidence(mat, tp_map = tp_map, thr = lfc_thr, direction = direction)
  
  if (is.null(inc)) {
    message(species, ": no genes pass threshold for ", direction, " at |LFC|>", lfc_thr)
    return(invisible(NULL))
  }
  
  tp_cols <- setdiff(colnames(inc), "gene_id")
  
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    df <- tibble::as_tibble(inc) %>%
      dplyr::select(dplyr::all_of(tp_cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
    
    p <- ComplexUpset::upset(
      df,
      intersect = tp_cols,
      name = paste0(species, " ", direction, " (|LFC|>", lfc_thr, ")"),
      width_ratio = 0.15
    ) + ggplot2::ggtitle(paste0(species, " ", direction, " UpSet (standardized)"))
    
    print(p)
  } else {
    UpSetR::upset(
      inc[, tp_cols, drop = FALSE],
      nsets = min(max_sets, length(tp_cols)),
      nintersects = max_intersections,
      mainbar.y.label = paste0("Intersections (", direction, ")"),
      sets.x.label    = paste0("Genes per timepoint (", direction, ")")
    )
  }
  
  invisible(inc)
}

plot_upset_grid <- function(mats, tp_maps, lfc_thr = 1, direction = c("up","down")) {
  direction <- match.arg(direction)
  for (sp in names(mats)) {
    plot_upset_one(mats[[sp]], tp_map = tp_maps[[sp]], species = sp, lfc_thr = lfc_thr, direction = direction)
  }
  invisible(NULL)
}

## =======================

# Load (and automatically filter + standardize timepoints)
loaded <- load_lfc_mats()
mats    <- loaded$mats      # filtered for comparability
tp_maps <- loaded$tp_maps   # per-species mapping tXXXX -> T01..

# 1) Boxplots (T01..)
plot_boxplots_grid(mats, tp_maps = tp_maps)
# plot_boxplot_one(mats, tp_maps = tp_maps, "Cgla")

# 2) DEG bars (T01..; Calb/Klac t0150+t0210 dropped)
plot_deg_bars_grid(mats, lfc_thr = 1, tp_maps = tp_maps)

# 3) PHO84 overlay using shared comparable timepoints (T01..)
ids <- c(Scer="YML123C", Cgla="GWK60_B02321", Calb="orf19.655", Klac="KDRO_C02880")
plot_gene_overlay(mats, tp_maps = tp_maps, ids = ids)

# 4) UpSet (uses T01..)
# plot_upset_one(mats$Scer, tp_map = tp_maps$Scer, species="Scer", lfc_thr=1, direction="up")
# plot_upset_one(mats$Scer, tp_map = tp_maps$Scer, species="Scer", lfc_thr=1, direction="down")
plot_upset_grid(mats, tp_maps, lfc_thr = 1, direction = "up")
