
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

load_lfc_mats <- function(
    fn_Scer = here("06.results","Scer","paired","combined","20251019__Scer__paired__LFC_matrix.csv"),
    fn_Cgla = here("06.results","Cgla","paired","combined","20251019__Cgla__paired__LFC_matrix.csv"),
    fn_Calb = here("06.results","Calb","paired","combined","20251019__Calb__paired__LFC_matrix.csv"),
    fn_Klac = here("06.results","Klac","paired","combined","20251019__Klac__paired__LFC_matrix.csv")
) {
  to_mat <- function(df) {
    df2 <- as.data.frame(df, stringsAsFactors = FALSE)
    stopifnot("gene_id" %in% names(df2))
    rownames(df2) <- df2$gene_id
    
    tp_cols <- grep("^t\\d{4}$", names(df2), value = TRUE)
    if (length(tp_cols) == 0) stop("No tXXXX columns found in the input table.")
    
    tp_cols <- tp_cols[order(as.integer(sub("^t", "", tp_cols)))]
    as.matrix(df2[, tp_cols, drop = FALSE])
  }
  
  Scer <- readr::read_csv(fn_Scer, show_col_types = FALSE)
  Cgla <- readr::read_csv(fn_Cgla, show_col_types = FALSE)
  Calb <- readr::read_csv(fn_Calb, show_col_types = FALSE)
  Klac <- readr::read_csv(fn_Klac, show_col_types = FALSE)
  
  mats <- list(
    Scer = to_mat(Scer),
    Cgla = to_mat(Cgla),
    Calb = to_mat(Calb),
    Klac = to_mat(Klac)
  )
  mats
}

count_deg_by_timepoint <- function(mat, thr = 1) {
  tibble::tibble(
    timepoint = colnames(mat),
    up   = colSums(mat >=  thr, na.rm = TRUE),
    down = colSums(mat <= -thr, na.rm = TRUE)
  )
}

prep_deg_counts <- function(mats, lfc_thr = 1) {
  deg_counts <- lapply(mats, count_deg_by_timepoint, thr = lfc_thr)
  ylim_max <- max(unlist(lapply(deg_counts, \(d) c(d$up, d$down))), na.rm = TRUE)
  list(deg_counts = deg_counts, ylim_max = ylim_max)
}


## =======================
## 1) Boxplots (per species)
## =======================

plot_boxplots_grid <- function(mats) {
  op <- par(mfrow = c(2,2), mar = c(6,4,3,1) + 0.1)
  on.exit(par(op), add = TRUE)
  
  for (sp in names(mats)) {
    boxplot(as.data.frame(mats[[sp]]),
            las = 2,
            xlab = "Timepoint",
            ylab = "log2FC",
            main = paste0(sp, " LFC distribution"))
    abline(h = 0, col = "grey70", lty = 2)
  }
  invisible(NULL)
}

plot_boxplot_one <- function(mats, species = c("Scer","Cgla","Calb","Klac")) {
  species <- match.arg(species)
  boxplot(as.data.frame(mats[[species]]),
          las = 2,
          xlab = "Timepoint",
          ylab = "log2FC",
          main = paste0(species, " LFC distribution"))
  abline(h = 0, col = "grey70", lty = 2)
  invisible(NULL)
}


## =======================
## 2) Up/Down DEG barplots
## =======================

plot_deg_bars_one <- function(df_counts, species, lfc_thr = 1, ylim_max = NULL,
                              col_up = "firebrick3", col_down = "dodgerblue3") {
  m <- rbind(df_counts$up, df_counts$down)
  colnames(m) <- df_counts$timepoint
  if (is.null(ylim_max)) ylim_max <- max(m, na.rm = TRUE)
  
  op <- par(mar = c(7,4,3,1) + 0.1)
  on.exit(par(op), add = TRUE)
  
  barplot(m, beside = TRUE, col = c(col_up, col_down), border = NA,
          ylim = c(0, ylim_max * 1.05), las = 2,
          ylab = "DEG count", xlab = "Timepoint",
          main = paste0(species, " DEGs (|LFC| > ", lfc_thr, ")"))
  legend("topleft", fill = c(col_up, col_down), legend = c("Up", "Down"), bty = "n")
  abline(h = 0, col = "grey70", lty = 2)
  invisible(NULL)
}

plot_deg_bars_grid <- function(mats, lfc_thr = 1) {
  prep <- prep_deg_counts(mats, lfc_thr = lfc_thr)
  deg_counts <- prep$deg_counts
  ylim_max   <- prep$ylim_max
  
  op <- par(mfrow = c(2,2), mar = c(7,4,3,1) + 0.1)
  on.exit(par(op), add = TRUE)
  
  for (sp in names(mats)) {
    plot_deg_bars_one(deg_counts[[sp]], sp, lfc_thr = lfc_thr, ylim_max = ylim_max)
  }
  invisible(prep)
}


## =======================
## 3) PHO84 ortholog overlay
## =======================

plot_gene_overlay <- function(mats, ids, main = "PHO84 ortholog across -Pi time points") {
  tp_common <- Reduce(intersect, lapply(mats, colnames)) |> sort()
  ys <- lapply(names(mats), function(sp) {
    mat <- mats[[sp]]
    gid <- ids[[sp]]
    if (is.null(gid) || !gid %in% rownames(mat)) return(NULL)
    suppressWarnings(as.numeric(mat[gid, tp_common]))
  })
  names(ys) <- names(mats)
  ys <- ys[!vapply(ys, is.null, logical(1))]
  if (length(ys) == 0) stop("None of the provided gene IDs were found in the matrices.")
  
  y_range <- range(unlist(ys), na.rm = TRUE)
  x_vals <- seq_along(tp_common)
  cols <- c(Scer="black", Cgla="dodgerblue3", Calb="firebrick3", Klac="darkgreen")
  pch  <- c(Scer=16, Cgla=17, Calb=15, Klac=18)
  
  plot(x_vals, rep(NA_real_, length(x_vals)), type="n",
       xaxt="n", xlab="Timepoint (vs t0000)", ylab="log2FC",
       ylim=y_range, main=main)
  axis(1, at = x_vals, labels = tp_common, las = 2)
  abline(h=0, col="grey60", lty=2)
  
  for (sp in names(ys)) {
    lines(x_vals, ys[[sp]], type="b", pch=pch[sp], col=cols[sp])
  }
  legend("topleft", legend = names(ys), col = cols[names(ys)], pch = pch[names(ys)], bty="n")
  invisible(NULL)
}


## =======================
## 4) UpSet plots (Up/Down)
## =======================

.order_tp_cols <- function(tp_cols) {
  # expects names like "t0000", "t0015", "t0120", etc.
  tp_num <- suppressWarnings(as.integer(gsub("^t", "", tp_cols)))
  tp_cols[order(tp_num, na.last = TRUE)]
}

# direction = up/down/both
make_upset_incidence <- function(mat, thr = 1, direction = c("up","down","both")) {
  direction <- match.arg(direction)
  tp <- colnames(mat)
  
  sets <- lapply(tp, function(t) {
    if (direction == "up") {
      rownames(mat)[which(mat[, t] >= thr)]
    } else if (direction == "down") {
      rownames(mat)[which(mat[, t] <= -thr)]
    } else { # both
      rownames(mat)[which(abs(mat[, t]) >= thr)]
    }
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

plot_upset_one <- function(mat,
                           species = "Scer",
                           lfc_thr = 1,
                           direction = c("up","down","both"),
                           max_sets = 8,
                           max_intersections = 40,
                           save_width = 16,
                           save_height = 6) {
  
  direction <- match.arg(direction)
  inc <- make_upset_incidence(mat, thr = lfc_thr, direction = direction)
  
  if (is.null(inc)) {
    message(species, ": no genes pass threshold for ", direction, " at |LFC|>", lfc_thr)
    return(invisible(NULL))
  }
  
  tp_cols <- colnames(inc)[grepl("^t\\d{4}$", colnames(inc))]
  
  # force chronological order
  tp_cols <- .order_tp_cols(tp_cols)
  
  # Prefer ComplexUpset if installed (ggplot-based)
  if (requireNamespace("ComplexUpset", quietly = TRUE)) {
    
    df <- tibble::as_tibble(inc) %>%
      dplyr::select(dplyr::all_of(tp_cols)) %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
    
    dir_label <- switch(direction,
                        up   = "Up",
                        down = "Down",
                        both = "Up+Down (|LFC| threshold)")
    
    # Wider bars, spacing, numbers, no intersection labels at bottom
    p <- ComplexUpset::upset(
      df,
      intersect = tp_cols,
      name = paste0(species, " ", dir_label, " (|LFC|>", lfc_thr, ")"),
      width_ratio = 0.4,
      sort_sets = FALSE,
      sort_intersections_by = "cardinality",
      n_intersections = max_intersections,
      base_annotations = list(
        'Intersection size' = ComplexUpset::intersection_size(
          counts = TRUE,
          bar_number_threshold = 1,
          width = 0.6,                        # 0.6 = wider spacing, thinner bars
          text = list(size = 3.5)
        )
      )
    ) +
      ggplot2::ggtitle(paste0(species, " UpSet: ", dir_label)) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 14),
        axis.text.x = ggplot2::element_blank(),      # removes intersection labels
        axis.ticks.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(face = "bold", size = 11),
        axis.title.x = ggplot2::element_text(face = "bold", size = 12),
        axis.title.y = ggplot2::element_text(face = "bold", size = 12),
        strip.text = ggplot2::element_text(face = "bold", size = 14),  # timepoint labels
        panel.spacing.x = grid::unit(0.6, "cm"),
        panel.spacing.y = grid::unit(0.4, "cm")
      )
    
    print(p)
    return(invisible(p))
    
  } else {
    # UpSetR fallback
    UpSetR::upset(
      inc[, tp_cols, drop = FALSE],
      nsets = min(max_sets, length(tp_cols)),
      nintersects = max_intersections,
      mainbar.y.label = paste0("Intersections (", direction, ")"),
      sets.x.label    = paste0("Genes per timepoint (", direction, ")"),
      mb.ratio = c(0.65, 0.35),
      text.scale = c(1.6, 1.6, 1.4, 1.4, 1.6, 1.4),
      show.numbers = "yes"
    )
  }
  
  invisible(inc)
}

# convenience wrapper: produce all 3 versions (up, down, combined)
plot_upset_all <- function(mat, species = "Scer", lfc_thr = 1,
                           max_sets = 8, max_intersections = 40) {
  plot_upset_one(mat, species = species, lfc_thr = lfc_thr, direction = "up",
                 max_sets = max_sets, max_intersections = max_intersections)
  plot_upset_one(mat, species = species, lfc_thr = lfc_thr, direction = "down",
                 max_sets = max_sets, max_intersections = max_intersections)
  plot_upset_one(mat, species = species, lfc_thr = lfc_thr, direction = "both",
                 max_sets = max_sets, max_intersections = max_intersections)
  invisible(NULL)
}




plot_upset_one(mats$Scer, species = "Scer", lfc_thr = 1, direction = "up")

plot_upset_one(mats$Scer, species = "Scer", lfc_thr = 1, direction = "down")
plot_upset_one(mats$Scer, species = "Scer", lfc_thr = 1, direction = "both")

# --- C. glabrata ---
plot_upset_one(mats$Cgla, species = "Cgla", lfc_thr = 1, direction = "up")
plot_upset_one(mats$Cgla, species = "Cgla", lfc_thr = 1, direction = "down")
plot_upset_one(mats$Cgla, species = "Cgla", lfc_thr = 1, direction = "both")

# --- C. albicans ---
plot_upset_one(mats$Calb, species = "Calb", lfc_thr = 1, direction = "up")
plot_upset_one(mats$Calb, species = "Calb", lfc_thr = 1, direction = "down")
plot_upset_one(mats$Calb, species = "Calb", lfc_thr = 1, direction = "both")

# --- K. lactis ---
plot_upset_one(mats$Klac, species = "Klac", lfc_thr = 1, direction = "up")
plot_upset_one(mats$Klac, species = "Klac", lfc_thr = 1, direction = "down")
plot_upset_one(mats$Klac, species = "Klac", lfc_thr = 1, direction = "both")

# Save with wide dimensions for proper bar width + spacing
p <- plot_upset_one(mats$Scer, species = "Scer", lfc_thr = 1, direction = "up")
ggsave("Scer_upset_up.png", plot = p, width = 16, height = 6, dpi = 300)
# --- S. cerevisiae ---