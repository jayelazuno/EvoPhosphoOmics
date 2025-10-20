# Author: Joshua Ayelazuno
# 2025-10-19
# Title: Comparative Time-Course Transcriptome Analysis of Yeast Species under Phosphate Starvation (-Pi)
# ==============================================================================
# Description:
#   This pipeline performs DESeq2-based differential expression analysis across
#   four yeast species (S. cerevisiae, C. glabrata, C. albicans, and K. lactis)
#   sampled during a phosphate starvation (-Pi) time course.
#
#   Technical replicates are collapsed per sample, and dispersion is estimated
#   using paired biological replicates (~ bio + condition). For each species and
#   timepoint (tXXXX vs t0000), the pipeline generates:
#     • Normalized and collapsed count matrices
#     • Differential expression tables (FDR < 0.01, |log2FC| > 1)
#     • Diagnostic plots (PCA, MA, Volcano, p-value, and FDR histograms)
#     • Combined multi-page summary PDFs for each plot type
#
# Output structure:
#   06.results/<Species>/paired/   – DEG tables and raw counts
#   08.plots/<Species>/paired/     – Per-timepoint and combined diagnostic PDFs
#
# Dependencies: DESeq2, apeglm, ggplot2, dplyr, readr, stringr, tibble, pdftools
# ==============================================================================



run_all_timecourses <- function(
    counts_list,                         # named list: list(Scer=Scer, Cgla=Cgla, Calb=Calb, Klac=Klac)
    species_list   = names(counts_list), # defaults to names in counts_list
    results_root   = here::here("06.results"),
    plots_root     = here::here("08.plots"),
    lfc_thr        = 1,
    fdr_cutoff     = 0.01,
    use_abs        = TRUE,
    date_tag       = format(Sys.Date(), "%Y%m%d")
) {
  ## ---- fail fast if missing ----
  req <- c("DESeq2","apeglm","ggplot2","dplyr","readr","stringr","tibble")
  miss <- req[!vapply(req, requireNamespace, FUN.VALUE = logical(1), quietly = TRUE)]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse=", "),
                         "\nPlease install before running.")
  suppressPackageStartupMessages({
    library(DESeq2); library(apeglm); library(ggplot2)
    library(dplyr);  library(readr);  library(stringr); library(tibble)
  })
  
  ## ---- ----
  fmt_thr <- function(x) gsub("\\D", "", sprintf("%.3f", x))  # 0.01 -> "001"
  
  ensure_layout <- function() {
    dirs <- c(
      file.path(results_root, species_list, "paired"),
      file.path(plots_root,   species_list, "paired"),
      file.path(results_root, "combined"),
      file.path(plots_root,   "combined")
    )
    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  }
  
  make_paths <- function(species, tp, fdr = fdr_cutoff, lfc = lfc_thr) {
    outdir <- file.path(results_root, species, "paired")
    base <- file.path(outdir, sprintf("%s__%s__paired__%s_vs_t0000__", date_tag, species, tp))
    list(
      raw_counts = paste0(base, "raw_counts.csv"),
      de_full    = paste0(base, "DE_full.csv"),
      deg_all    = paste0(base, "DEGs_all__fdr", fmt_thr(fdr), "_lfc", fmt_thr(lfc), ".csv"),
      deg_up     = paste0(base, "DEGs_up__fdr",  fmt_thr(fdr), "_lfc", fmt_thr(lfc), ".csv"),
      deg_down   = paste0(base, "DEGs_down__fdr",fmt_thr(fdr), "_lfc", fmt_thr(lfc), ".csv")
    )
  }
  make_plot_paths <- function(species, baseline, tp) {
    outdir <- file.path(plots_root, species, "paired")
    prefix <- file.path(outdir, sprintf("%s__%s__paired__", date_tag, species))
    list(
      pca      = sprintf("%sPCA_%s_vs_%s.pdf",         prefix, baseline, tp),
      ma       = sprintf("%sMA_%s_vs_%s.pdf",          prefix, baseline, tp),
      volcano  = sprintf("%sVolcano_%s_vs_%s.pdf",     prefix, baseline, tp),
      hist_p   = sprintf("%sHist_pvalue_%s_vs_%s.pdf", prefix, baseline, tp),
      hist_fdr = sprintf("%sHist_fdr_%s_vs_%s.pdf",    prefix, baseline, tp)
    )
  }
  
  run_species_paired <- function(counts_df, species) {
    if (!"gene_id" %in% names(counts_df)) stop("[", species, "] gene_id column missing.")
    rownames(counts_df) <- counts_df$gene_id
    
    # extract only this species’ columns
    count_mat <- counts_df[, setdiff(colnames(counts_df), "gene_id"), drop = FALSE]
    sel_sp <- grepl(sprintf("^%s\\.", species), colnames(count_mat))
    count_mat <- count_mat[, sel_sp, drop = FALSE]
    if (ncol(count_mat) == 0L) stop("[", species, "] No columns for this species.")
    
    # integerize counts
    count_mat <- as.data.frame(lapply(count_mat, function(x) as.integer(round(as.numeric(x)))), check.names = FALSE)
    rownames(count_mat) <- rownames(counts_df)
    
    # parse timepoints
    time_points <- gsub("^.*?(t\\d{4}).*$", "\\1", colnames(count_mat))
    if (any(!grepl("^t\\d{4}$", time_points)))
      stop("[", species, "] Unexpected timepoint token in some columns.")
    uniq_times <- unique(time_points)
    baseline <- "t0000"
    if (!baseline %in% uniq_times) stop("[", species, "] baseline t0000 not found.")
    other_times <- setdiff(sort(uniq_times), baseline)
    
    for (tp in other_times) {
      # select BOTH bios (b1 & b2) at baseline and tp (incl. both tech reps)
      sel_cols <- grep(
        paste0("^", species, "\\.(?:", baseline, "|", tp, ")\\.b[12]\\.[12]$"),
        colnames(count_mat),
        value = TRUE
      )
      if (length(sel_cols) != 8L) {
        message(sprintf("[%s %s] Expected 8 columns (2 time x 2 bio x 2 tech), got %d.",
                        species, tp, length(sel_cols)))
      }
      sub_counts <- count_mat[, sel_cols, drop = FALSE]
      
      # derive condition/bio from names; collapse tech reps via sample_id
      condition <- ifelse(grepl(paste0("\\.", baseline, "\\."), sel_cols), "untreated", "treated")
      bio       <- sub(".*\\.b([12])\\..*$", "b\\1", sel_cols) # "b1"/"b2"
      sample_id <- paste(condition, bio, sep = "_")             # untreated_b1, treated_b1, untreated_b2, treated_b2
      
      # ---- pick columns & build colData ----
      sel_cols <- grep(paste0("^", species, "\\.(?:", baseline, "|", tp, ")\\.b[12]\\.[12]$"),
                       colnames(count_mat), value = TRUE)
      sub_counts <- count_mat[, sel_cols, drop = FALSE]
      
      condition <- ifelse(grepl(paste0("\\.", baseline, "\\."), sel_cols), "untreated", "treated")
      bio       <- sub(".*\\.b([12])\\..*$", "b\\1", sel_cols)
      sample_id <- paste(condition, bio, sep = "_")  # collapse key
      
      coldata <- data.frame(
        row.names = sel_cols,
        condition = factor(condition, levels = c("untreated","treated")),
        bio       = factor(bio, levels = c("b1","b2")),
        sample    = sample_id,
        run       = sel_cols,
        stringsAsFactors = FALSE
      )
      
      # ---- assess availability & choose design ----
      present_samples <- unique(coldata$sample)                # e.g., untreated_b1, treated_b1, ...
      n_present       <- length(present_samples)
      conds_present   <- unique(coldata$condition)
      both_conditions <- length(conds_present) == 2
      
      # are both bios present in both conditions?
      has_paired <- all(c("untreated_b1","treated_b1","untreated_b2","treated_b2") %in% present_samples)
      
      if (!both_conditions || n_present < 3) {
        message(sprintf("[%s %s] Skipping: insufficient samples (conditions=%s, unique samples=%d).",
                        species, tp, paste(conds_present, collapse="/"), n_present))
        next
      }
      design_formula <- if (has_paired) ~ bio + condition else ~ condition
      if (!has_paired) {
        message(sprintf("[%s %s] Using unpaired design (~ condition). Present samples: %s",
                        species, tp, paste(sort(present_samples), collapse=", ")))
      } else {
        message(sprintf("[%s %s] Using paired design (~ bio + condition).", species, tp))
      }
      
      # ---- DESeq with tech collapse ----
      dds <- DESeqDataSetFromMatrix(countData = sub_counts, colData = coldata, design = design_formula)
      dds_col <- collapseReplicates(dds, groupby = coldata$sample, run = coldata$run)
      
      # set colData on collapsed object (S4 DataFrame, ordered to columns)
      library(S4Vectors)
      sample_map <- unique(coldata[, c("sample","condition","bio")]); rownames(sample_map) <- sample_map$sample
      sample_map2 <- sample_map[colnames(dds_col), c("condition","bio")]
      sample_map2$condition <- factor(sample_map2$condition, levels = c("untreated","treated"))
      sample_map2$bio       <- factor(sample_map2$bio,       levels = c("b1","b2"))
      colData(dds_col) <- S4Vectors::DataFrame(sample_map2)
      
      dds_col <- DESeq(dds_col)
      coef_name <- "condition_treated_vs_untreated"
      res_lfc <- lfcShrink(dds_col, coef = coef_name, type = "apeglm")
      
      # write results
      paths <- make_paths(species, tp)
      raw_counts <- as.data.frame(counts(dds_col)) %>% rownames_to_column("gene_id")
      readr::write_csv(raw_counts, paths$raw_counts)
      
      de_full <- as.data.frame(res_lfc) %>%
        rownames_to_column("gene_id") %>%
        mutate(gene_id = sub("^gene-", "", gene_id))
      readr::write_csv(de_full, paths$de_full)
      
      is_sig <- with(de_full,
                     !is.na(padj) & padj < fdr_cutoff &
                       if (use_abs) abs(log2FoldChange) > lfc_thr else log2FoldChange > lfc_thr)
      deg_all  <- de_full[is_sig, , drop = FALSE]
      deg_up   <- deg_all[deg_all$log2FoldChange >  0, , drop = FALSE]
      deg_down <- deg_all[deg_all$log2FoldChange <  0, , drop = FALSE]
      readr::write_csv(deg_all,  paths$deg_all)
      readr::write_csv(deg_up,   paths$deg_up)
      readr::write_csv(deg_down, paths$deg_down)
      
      # plots
      ppaths <- make_plot_paths(species, baseline, tp)
      rld <- rlog(dds_col, blind = TRUE)
      rldmat <- assay(rld)
      pca <- prcomp(t(rldmat))
      pct <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)))
      sample_info <- as.data.frame(colData(dds_col))
      sample_info$sample <- rownames(sample_info)
      pca_df <- data.frame(
        PC1    = pca$x[, 1],
        PC2    = pca$x[, 2],
        sample = rownames(pca$x),
        group  = sample_info$condition
      )
      p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
        geom_point(size = 5, alpha = 0.85) +
        labs(title = sprintf("%s paired: %s vs %s", species, baseline, tp),
             x = paste0("PC1 (", pct[1], "%)"),
             y = paste0("PC2 (", pct[2], "%)")) +
        theme_minimal(base_size = 16) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
      ggsave(ppaths$pca, p_pca, width = 8, height = 6)
      
      pdf(ppaths$ma, width = 7, height = 5)
      plotMA(res_lfc, main = sprintf("MA - %s paired: %s vs %s", species, baseline, tp), cex = 0.5)
      abline(h = c(-lfc_thr, lfc_thr), col = "dodgerblue3", lty = 2)
      dev.off()
      
      tab <- data.frame(logFC = res_lfc$log2FoldChange, negLogFDR = -log10(res_lfc$padj))
      sig_flag <- with(tab, !is.na(negLogFDR) &
                         (if (use_abs) abs(logFC) > lfc_thr else logFC > lfc_thr) &
                         negLogFDR > -log10(fdr_cutoff))
      pdf(ppaths$volcano, width = 7, height = 5)
      plot(tab$logFC, tab$negLogFDR, pch = 16, cex = 0.5,
           xlab = expression(bold(log[2]~fold~change)),
           ylab = expression(bold(-log[10]~FDR)),
           main = sprintf("Volcano - %s paired: %s vs %s", species, baseline, tp))
      points(tab$logFC[sig_flag], tab$negLogFDR[sig_flag], col = "red", pch = 16, cex = 0.6)
      abline(h = -log10(fdr_cutoff), col = "darkgreen", lty = 2)
      abline(v = c(-lfc_thr, lfc_thr), col = "dodgerblue3", lty = 2)
      dev.off()
      
      pdf(ppaths$hist_p, width = 7, height = 5)
      hist(res_lfc$pvalue, breaks = 50, col = "grey",
           main = sprintf("Raw p-values: %s paired %s vs %s", species, baseline, tp),
           xlab = "p-value")
      dev.off()
      
      pdf(ppaths$hist_fdr, width = 7, height = 5)
      hist(res_lfc$padj, breaks = 50, col = "grey",
           main = sprintf("Adjusted p-values (FDR): %s paired %s vs %s", species, baseline, tp),
           xlab = "FDR (padj)")
      dev.off()
      
      message(sprintf("Finished %s paired: %s vs %s", species, baseline, tp))
    }
    
    invisible(TRUE)
  }
  
  ## ---- ensure dirs, then run ----
  ensure_layout()
  for (sp in species_list) {
    if (!sp %in% names(counts_list)) {
      warning(sprintf("Skipping %s: not found in counts_list"), immediate. = TRUE)
      next
    }
    message(sprintf("==> Running %s ...", sp))
    run_species_paired(counts_df = counts_list[[sp]], species = sp)
  }
  invisible(TRUE)
}


# 1) Read each species matrix (already renamed to .b1/.b2 with .1/.2)
Scer <- read_tsv(here("06.results","Scer","20251014_Scer_counts.tsv"), show_col_types = FALSE)
Cgla <- read_tsv(here("06.results","Cgla","20251014_Cgla_counts.tsv"), show_col_types = FALSE)
Calb <- read_tsv(here("06.results","Calb","20251014_Calb_counts.tsv"), show_col_types = FALSE)
Klac <- read_tsv(here("06.results","Klac","20251014_Klac_counts.tsv"), show_col_types = FALSE)

# 2) Put them in a named list (names must match species codes)
counts_list <- list(Scer = Scer, Cgla = Cgla, Calb = Calb, Klac = Klac)

# 3) Run (species_list defaults to names(counts_list))
run_all_timecourses(
  counts_list  = counts_list,
  results_root = here("06.results"),
  plots_root   = here("08.plots"),
  lfc_thr      = 1,
  fdr_cutoff   = 0.01,
  use_abs      = TRUE
)


