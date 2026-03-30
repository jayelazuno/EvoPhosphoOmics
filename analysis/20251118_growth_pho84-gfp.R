
## ================= Load Packages =======================
library(growthcurver)
library(tidyverse)
library(cowplot)
library(here)
library(gt)
library(hms)  # Explicitly load hms package

## ================= Custom Theme & Colors =======================
custom_colors <- c("S.cerevisiae" = "darkorchid4", 
                   "C.glabrata" = "chartreuse4", 
                   "C.albicans" = "brown3", 
                   "K.lactis" = "aquamarine3")

publication_theme <- theme_minimal() +
  theme(
    legend.text = element_text(face = "italic", size = 14),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

## ================= Load & Process Data =======================
# Read data with duplicate column handling
od_raw <- read_csv(here("05.metadata","20230923-growth-PHO84-gfp-data", "20230907_species_growth.csv"))

od <- od_raw %>%
mutate(
  # Handle H:MM:SS format (e.g., "0:14:10")
  # Convert to character first if needed, then parse as hms
  Time_char = as.character(Time),
  Time_hms = hms::as_hms(Time_char),
  Time_hrs = as.numeric(Time_hms) / 3600  # hms stores as seconds, convert to hours
) %>%
  select(-Time, -Time_char, -Time_hms) %>%   
  
pivot_longer(
  cols = -Time_hrs,
  names_to = "variable",
  values_to = "OD600"
) %>%
  
mutate(
  variable = str_remove(variable, "_$"),
  Species  = str_extract(variable, "^[A-Za-z]+\\.[a-z]+"),
  Pi_mM    = str_extract(variable, "[0-9.]+mM"),
  Pi_mM    = as.numeric(str_remove(Pi_mM, "mM"))
) %>%
  
  select(Time_hrs, Species, Pi_mM, OD600) %>%
  filter(!is.na(Species), !is.na(Pi_mM))


# Set species order
species_order <- c("S.cerevisiae", "C.glabrata", "C.albicans", "K.lactis")
od$Species <- factor(od$Species, levels = species_order)


#### ===== Growth Curves =========
growth_curve_plot <- od %>%
  mutate(Pi_Label = paste0(Pi_mM, " mM Pi")) %>%
  ggplot(aes(x = Time_hrs, y = OD600, color = Species)) +
  geom_point(size = 0.7, alpha = 0.6, position = position_jitter(width = 0.2)) +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = " ",
    x = "Time (hours)",
    y = expression(OD[600]),
    color = "Species"
  ) +
  facet_wrap(~ Pi_Label) +
  publication_theme +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"))

growth_curve_plot

## ================= Fit Growth Models =======================
growth_params <- od %>%
  group_by(Species, Pi_mM) %>%
  group_modify(~ {
    time_hrs <- .x$Time_hrs
    OD_values <- .x$OD600
    
    tryCatch({
      gc_fit <- SummarizeGrowth(time_hrs, OD_values)
      
      tibble(
        Carrying_Capacity = gc_fit$vals$k,
        Growth_Rate = gc_fit$vals$r,
        Doubling_Time = gc_fit$vals$t_gen,
        Fit_Quality = gc_fit$vals$sigma,
        Note = gc_fit$vals$note
      )
    }, error = function(e) {
      tibble(
        Carrying_Capacity = NA, Growth_Rate = NA, Doubling_Time = NA,
        Fit_Quality = NA, Note = paste("Error:", e$message)
      )
    })
  }) %>%
  ungroup()

# Set species order in growth parameters
growth_params$Species <- factor(growth_params$Species, levels = species_order)

## ================= Publication Parameter Plots =======================
# Carrying Capacity Plot
k_plot <- ggplot(growth_params, aes(x = Species, y = Carrying_Capacity, fill = Species)) +
  geom_col(position = "dodge", alpha = 0.8, show.legend = FALSE, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Carrying Capacity",
    y = expression(bold(OD[600])),
    x = ""
  ) +
  facet_wrap(~ paste0(Pi_mM, " mM Pi")) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Growth Rate Plot
r_plot <- ggplot(growth_params, aes(x = Species, y = Growth_Rate, fill = Species)) +
  geom_col(position = "dodge", alpha = 0.8, show.legend = FALSE, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Growth Rate",
    y = expression(bold("r (h"^{-1}*")")),
    x = ""
  ) +
  facet_wrap(~ paste0(Pi_mM, " mM Pi")) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Doubling Time Plot
doubling_plot <- ggplot(growth_params, aes(x = Species, y = Doubling_Time, fill = Species)) +
  geom_col(position = "dodge", alpha = 0.8, show.legend = FALSE, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(
    title = "Doubling Time",
    y = "Time (hours)",
    x = ""
  ) +
  facet_wrap(~ paste0(Pi_mM, " mM Pi")) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Display parameter plots
k_plot
r_plot
doubling_plot

## ================= Publication Quality Table =======================
growth_table <- growth_params %>%
  mutate(across(c(Carrying_Capacity, Growth_Rate, Doubling_Time), 
                ~ round(., 3))) %>%
  select(Species, Pi_mM, Carrying_Capacity, Growth_Rate, Doubling_Time) %>%
  pivot_wider(
    names_from = Pi_mM,
    values_from = c(Carrying_Capacity, Growth_Rate, Doubling_Time),
    names_glue = "{.value}_{Pi_mM}mM"
  )


# Arrange all four plots in 2x2 grid
combined_plot <- plot_grid(
  growth_curve_plot,
  k_plot,
  r_plot, 
  doubling_plot,
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  label_size = 14
)

# Display combined plot
combined_plot

# Create formatted table
formatted_table <- growth_table %>%
  select(-Doubling_Time_0mM, -Doubling_Time_7.3mM) %>%  # Remove doubling time columns
  gt() %>%
  tab_header(
    title = "Growth Parameters Across Species and Phosphate Conditions",
    subtitle = "Carrying capacity (OD₆₀₀) and growth rate (r, h⁻¹)"
  ) %>%
  fmt_number(decimals = 3) %>%
  cols_label(
    Species = "Species",
    Carrying_Capacity_0mM = "K (0 mM)",
    Carrying_Capacity_7.3mM = "K (7.3 mM)", 
    Growth_Rate_0mM = "r (0 mM)",
    Growth_Rate_7.3mM = "r (7.3 mM)"
  ) %>%
  tab_spanner(
    label = "Carrying Capacity",
    columns = c(Carrying_Capacity_0mM, Carrying_Capacity_7.3mM)
  ) %>%
  tab_spanner(
    label = "Growth Rate (h⁻¹)", 
    columns = c(Growth_Rate_0mM, Growth_Rate_7.3mM)
  ) %>%
  # Add bold headers and bold italic species names
  tab_style(
    style = cell_text(weight = "bold"),
    locations = list(
      cells_column_labels(),
      cells_column_spanners(),
      cells_title()
    )
  ) %>%
  tab_style(
    style = cell_text(style = "italic", weight = "bold"),
    locations = cells_body(columns = Species)
  ) %>%
  tab_options(
    table.font.size = 12,
    heading.title.font.size = 16,
    heading.subtitle.font.size = 14
  )

# Display table
formatted_table
##############################

flow <- read_csv(here("05.metadata","20230923-growth-PHO84-gfp-data","20230907_PHO84-gfp_flow.csv"), 
                 name_repair = "unique") %>%
  # Create unique column names by combining with index
  rename_with(~ {
    base_names <- gsub("\\.\\.\\d+", "", .x)
    # Add replicate numbers to make them unique
    make.unique(base_names, sep = "_rep")
  }) %>%
  # Reshape using pivot_longer
  pivot_longer(
    cols = -`Time(hrs)`,
    names_to = "variable",
    values_to = "value"
  ) %>%
  na.omit() %>%
  # Extract species name (everything before "PHO")
  mutate(Species = str_extract(variable, "^[A-Za-z]+\\.[a-z]+")) %>%
  select(Time_hrs = `Time(hrs)`, Species, `Log10(PHO84pr-EGFP(a.u.))` = value) %>%
  filter(!is.na(Species))

# Set species order to match your other plots
flow$Species <- factor(flow$Species, levels = species_order)

# Create the flow cytometry plot with custom theme and colors
fl <- flow %>%
  ggplot(aes(x = Time_hrs, y = `Log10(PHO84pr-EGFP(a.u.))`, group = Species)) +
  geom_point(aes(color = Species), size = 1.5, alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", aes(color = Species), linewidth = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "PHO84 Promoter Activity Over Time",
    x = "Time (hours)",
    y = expression(bold("Log"[10]*"(PHO84pr-EGFP (a.u.))")),
    color = "Species"
  ) +
  publication_theme +
  theme(axis.text.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13))
fl

#################################################
# normalizing PHO84 - expression 
##################################################
# nomralized to 1 
flow_normalized <- flow %>%
  group_by(Species) %>%
  mutate(Normalized_GFP = (`Log10(PHO84pr-EGFP(a.u.))` - min(`Log10(PHO84pr-EGFP(a.u.))`)) / 
           (max(`Log10(PHO84pr-EGFP(a.u.))`) - min(`Log10(PHO84pr-EGFP(a.u.))`))) %>%
  ungroup()

fl_norm <- flow_normalized %>%
  ggplot(aes(x = Time_hrs, y = Normalized_GFP, group = Species)) +
  geom_point(aes(color = Species), size = 1.5, alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", aes(color = Species), linewidth = 1.2) +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "Normalized PHO84 Promoter Activity Over Time",
    x = "Time (hours)",
    y = "Normalized PHO84pr-EGFP",
    color = "Species"
  ) +
  publication_theme +
  theme(axis.text.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13))
fl_norm

flow_foldchange <- flow %>%
  group_by(Species) %>%
  mutate(Baseline = mean(`Log10(PHO84pr-EGFP(a.u.))`[Time_hrs == 0], na.rm = TRUE),
         Fold_Change = `Log10(PHO84pr-EGFP(a.u.))` - Baseline) %>%
  ungroup()

######## normalized to fold change 
fl_fc <- flow_foldchange %>%
  ggplot(aes(x = Time_hrs, y = Fold_Change, group = Species)) +
  geom_point(aes(color = Species), size = 1.5, alpha = 0.6) +
  stat_summary(fun = "mean", geom = "line", aes(color = Species), linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = custom_colors) +
  labs(
    title = "PHO84 Promoter Induction Over Time",
    x = "Time (hours)",
    y = expression(bold(Delta*"Log"[10]*"(PHO84pr-EGFP)")),
    color = "Species"
  ) +
  publication_theme +
  theme(axis.text.x = element_text(face = "bold", size = 13),
        axis.text.y = element_text(face = "bold", size = 13))
fl_fc

## ================= PHO84 Expression Plot Script =======================

# Define file paths
fn_Scer <- here("06.results","Scer","paired","combined","20251019__Scer__paired__LFC_matrix.csv")
fn_Cgla <- here("06.results","Cgla","paired","combined","20251019__Cgla__paired__LFC_matrix.csv")
fn_Calb <- here("06.results","Calb","paired","combined","20251019__Calb__paired__LFC_matrix.csv")
fn_Klac <- here("06.results","Klac","paired","combined","20251019__Klac__paired__LFC_matrix.csv")

# Define PHO84 ortholog IDs
ids <- c(Scer = "YML123C", 
         Cgla = "GWK60_B02321", 
         Calb = "orf19.655", 
         Klac = "KDRO_C02880")

# Load data
Scer <- read_csv(fn_Scer, show_col_types = FALSE)
Cgla <- read_csv(fn_Cgla, show_col_types = FALSE)
Calb <- read_csv(fn_Calb, show_col_types = FALSE)
Klac <- read_csv(fn_Klac, show_col_types = FALSE)

# Convert to matrices
to_mat <- function(df) {
  df2 <- as.data.frame(df, stringsAsFactors = FALSE)
  rownames(df2) <- df2$gene_id
  tp_cols <- grep("^t\\d{4}$", names(df2), value = TRUE)
  tp_cols <- tp_cols[order(as.integer(sub("^t", "", tp_cols)))]
  as.matrix(df2[, tp_cols, drop = FALSE])
}

Scer_mat <- to_mat(Scer)
Cgla_mat <- to_mat(Cgla)
Calb_mat <- to_mat(Calb)
Klac_mat <- to_mat(Klac)

mats <- list(Scer = Scer_mat, Cgla = Cgla_mat, Calb = Calb_mat, Klac = Klac_mat)

# Get ALL time points from each species (don't restrict to common ones)
pho84_data <- map_dfr(names(mats), function(sp) {
  mat <- mats[[sp]]
  gid <- ids[[sp]]
  
  if (gid %in% rownames(mat)) {
    # Get all timepoints for this species
    all_tps <- colnames(mat)
    values <- suppressWarnings(as.numeric(mat[gid, all_tps]))
    
    # Convert species codes to full names for coloring
    species_name <- case_when(
      sp == "Scer" ~ "S.cerevisiae",
      sp == "Cgla" ~ "C.glabrata", 
      sp == "Calb" ~ "C.albicans",
      sp == "Klac" ~ "K.lactis"
    )
    
    tibble(
      Species = species_name,
      Timepoint = all_tps,
      LFC = values,
      Time_numeric = as.numeric(sub("^t", "", all_tps))
    )
  } else {
    NULL
  }
})

# MANUALLY ADD TIME 0 POINT WITH LFC = 0 FOR EACH SPECIES
pho84_data_with_zero <- pho84_data %>%
  group_by(Species) %>%
  group_modify(~ {
    # Check if t0000 already exists
    if (!"t0000" %in% .x$Timepoint) {
      # Add manual point at time 0 with LFC = 0
      zero_point <- tibble(
        Species = .x$Species[1],
        Timepoint = "t0000",
        LFC = 0,
        Time_numeric = 0
      )
      bind_rows(zero_point, .x)
    } else {
      .x
    }
  }) %>%
  ungroup()

# Now use pho84_data_with_zero for plotting
pho84_plot <- ggplot(pho84_data_with_zero, aes(x = Time_numeric, y = LFC, color = Species)) +
  geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = custom_colors) +
  scale_x_continuous(
    breaks = unique(pho84_data_with_zero$Time_numeric),
    labels = unique(pho84_data_with_zero$Timepoint)
  ) +
  labs(
    title = "PHO84 Ortholog Expression",
    x = "Timepoint (mins)",
    y = "log2 Fold Change",
    color = "Species"
  ) +
  theme_minimal() +
  publication_theme +
  theme(
    axis.text.x = element_text(face = "bold", size = 13, angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold", size = 13),
    axis.title.y = element_text(face = "bold", size = 14),
    legend.position = "top"
  )

# Display the plot
pho84_plot
