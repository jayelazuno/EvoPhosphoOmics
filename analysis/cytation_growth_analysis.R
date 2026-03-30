## ================= Load Packages =======================
library(growthcurver)
library(tidyverse)
library(cowplot)
library(here)
library(gt)

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

## ================= Load & Process Cytation Data =======================
# Read your cytation data - UPDATE THIS PATH
cytation_raw <- read_csv(here("05.metadata","20230923-growth-PHO84-gfp-data", "20250709_log_growth_eo.csv"))

# Convert time format (H:MM:SS) to decimal hours
cytation_data <- cytation_raw %>%
  mutate(
    Time_hrs = sapply(Time, function(t) {
      parts <- as.numeric(unlist(strsplit(as.character(t), ":")))
      if (length(parts) == 3) {
        parts[1] + parts[2]/60 + parts[3]/3600
      } else {
        NA
      }
    })
  ) %>%
  select(-Time) %>%
  select(Time_hrs, everything())

# Reshape data and map wells to species/conditions
od <- cytation_data %>%
  pivot_longer(
    cols = -Time_hrs,
    names_to = "Well",
    values_to = "OD600"
  ) %>%
  mutate(
    # Extract row letter from well (A1 -> A, B12 -> B, etc.)
    Row = str_extract(Well, "^[A-H]"),
    
    # Map rows to species and phosphate conditions
    Species = case_when(
      Row %in% c("A", "E") ~ "S.cerevisiae",
      Row %in% c("B", "F") ~ "C.glabrata",
      Row %in% c("C", "G") ~ "K.lactis",
      Row %in% c("D", "H") ~ "C.albicans"
    ),
    
    Pi_mM = case_when(
      Row %in% c("A", "B", "C", "D") ~ 7.3,
      Row %in% c("E", "F", "G", "H") ~ 0
    )
  ) %>%
  filter(!is.na(OD600)) %>%
  # Extract well number and keep only wells 1-4 (4 replicates per condition)
  mutate(Well_Number = as.numeric(str_extract(Well, "\\d+"))) %>%
  filter(Well_Number %in% 1:4) %>%
  select(Time_hrs, Species, Pi_mM, Well, OD600)

# Set species order
species_order <- c("S.cerevisiae", "C.glabrata", "C.albicans", "K.lactis")
od$Species <- factor(od$Species, levels = species_order)

cat("✓ Data loaded:", nrow(od), "measurements\n")
cat("✓ Species:", paste(unique(od$Species), collapse = ", "), "\n")
cat("✓ Conditions:", paste(unique(od$Pi_mM), "mM", collapse = ", "), "\n\n")

## ================= Growth Curves =======================
growth_curve_plot <- od %>%
  mutate(Pi_Label = paste0(Pi_mM, " mM Pi")) %>%
  ggplot(aes(x = Time_hrs, y = OD600, color = Species)) +
  geom_point(size = 0.7, alpha = 0.6, position = position_jitter(width = 0.1)) +
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
cat("Fitting growth models...\n")

growth_params <- od %>%
  group_by(Species, Pi_mM) %>%
  group_modify(~ {
    time_hrs <- .x$Time_hrs
    OD_values <- .x$OD600
    
    # Check if there's adequate growth
    if (max(OD_values, na.rm = TRUE) < 0.2) {
      return(tibble(
        Carrying_Capacity = NA, 
        Growth_Rate = NA, 
        Doubling_Time = NA,
        Fit_Quality = NA, 
        Note = "Insufficient growth"
      ))
    }
    
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
        Carrying_Capacity = NA, 
        Growth_Rate = NA, 
        Doubling_Time = NA,
        Fit_Quality = NA, 
        Note = paste("Error:", e$message)
      )
    })
  }) %>%
  ungroup()

growth_params$Species <- factor(growth_params$Species, levels = species_order)

cat("✓ Models fitted\n\n")

growth_params

## ================= Parameter Plots =======================
# Carrying Capacity
k_plot <- ggplot(growth_params, aes(x = Species, y = Carrying_Capacity, fill = Species)) +
  geom_col(alpha = 0.8, show.legend = FALSE, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Carrying Capacity", y = expression(bold(OD[600])), x = "") +
  facet_wrap(~ paste0(Pi_mM, " mM Pi")) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Growth Rate
r_plot <- ggplot(growth_params, aes(x = Species, y = Growth_Rate, fill = Species)) +
  geom_col(alpha = 0.8, show.legend = FALSE, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Growth Rate", y = expression(bold("r (h"^{-1}*")")), x = "") +
  facet_wrap(~ paste0(Pi_mM, " mM Pi")) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

# Doubling Time
doubling_plot <- ggplot(growth_params, aes(x = Species, y = Doubling_Time, fill = Species)) +
  geom_col(alpha = 0.8, show.legend = FALSE, width = 0.6) +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Doubling Time", y = "Time (hours)", x = "") +
  facet_wrap(~ paste0(Pi_mM, " mM Pi")) +
  publication_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"))

k_plot
r_plot
doubling_plot

## ================= Combined Plot =======================
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

combined_plot

## ================= Table =======================
growth_table <- growth_params %>%
  mutate(across(c(Carrying_Capacity, Growth_Rate, Doubling_Time), ~round(., 3))) %>%
  select(Species, Pi_mM, Carrying_Capacity, Growth_Rate, Doubling_Time) %>%
  pivot_wider(
    names_from = Pi_mM,
    values_from = c(Carrying_Capacity, Growth_Rate, Doubling_Time),
    names_glue = "{.value}_{Pi_mM}mM"
  )

formatted_table <- growth_table %>%
  gt() %>%
  tab_header(
    title = "Growth Parameters Across Species and Phosphate Conditions",
    subtitle = "Carrying capacity (OD₆₀₀), growth rate (r), and doubling time (hours)"
  ) %>%
  fmt_number(decimals = 3) %>%
  cols_label(
    Species = "Species",
    Carrying_Capacity_0mM = "K (0 mM)",
    Carrying_Capacity_7.3mM = "K (7.3 mM)", 
    Growth_Rate_0mM = "r (0 mM)",
    Growth_Rate_7.3mM = "r (7.3 mM)",
    Doubling_Time_0mM = "Td (0 mM)", 
    Doubling_Time_7.3mM = "Td (7.3 mM)"
  ) %>%
  tab_spanner(label = "Carrying Capacity",
              columns = c(Carrying_Capacity_0mM, Carrying_Capacity_7.3mM)) %>%
  tab_spanner(label = "Growth Rate (h⁻¹)", 
              columns = c(Growth_Rate_0mM, Growth_Rate_7.3mM)) %>%
  tab_spanner(label = "Doubling Time (h)",
              columns = c(Doubling_Time_0mM, Doubling_Time_7.3mM)) %>%
  tab_style(style = cell_text(weight = "bold"),
            locations = list(cells_column_labels(), cells_column_spanners(), cells_title())) %>%
  tab_style(style = cell_text(style = "italic", weight = "bold"),
            locations = cells_body(columns = Species)) %>%
  tab_options(table.font.size = 12, heading.title.font.size = 16, 
              heading.subtitle.font.size = 14)

formatted_table

