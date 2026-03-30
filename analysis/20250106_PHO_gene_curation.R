---
title: "PHO regulon literature curation and support scoring"
author: "Joshua Ayelazuno"
date: "2025-01-06"
output: html_document
---
  
library(tidyverse)
library(readxl)
library(here)
library(ggvenn)

pho_file <- here("06.results", "Scer", "20251229_PHO_Regulon.xlsx")

# sheets 2–6 
paper_sheets <- tibble(
  sheet = c(
    "Gurvich_et_al_2017",
    "He_et_al_2017",
    "Zhou_OShea_2011",
    "Springer_et_al_2003",
    "Ogawa_et_al_2000"
  ),
  weight = c(
    1.0,  # Gurvich 2017
    0.5,  # He 2017 (derived from Zhou)
    0.5,  # Zhou & O'Shea 2011
    1.0,  # Springer 2003
    1.0   # Ogawa 2000
  )
)

# ---- read sheets into a nested tibble ----
pho_by_paper <- paper_sheets %>%
  mutate(data = map(sheet, ~
                      read_excel(
                        pho_file,
                        sheet = .x,
                        col_types = "text"
                      ) %>%
                      transmute(
                        gene_id   = str_trim(as.character(gene_id)),
                        gene_name = str_trim(as.character(gene_name)),
                        paper     = .x
                      ) %>%
                      # drop totally empty rows
                      filter(!(is.na(gene_id) & is.na(gene_name))) %>%
                      filter(!(gene_id == "" & gene_name == ""))
  ))

# ---- LONG table (needed for scoring + Venn) ----
pho_long <- pho_by_paper %>%
  select(sheet, weight, data) %>%
  unnest(data) %>%
  mutate(
    gene = if_else(is.na(gene_name) | gene_name == "", gene_id, gene_name),
    gene = str_trim(as.character(gene))
  ) %>%
  select(gene_id, gene_name, gene, paper = sheet, weight)

# ---- named list (for Venn) ----
pho_list <- pho_long %>%
  distinct(paper, gene) %>%
  group_by(paper) %>%
  summarise(genes = list(sort(unique(gene))), .groups = "drop") %>%
  deframe()

# ---- support scoring (unique gene x paper) ----
pho_support <- pho_long %>%
  distinct(gene, paper, .keep_all = TRUE) %>%
  group_by(gene) %>%
  summarise(
    n_papers = n(),
    support_score = sum(weight),
    supporting_papers = paste(paper, collapse = "; "),
    .groups = "drop"
  )

# exclude genes supported by only ONE paper
pho_supported_final <- pho_support %>%
  filter(n_papers > 1) %>%
  arrange(desc(support_score), desc(n_papers), gene)

# attach gene_id / gene_name back
pho_final <- pho_supported_final %>%
  left_join(
    pho_long %>%
      select(gene, gene_id, gene_name) %>%
      distinct(),
    by = "gene"
  ) %>%
  select(
    gene_id,
    gene_name,
    n_papers,
    support_score,
    supporting_papers
  ) %>%
  arrange(desc(support_score), desc(n_papers), gene_name %||% gene_id)

write_csv(
  pho_final,
  here("06.results", "Scer", "20250106_curated_PHO_genes.csv")
)

#Venn diagram 
ggvenn(
  pho_list,
  fill_color = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"),
  stroke_size = 0.5,
  set_name_size = 5
) +
  theme_void() +
  ggtitle("Literature support for PHO regulon genes")


