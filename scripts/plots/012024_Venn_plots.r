---
title: "Venn plot Lagoa dos Ingleses "
author: "Gabriel Mendes"
date: "01/2024"
---

# Carregando bibliotecas ----
{
  library(Biostrings)
  library(DECIPHER)
  library(factoextra)
  library(future)
  library(ggh4x)
  library(ggord)
  library(ggpubr)
  library(ggvegan)
  library(Matrix)
  library(phyloseq)
  library(ShortRead)
  library(stringr)
  library(vegan)
  library(tidyverse)
  library(dplyr)
}

# Obtencao de dados ----

# All the data was generated by the 012024_post_plot.r code 

# Filtrando ----

resume_venn <- fish_ID_tbl %>%
  filter(expedition %in% c("Novembro 2021")) %>%
  select(`Curated ID`, `Curated genus`, `Family (BLASTn)`, `Curated Order (BLASTn)`, filter) %>%
  unique()

# Especies em comum MCE e Sterivex
spp_comm <- resume_venn %>% 
  group_by(`Curated ID`) %>% # Trocar para spp, genero, familia e ordem
  summarise(n_filtros = n_distinct(filter)) %>%
  filter(n_filtros >= 2) %>%
  select(`Curated ID`) # Trocar para spp, genero, familia e ordem
View(spp_comm)
  
# Especies exclusivas
spp_diff <- resume_venn %>%
  group_by(`Curated ID`) %>% # Trocar para spp, genero, familia e ordem
  summarise(MCE_count = sum(as.integer(filter == "MCE")),
            STX_count = sum(as.integer(filter == "Sterivex")),
            MCE_only = MCE_count >= 1 & STX_count == 0,
            STX_only = MCE_count == 0 & STX_count >= 1) %>% 
  ungroup()
View(spp_diff)
