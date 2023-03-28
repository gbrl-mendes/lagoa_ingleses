---
title: "Novo Tile Plot"
author: "Gabriel Mendes"
date: "22/03/2023"
---

# Carregando bibliotecas ----
{
  library(adespatial)
  library(base)
  library(Biostrings)
  library(DECIPHER)
  library(factoextra)
  library(future)
  library(ggforce)
  library(ggh4x)
  library(ggord)
  library(ggplot2)
  library(ggpubr)
  library(ggvegan)
  library(Matrix)
  library(phyloseq)
  library(ShortRead)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(vegan)
  library(dplyr)
}

# Caminhos ----
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figuras")
  tbl_path <- paste0(prjct_path,"/tabelas/raw/run_2_4_5")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

# Obtencao dos dados ----
{
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v2023.csv"), sep = ",", check.names = FALSE) %>% tibble()
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Oreochromis niloticus")] <- "Tilapia rendalli" # Oreochromis niloticus é Tilapia
}

# Tile Plots ----

# Criacao das tabelas all_ids_tbl e few_ids_tbl
{
  # Criacao da lista com os possiveis nomes atribuidos as ASVs
  {
    raw_results_tbl %>% colnames()
    raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  }
  
  # Agrupamento da ASVs que possuem os mesmos atributos abaixo
  {
    grouped_by_ID_BLASTid <- raw_results_tbl %>%
      select(c(
        "Sample",
        "Expedition",
        "Year",
        "Sample.Name",
        "File_name",
        "Curated ID",
        "Relative abundance on sample",
        "Point",
        "Class",
        "Order"
      )) %>% 
      group_by(Sample, `Curated ID`, Expedition, Point, Sample.Name, File_name, ) %>% 
      summarize(`Sample` = unique(Sample),
                `Curated ID` = unique(`Curated ID`),
                `Class`,
                `Order`,
                # `Class` = unique(`Class`),
                # `Order` = unique(`Order`),
                `Expedition` = unique(Expedition),
                `Year` = unique(Year),
                `Point` = unique(Point),
                `Sample.Name` = unique(Sample.Name),
                `File_name` = unique(File_name),
                `RRA` = sum(`Relative abundance on sample`)) %>%
      ungroup()
  }
  
  # Organizar as especies
  {
    grouped_by_ID_BLASTid$`Curated ID` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
    grouped_by_ID_BLASTid$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  }
  
  # Organizar ordem
  {
    fish_ordens <- c("Characiformes", # ordem das ordens de peixes 
                     "Siluriformes", 
                     "Cichliformes", 
                     "Gymnotiformes", 
                     "Cyprinodontiformes", 
                     "Salmoniformes")
  }  
  
  {
    nfish_classes <- c("Aves", # ordem das classes de nao peixes
                       "Actinobacteria",
                       "Mammalia")
  }
  
  {
    nfish_ordens <- c("Actinomycetales", # ordem das ordens de nao peixes 
                      "Artiodactyla", 
                      "Carnivora",
                      "Didelphimorphia",
                      "Lagomorpha",
                      "Primates",
                      "Rodentia",
                      "Suliformes")
  }
  
  # Definir quais serao as especies e ordenar
  
  {
    spp_fish <- raw_results_tbl %>%  # especies de peixes
      filter(Class == "Actinopteri") %>%
      # mutate(words_count = str_count(`Curated ID`, "\\w+")) %>% 
      # filter(words_count == 2) %>% 
      # filter(`Curated ID` != "") %>% 
      select(`Curated ID`) %>% 
      unique() %>% 
      arrange(`Curated ID`) %>% 
      pull()
  }
  
  {
    no_fish <- raw_results_tbl %>%  # outras especies
      filter(Class != "Actinopteri") %>% 
      filter(`Curated ID` != "") %>% 
      select(`Curated ID`) %>% 
      unique() %>% 
      arrange(`Curated ID`) %>% 
      pull()
  }  
  
  # Definir quais serao as amostras e ordenar
  {
    samples <- c("L1_nov21",
                 "L1_out21",
                 "L2_nov21",
                 "L2_out21",
                 "L3_nov21",
                 "L3_out21",
                 "L4_nov21",
                 "L4_out21",
                 "L1-neo-mi",
                 "L2 Dez/20",
                 "L2 Nov/20")
  }
  
  # Definir quais serao as campanhas e ordenar
  {
    expeditions <- c("Nov_Dec/20",
                       "Nov/20",
                       "Dec/20",
                       "out/21",
                       "Nov/21")
  }
  
  # tabela so peixes
  {
    fish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids apenas os peixes a nivel de spp
      mutate(`Curated ID` = factor(`Curated ID`,
                                     levels = rev(spp_fish))) %>% 
      mutate(`Order` = factor(`Order`,
                                       levels = fish_ordens))
  }
  
  # tabela nao peixes
  {
    nfish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids nao peixes a nivel de spp
      mutate(`Curated ID` = factor(`Curated ID`,
                                     levels = rev(no_fish))) %>% 
      mutate(`Order` = factor(`Order`,
                                       levels = nfish_ordens))
  }
  
  # tabela peixes & outros
  {
    all_ID_tbl <- grouped_by_ID_BLASTid %>% # ids nao peixes a nivel de spp
      # mutate(`Curated ID` = factor(`Curated ID`,
      #                                levels = rev(no_fish))) %>% 
      # mutate(`Order` = factor(`Order`,
      #                                  levels = nfish_ordens)) %>% 
      mutate('classificacao' = ifelse(`Curated ID` == "" | `Class` != "Actinopteri", "outros", `Curated ID`))
  }
}

  # ordem das ids da tabela all_ID_tbl DESCOBRIR COMO QUE COLOCA O "outros" NO FINAL DA ORDEM DAS ESPECIES
{
  all_ids <- all_ID_tbl %>%
    select(classificacao) %>%
    mutate(classificacao = ifelse(classificacao == "outros", "zzz_outros", classificacao)) %>%
    unique() %>%
    arrange(classificacao) %>%
    slice(-which(classificacao == "zzz_outros")) %>% 
    pull()

}

## Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5

  # Por campanha peixes
{
  tile_plot_campanhas <-
    fish_ID_tbl %>% 
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>% 
    mutate(Sample = factor(Sample, levels = samples)) %>% 
    mutate(Point = factor(Point)) %>% 
    mutate(File_name = factor(File_name)) %>% 
    filter(!is.na(`Curated ID`)) %>%
    mutate(`Order` = factor(`Order`)) %>% 
    ggplot(aes(y = `Curated ID`, 
               group = `Curated ID`, 
               x = Point, 
               fill = RRA)) +
    geom_tile() +
    geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
              colour = "black", size = 3) +
    facet_grid(cols = vars(Expedition), 
               rows = vars(`Order`),
               space = "free_y", drop = TRUE, 
               scales = "free_y") +
    labs(fill ='Abundância \nrelativa (%)',
         x = "Pontos",
         y= "Espécies") +
    ggtitle(label = "Espécies identificadas") +
    scale_fill_continuous(
      trans = "log10",
      breaks = c(0.001, 0.01, 0.1, 1, 10, 75),
      type = "viridis") +
    theme(text=element_text(size = 13)) +
    theme(strip.text.y = element_text(angle=0))
}

  ggsave(plot = tile_plot_campanhas, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/tileplots/tile_plot_campanhas.pdf",
         device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)

# Por campanha nao-peixes
{
  tile_plot_campanhas <-
    nfish_ID_tbl %>%
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>%
    mutate(Sample = factor(Sample, levels = samples)) %>%
    mutate(Point = factor(Point)) %>%
    mutate(File_name = factor(File_name)) %>%
    filter(!is.na(`Curated ID`)) %>%
    mutate(`Order` = factor(`Order`)) %>%
    mutate(`Class` = factor(`Class`)) %>%
    ggplot(aes(y = `Curated ID`,
               group = `Curated ID`,
               x = Point,
               fill = RRA)) +
    geom_tile() +
    geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
              colour = "black", size = 3) +
    facet_grid(cols = vars(Expedition),
               rows = vars(`Class`, `Order`),
               space = "free_y", drop = TRUE,
               scales = "free_y") +
    labs(fill ='Abundância \nrelativa (%)',
         x = "Pontos",
         y= "Espécies") +
    ggtitle(label = "Espécies identificadas") +
    scale_fill_continuous(
      trans = "log10",
      breaks = c(0.001, 0.01, 0.1, 1, 10, 75),
      type = "viridis") +
    theme(text=element_text(size = 13)) +
    theme(strip.text.y = element_text(angle=0))
  }
  
  # Por campanha tudo
{
  tile_plot_campanhas <-
    al %>%
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>%
    mutate(Sample = factor(Sample, levels = samples)) %>%
    mutate(Point = factor(Point)) %>%
    mutate(File_name = factor(File_name)) %>%
    mutate(`Order` = factor(`Order`)) %>%
    mutate(`Class` = factor(`Class`)) %>%
    ggplot(aes(y = `Curated ID`,
               group = `Curated ID`,
               x = Point,
               fill = RRA)) +
    geom_tile() +
    geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
              colour = "black", size = 3) +
    facet_grid(cols = vars(Expedition),
               rows = vars(`Class`, `Order`),
               space = "free_y", drop = TRUE,
               scales = "free_y") +
    labs(fill ='Abundância \nrelativa (%)',
         x = "Pontos",
         y= "Espécies") +
    ggtitle(label = "Espécies identificadas") +
    scale_fill_continuous(
      trans = "log10",
      breaks = c(0.001, 0.01, 0.1, 1, 10, 75),
      type = "viridis") +
    theme(text=element_text(size = 13)) +
    theme(strip.text.y = element_text(angle=0))
  }
  
  ggsave(plot = tile_plot_campanhas, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/tileplots/tile_plot_campanhas.pdf",
         device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)

