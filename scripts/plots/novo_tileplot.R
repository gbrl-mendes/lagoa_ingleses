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
      mutate(words_count = str_count(`Curated ID`, "\\w+")) %>%
      filter(words_count == 2) %>%
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
    expeditions <- c(
                    "Nov_Dec/20",
                     "Nov/20",
                     "Dec/20",
                     "out/21",
                     "Nov/21")
  }
  
  # tabela so peixes
  {
    fish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids apenas os peixes a nivel de spp
      mutate(`Curated ID` = factor(`Curated ID`,
                                   levels = spp_fish)) %>% 
      mutate(`Order` = factor(`Order`,
                              levels = fish_ordens)) %>% 
      filter(`Curated ID` != "") %>%
      filter(!is.na(`Curated ID`))
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
      mutate('classificacao' = ifelse(`Curated ID` == "" | `Class` != "Actinopteri",
                                      "outros", `Curated ID`))
 
      all_ID_tbl <- all_ID_tbl %>% # ordem das ids da tabela all_ID_tbl
        mutate('classificacao' = factor(classificacao,
                                        levels = c("Acinocheirodon melanogramma",
                                                   "Astyanax fasciatus",
                                                   "Astyanax lacustris",
                                                   "Astyanax spp.",
                                                   "Brycon orthotaenia",
                                                   "Bryconamericus stramineus",
                                                   "Cichla spp.",
                                                   "Colossoma macropomum",
                                                   "Coptodon zillii",
                                                   "Eigenmannia virescens",
                                                   "Gymnotus carapo",
                                                   "Hemigrammus gracilis",
                                                   "Hemigrammus marginatus",
                                                   "Hoplias intermedius",
                                                   "Hoplias malabaricus",
                                                   "Hoplias spp.",
                                                   "Hypomasticus steindachneri",
                                                   "Leporellus vittatus",
                                                   "Leporinus piau",
                                                   "Leporinus reinhardti",
                                                   "Leporinus taeniatus",
                                                   "Megaleporinus elongatus",
                                                   "Megaleporinus garmani",
                                                   "Moenkhausia costae",
                                                   "Myleus micans",
                                                   "Orthospinus franciscensis",
                                                   "Pimelodus fur",
                                                   "Pimelodus maculatus",
                                                   "Pimelodus pohli",
                                                   "Pimelodus spp.",
                                                   "Planaltina myersi",
                                                   "Poecilia reticulata",
                                                   "Prochilodus costatus",
                                                   "Pseudoplatystoma corruscans",
                                                   "Pygocentrus piraya",
                                                   "Rhamdia quelen",
                                                   "Salmo salar",
                                                   "Serrasalmus brandtii",
                                                   "Siluriformes",
                                                   "Tilapia rendalli",
                                                   "Wertheimeria maculata",
                                                   "outros")))
  }
}


# Tile Plots ----

## Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5

  # Por campanha peixes
{
  tile_plot_campanhas <-
    fish_ID_tbl %>% 
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>% 
    mutate(Sample = factor(Sample, levels = samples)) %>% 
    mutate(Point = factor(Point)) %>% 
    mutate(File_name = factor(File_name)) %>% 
    mutate(`Order` = factor(`Order`)) %>% 
    ggplot(aes(y = `Curated ID`, 
               group = `Curated ID`, 
               x = Point, 
               fill = RRA)) +
    geom_tile() +
    # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
    #           colour = "black", size = 3) +
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

  ggsave(plot = tile_plot_campanhas, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/tileplots/tile_plot_campanhas_RRA.pdf",
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
    all_ID_tbl %>%
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>%
    mutate(Sample = factor(Sample, levels = samples)) %>%
    mutate(Point = factor(Point)) %>%
    mutate(File_name = factor(File_name)) %>%
    mutate(`Order` = factor(`Order`)) %>%
    mutate(`Class` = factor(`Class`)) %>%
    ggplot(aes(y = classificacao,
               group = classificacao,
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


  
# Alpha Plots ----
  
  # Paleta para o plot
  {color1 <- c("#D55E00", #(vermelho),
               "#0072B2", #(azul) 
               "#E69F00", #(laranja) 
               "#009E73", #(verde) 
               "#F0E442", #(amarelo) 
               "#56B4E9", #(azul claro) 
               "#CC79A7" #(rosa)
               )
  
    color2 <- c("#dc562a", # fonte: https://medialab.github.io/iwanthue/
              "#8e9fd9", 
              "#ca833f", 
              "#4d5a91", 
              "#cbcd9f", 
              "#2f2641", 
              "#d87d70", 
              "#313a25", 
              "#c94576", 
              "#4e7989", 
              "#8e3526", 
              "#d5acad", 
              "#591f2d", 
              "#550037", 
              "#8d6571",
              "#67d6b1", 
              "#6639c6", 
              "#6cda4e", 
              "#c24fd6", 
              "#c9db3e", 
              "#422478", 
              "#a7df83", 
              "#dc42ab", 
              "#519a40", 
              "#736dd8", 
              "#d9b647", 
              "#853a7d", 
              "#858438", 
              "#d48fce", 
              "#4e7b53", 
              "#db3c4e", 
              "#8acad8",
              "#4c8649",
              "#6e67b4",
              "#8f752d",
              "#b2523d",
              "#b04c7f",
              "#FF7F00", 
              "#008080", 
              "#FF1493", 
              "#00FA9A", 
              "#FF00FF" 
              )
  }
  
  
  # Ordem das ordens
  {
    ordem_ordens <- c("Outros", 
                     "Characiformes", 
                     "Cichliformes", 
                     "Cyprinodontiformes",
                     "Gymnotiformes", 
                     "Salmoniformes",
                     "Siluriformes")
  }
  

  # Tabela alpha_all_ID_tbl 
  {
    alpha_all_ID_tbl <- 
      raw_results_tbl %>% 
      mutate('classificacao' = ifelse(`Order` == "" | `Class` != "Actinopteri",
                                      "Outros", `Order`)) %>%
      mutate(classificacao = factor(classificacao,
                                    levels = ordem_ordens)) %>%
      mutate("Mês" = str_split_fixed(string = .$Sample, # criando uma nova coluna quebrando as infos da coluna Sample
                                     pattern = "_",
                                     n = 2)[,2]) %>%
      mutate(Mês = factor(Mês, levels = c("out21",
                                          "nov21",
                                          "dez20",
                                          "nov20"))) %>%
      mutate(Point = factor(Point, levels = c("L4",
                                              "L3",
                                              "L2",
                                              "L1"
                                              ))) %>%
      mutate(Expedition = factor(Expedition, levels = expeditions)) %>%
      group_by(Sample,
               Read_origin,
               Primer,
               `Curated ID`,
               classificacao,
               Year,
               Point,
               Expedition,
               Mês
      ) %>%
      summarize(`Num ASVs` = length(unique(`ASV (Sequence)`)),
                `Num OTUs` = length(unique(`OTU`)),
                `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
                `Ponto` = Point,
                `Sample` = Sample) %>%
      ungroup() %>%
      unique() 
  }
  
  # Alpha plot 
  {
    alpha_plot_id_sample <- 
      alpha_all_ID_tbl %>% 
      mutate(classificacao = factor(classificacao)) %>%
            ggplot(aes(x = Point, 
                 y = `ID Abundance on sample (%)`,
                 order = classificacao,
                 group = Sample, 
                 fill = as.factor(classificacao))) +
      ggtitle(label = "Proporção de ordens identificadas") + 
      xlab("Ponto") + 
      ylab("Proporção %") + 
      guides(fill = guide_legend("Espécie")) + ## alterar o titulo da legenda 
      geom_bar(stat = "identity", position = "stack") +
      coord_flip() +
      facet_grid(cols = vars(Expedition)) +
      scale_fill_manual(values = color1)
      
    
    
    ## Plotando
    ggsave(plot = alpha_plot_id_sample, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/novembro/barplots/alpha_plot_id_sample.pdf",
           device = "pdf", units = "cm", height = 15, width = 28, dpi = 600)
  }