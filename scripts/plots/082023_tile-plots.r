---
title: "Tile plots Lagoa dos Ingleses"
author: "Gabriel Mendes"
date: "22/03/2023"
---

## Carregando bibliotecas ----
{
  library(adespatial)
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
}

## Caminhos ----
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figuras")
  tbl_path <- paste0(prjct_path,"/tabelas/raw/run_2_4_5")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

## Obtencao dos dados ----

  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v062023.csv"), sep = ",", check.names = FALSE) %>% tibble()
  
## Alteracoes pos-pipeline ----
{ 
  # Correcao de nomes 
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Oreochromis niloticus")] <- "Coptodon sp."#"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Coptodon zillii")] <- "Coptodon sp."#"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Tilapia rendalli")] <- "Coptodon sp." #"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Cichla spp.")] <- "Coptodon sp." #"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
    
  # Alteracoes pos construcao das arvores
  
  # Remocao das ASVs de Teleo. Como havia cabecalhos (ASV header) com nome igual 
  # de ASVs que eram pra ficar, optei por retirar pela sequencia 
  teleo_ASVs_seq <- c(
    # ASV_384_175bp|Rhamdia_quelen
    "GTGAAGGACTAGCAGTAAGTAAAATTGGTACAACCAAAAACGTCAGGTCGAGGTGTAGCGAATGAAGTGGAAAGAAATGGGCTACATTTTCTATTCAATAGAATAACACGAATGACACCCTGAAACATATGTCTGAAGGTGGATTTAGTAGTAAAAAGCAAATAGTGTGTCCTTTTGAATTAGGCTCTGAGACGC",
    #ASV_15_167bp|Moenkhausia_costae
    "GTGAAGGAATAAAAGTAAGCAAAATGGGTCTCACCCAAGAACGTCAGGTCGAGGTGTAGCTTACGGAGTGGAAAGAAATGGGCTACATTTTCTATCCCAGAAAACTACGAACGGCACTGTGAAATCAGTGCCCGAAGGAGGATTTAGCAGTAAAAATCAAATAGAGCGTGATTTTGAAGCAGGCTCTGAGAAGC",
    #ASV_326_175bp|Hoplias_spp.#
    "GTGAAGGTACTACAGTAAGCAAAATGGGTAAAACCCAAAACGTCAGGTCGAGGTGTAGCTTATGGGGTGGAAAGAAATGGGCTACATTTCCTACATAAAAGGATATTACGAACGACCTCATGAAACTGAAGTCTTAAGGTGGATTTAGCAGTAAATATAAGTAGAGCATTCTTTTGAAGCCGGCTCTGAAACGC",
    #ASV_201_172bp|Hoplias_malabaricus
    "GTGAAGGGCCTACAGTAAGCAAAATGGGCAAGACCCAGAACGTCAGGTCGAGGTGTAGCTTATGGAGTGGAAAGAAATGGGCTACATTTCCTTCACTAAAGGATATCACGAACGACTTCATGAAACTGAAGTCCCAAGGTGGATTTAGCAGTAAACATTAATAGAGTATTATTTTGAAGCCGGCTCTGAAACGC",
    #ASV_146_170bp|Coptodon_sp    
    "GTGAAGGGTCTACAGTAAGCAAAACTAGTACAACTCAAAACGCCAGGTCGAGGTGTAGCATATGAGAGGGGAAGAAATGGGCTACATTCCCTACCGCAGGGAACACGAACAATGTAATGAAACGTACATTAAAAGGAGGATTTAGCAGTAAGCAGAAAATAGAGCGTTCCGCTGAAATCGGCCCTGAAGCGC")
  
  raw_results_tbl <- raw_results_tbl %>% 
    filter(!`ASV (Sequence)` %in% teleo_ASVs_seq)
  
  # Nova correcao de nomes
  # Ver o documento Verificacao das identificacoes do DADA2
  # Substituindo o nome de todas ASVs 
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Hoplias brasiliensis")] <- "Hoplias sp." # Ver artigo do Heron 2023
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Hoplias intermedius")] <- "Hoplias sp." # Ver artigo do Heron 2023
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Astyanax fasciatus")] <- "Psalidodon fasciatus" # A. fasciatus. por Psalidodon fasciatus
  
  # Substituindo ASVs especificas
  raw_results_tbl <- # Orthospinus por Moenkhausia
    raw_results_tbl %>% mutate(`Curated ID` =
                                 ifelse(`ASV header` == ">ASV_287_175bp" &
                                          `Curated ID` == "Orthospinus franciscensis",
                                        "Moenkhausia costae", `Curated ID`))
  
  raw_results_tbl <- # Astyanax spp. por A. lacustris
    raw_results_tbl %>% mutate(`Curated ID` = 
             ifelse((`ASV header` == ">ASV_38_168bp" | `ASV header` == ">ASV_34_168bp") &
                      str_detect(`Curated ID`, "Astyanax spp."), 
                    "Astyanax lacustris", `Curated ID`))
  
  raw_results_tbl <- # Astyanax spp. por A. fasciatus
    raw_results_tbl %>% mutate(`Curated ID` = 
             ifelse((`ASV header` == ">ASV_658_168bp" | `ASV header` == ">ASV_657_168bp") &
                      str_detect(`Curated ID`, "Astyanax spp."), 
                    "Astyanax fasciatus", `Curated ID`))  
  
  raw_results_tbl <- # Poecilia reticulata por A. melanogramma
    raw_results_tbl %>% mutate(`Curated ID` =
                                 ifelse(`ASV header` == ">ASV_530_170bp" &
                                          `Curated ID` == "Poecilia reticulata",
                                        "Acinocheirodon melanogramma", `Curated ID`))
  }

## Tabelas ----
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
        "New_name",
        "Class",
        "Order",
        "Abundance",
        "ASV header",
        "ASV (Sequence)"
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
                `New_name` = unique(New_name),
                `Sample.Name` = unique(Sample.Name),
                `File_name` = unique(File_name),
                `Abundance`,
                `RRA` = sum(`Relative abundance on sample`),
                `ASV header`,
                `ASV (Sequence)`
                ) %>%
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
      filter(words_count >= 2) %>%
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
    samples <- c("L1_out21",
                 "L2_out21",
                 "L3_out21",
                 "L4_out21",
                 "L1_nov21",
                 "L2_nov21",
                 "L3_nov21",
                 "L4_nov21",
                 "LI1-neo-mi",
                 "L2_nov20",
                 "L2_dez20")
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
                              levels = fish_ordens)) %>%
      mutate('Nível' = ifelse(`Year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                                      "Cheio", "Vazio")) %>%
      filter(`Curated ID` != "") %>%
      filter(!is.na(`Curated ID`)) %>% 
      arrange(`Curated ID`)
  }
  
  # tabela nao peixes
  {
    nfish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids nao peixes a nivel de spp
      mutate(`Curated ID` = factor(`Curated ID`,
                                   levels = (no_fish))) %>% 
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
                                                   "Tilapia rendalli/Coptodon zillii",
                                                   "Wertheimeria maculata",
                                                   "outros")))
  }
}

## Histograms ----

# Histogramas para visualizar as frequencias das ASVs e observar qual o 
# melhor threshold para eliminar ASVs espurias

# Peixes

library("plotly")
options(scipen = 100)

  hst_peixes  <- fish_ID_tbl %>%
    filter(RRA > 0.3) %>%
    filter(Expedition %in% c(
      # "Nov_Dec/20",
      # "Nov/20",
      # "Dec/20"
      "out/21",
      "Nov/21"
    )) %>%
    ggplot(aes(x = RRA, 
               fill = `Curated ID`,
               label = `Curated ID`)) +
    geom_histogram(binwidth = 0.25, color = "black") +
    labs(x = "Abundância", y = "Frequência") +
    theme_minimal() +
    scale_x_log10(breaks = c(0,0.001,0.01,0.1,1,5,10,50,100)) + 
    scale_y_log10() +
    # # facet_grid(rows = vars(Year)) +
    facet_grid(rows =  vars(Point))

  hst_peixes %>% ggplotly()
  
  
## nao-peixes
  
 hst_npeixes <- nfish_ID_tbl %>% 
    ggplot(aes(x = RRA, 
               fill = `Curated ID`,
               label = `Curated ID`)) +
    geom_histogram(binwidth = 0.25, color = "black") +
    labs(x = "Abundância", y = "Frequência") +
    theme_minimal() +
    scale_x_log10(breaks = c(0,0.001,0.01,0.1,1,5,10,50,100))+
    scale_y_log10()+
    facet_grid(rows = vars(Year))
  
 hst_npeixes %>% ggplotly()

## Tile Plots ----

# Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5

  ## Por campanha peixes
  
# 2020
{
  tile_plot_20 <-
    fish_ID_tbl %>% 
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>% 
    mutate(Sample = factor(Sample, levels = samples)) %>% 
    mutate(Point = factor(Point)) %>% 
    mutate(File_name = factor(File_name)) %>% 
    mutate(`Order` = factor(`Order`)) %>%
    filter(Expedition %in% c("Nov_Dec/20",
                             "Nov/20", # deixando apenas as amostras de 2020
                             "Dec/20"
    )) %>% 
    droplevels() %>% # remover os níveis não utilizados em cada fator após o filtro
    ggplot(aes(y = `Curated ID`, 
               group = `Curated ID`, 
               x = Point, 
               fill = RRA)) +
    geom_tile() +
    geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
              colour = "black", size = 5) +
    facet_grid(cols = vars(Expedition), 
               rows = vars(`Order`),
               space = "free",
               scales = "free",
               drop = TRUE) +
    labs(fill ='Abundância \nrelativa (%)',
         x = "Pontos",
         y= "Espécies") +
    ggtitle(label = "Nível normal - espécies identificadas (2020)") +
    scale_fill_continuous(trans = "log10", 
                          breaks = c(0.01, 1, 60),
                          # breaks = c( 0.001, 1, 60),
                          type = "viridis") +
    # theme(text=element_text(size = 13)) +
    theme(axis.text.x = element_text(size = 15, face = "bold"), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 19),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          strip.text = element_text(size = 13, face = "bold")) +
    theme(strip.text.y = element_text(angle=0))
}
ggsave(plot = tile_plot_20, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/tileplots/tile_plot_20.pdf",
       device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)

## 2021
{
  tile_plot_21 <-
    fish_ID_tbl %>%
    # filter(RRA >= 0.01) %>% 
    mutate(Expedition = factor(Expedition, levels = expeditions)) %>% 
    mutate(Sample = factor(Sample, levels = samples)) %>% 
    mutate(Point = factor(Point)) %>% 
    mutate(File_name = factor(File_name)) %>% 
    mutate(`Order` = factor(`Order`)) %>%
    # filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
    #                          "Nov/21"
    # )) %>% View()
    ggplot(aes(y = `Curated ID`, 
               group = `Curated ID`, 
               x = Point, 
               fill = RRA)) +
    geom_tile() +
    geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
              colour = "black", size = 4) +
    facet_grid(cols = vars(Expedition), 
               rows = vars(`Order`),
               space = "free", 
               scales = "free",
               drop = TRUE) +
    labs(fill ='Abundância \nrelativa (%)',
         x = "Pontos",
         y= "Espécies") +
    ggtitle(label = "Nível baixo - espécies identificadas (2021)") +
    scale_fill_continuous(
      trans = "log10",
      breaks = c( 0.001, 1, 75),
      type = "viridis") +
    theme(axis.text.x = element_text(size = 15, face = "bold"), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 15, face = "bold"),
          axis.title.x = element_text(size = 16, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"),
          plot.title = element_text(size = 19),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13),
          strip.text = element_text(size = 13, face = "bold")) +
    theme(strip.text.y = element_text(angle=0))
}
  ggsave(plot = tile_plot_21, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/tileplots/tile_plot_21.pdf",
         device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)

## 2020 (cheio) versus 2021 (vazio)
{
  {
    tile_plot_all <-
      fish_ID_tbl %>% 
      filter(RRA >= 0.01) %>%
      mutate(Expedition = factor(Expedition, levels = expeditions)) %>%
      mutate(Sample = factor(Sample, levels = samples)) %>%
      mutate(Point = factor(Point)) %>%
      mutate(File_name = factor(File_name)) %>%
      mutate(`Order` = factor(`Order`)) %>%
      group_by(Expedition, `Curated ID`, Point, New_name, Order, Nível) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(New_name %in% c("Ponte", "Fundacao")) %>%
      # filter(Sample %in% c("LI1-neo-mi", "L2_nov20", "L2_dez20", "L2_nov21", "L3_nov21", "L4_nov21")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`, 
                 x = 
                   # Point,
                   `New_name`,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      #           colour = "black", size = 4) +
      facet_grid(cols = vars(`Nível`), 
                 rows = vars(`Order`),
                 space = "free_y", 
                 scales = "free_y",
                 drop = TRUE) +
      # scale_x_discrete(limits = c("Fundacao", "Ponte"), expand = c(0, 0)) + # exibindo apenas o 2 pontos que existem em comum nas duas campanhas
      labs(fill ='Abundância \nrelativa (%)',
           x = "Pontos",
           y= "Espécies") +
      ggtitle(label = "Variação de espécies: reservatório cheio vs. vazio") +
      scale_fill_continuous(
        trans = "log10",
        breaks = c( 0.001, 1, 75),
        type = "viridis") +
      theme(axis.text = element_text(size = 15, face = "bold"), ## definindo o tamanho de cada texto
            axis.title = element_text(size = 16, face = "bold"),
            plot.title = element_text(size = 19, face = "bold", hjust = 0.5),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 13),
            strip.text = element_text(size = 13, face = "bold")) +
      theme(strip.text.y = element_text(angle=0))
    
    }  
  
  ggsave(plot = tile_plot_all, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/tile_plot/",
                          "tile_cheio_vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 20,
         width = 30,
         dpi = 600)
  
  # Spp com maior frequencia em cada ponto
  {
    top_seqs <-
      fish_ID_tbl %>%
      filter(Expedition %in% c(
        # "Nov_Dec/20",
        #                        "Nov/20", # deixando apenas as amostras de 2020
        #                        "Dec/20"
        "out/21", # deixando apenas as amostras de 2021
        "Nov/21"
      )) %>% 
      group_by(Point) %>%
      unique() %>% 
      top_n(3, RRA)
  }
  
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
  } 
  
  ggsave(plot = tile_plot_campanhas, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/tileplots/tile_plot_campanhas.pdf",
         device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)

## vazio versus cheio
  
  {
    cheio_vs_vazio <-
      fish_ID_tbl %>%
      filter(RRA >= 0.01) %>%
      mutate(Expedition = factor(Expedition, levels = expeditions)) %>% 
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(Point = factor(Point)) %>% 
      mutate(File_name = factor(File_name)) %>% 
      mutate(`Order` = factor(`Order`)) %>%
      group_by(Expedition, `Curated ID`, Point, New_name, Order, Nível) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(New_name %in% c("Ponte", "Fundacao")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`,
                 x = New_name, 
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
                # stat = "identity",
                # colour = "black", size = 4) +
      facet_grid(cols = vars(Nível),
                 labeller = labeller(Nível = c("Cheio" = "Full", "Vazio" = "Empty")),
                 rows = vars(`Order`),
                 space = "free", 
                 scales = "free",
                 drop = TRUE) +
      labs(fill ='Relative \nabundance (%)',
           x = "Points",
           y= "Species") +
      ggtitle(label = "Espécies detectadas: Lagoa cheia e vazia") +
      scale_fill_continuous(
        trans = "log10",
        breaks = c( 0.001, 1, 75),
        type = "viridis") +
      theme(plot.title = element_text(size = 25, face = "bold", hjust=0.7),
            axis.text.x = element_text(size = 13, face = "bold", angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 13, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 13),
            strip.text = element_text(size = 13, face = "bold")) +
      theme(strip.text.y = element_text(angle=0)) +
      theme(legend.position = "bottom") 
  }
  ggsave(plot = cheio_vs_vazio, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/tile_plot/",
                          "tile_plot_cheio_vs_vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf", units = "cm", height = 20, width = 20, dpi = 600)
  
 ## Alpha Plots (versao qualificacao) ----
  
  # Paleta para o plot
  {color1 <- c("#0072B2", #(azul) 
               "#E69F00", #(laranja)  
               "#CC79A7", #(rosa)
               "#F0E442", #(amarelo) 
               "#56B4E9", #(azul claro) 
               "#009E73", #(verde)               
               "#D55E00" #(vermelho)
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
    ordem_ordens <- c("Characiformes", 
                     "Cichliformes", 
                     "Cyprinodontiformes",
                     "Gymnotiformes", 
                     "Salmoniformes",
                     "Siluriformes",
                     "Outros")
  }
  

  # Tabela alpha_all_ID_tbl 
  {
    alpha_all_ID_tbl <- 
      raw_results_tbl %>% 
      mutate('classificacao' = ifelse(`Order` == "" | `Class` != "Actinopteri", ## com essa linha a gente inclui tudo que nao e' peixe
                                      "Outros", `Order`)) %>%                   ## e tudo que nao ta ao nivel de genero colocando na 
      mutate(classificacao = factor(classificacao,                              ## classificacao outros
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
               Class,
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
  ## esse plot inclui as ASVs que nao foram identificadas ate o nivel de Ordem
  
  #2021
  {
    alpha_plot_ordem <- ## barplot das ordens encontradas
      alpha_all_ID_tbl %>% 
      mutate(classificacao = factor(classificacao)) %>%
      filter(Expedition %in% c(
                               # "Nov_Dec/20",
                               # "Nov/20",
                               # "Dec/20"
                               "out/21", # deixando apenas as amostras de 2021
                               "Nov/21"
                               )) %>%
            ggplot(aes(x = Point, 
                 y = `ID Abundance on sample (%)`,
                 order = classificacao,
                 group = Sample, 
                 fill = as.factor(classificacao))) +
      xlab("Ponto") + 
      ylab("Proporção %") + 
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(cols = vars(Expedition)) +
      # facet_wrap(~ Expedition, ncol = 1) +
      coord_flip() + #virando 90o o grafico
      ggtitle(label = "Proporção de ordens identificadas (2021)") +  ## alterando o titulo do grafico
      guides(fill = guide_legend("Ordens")) + ## alterar o titulo da legenda 
      theme(axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 23),
            axis.title.y = element_text(size = 23),
            plot.title = element_text(size = 22),
            legend.text = element_text(size = 19),
            legend.title = element_text(size = 19),
            strip.text = element_text(size = 19),
            legend.position = "bottom",
            legend.background = element_rect(fill = "lightgray", color = NA)) +
      scale_fill_manual(values = color1)
      
    ## Plotando
    ggsave(plot = alpha_plot_ordem, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/alpha_plot_ordem_2021.pdf",
           device = "pdf", units = "cm", height = 20, width = 27.5, dpi = 600)
  }
  #2020
  {alpha_plot_ordem_2020 <- ## barplot das ordens encontradas
      alpha_all_ID_tbl %>% 
      mutate(classificacao = factor(classificacao)) %>%
      filter(Expedition %in% c(
        "Nov_Dec/20",
        "Nov/20",
        "Dec/20"
        # "out/21", # deixando apenas as amostras de 2021
        # "Nov/21"
      )) %>%
      ggplot(aes(x = Point, 
                 y = `ID Abundance on sample (%)`,
                 order = classificacao,
                 group = Sample, 
                 fill = as.factor(classificacao))) +
      xlab("Ponto") + 
      ylab("Proporção %") + 
      geom_bar(stat = "identity", position = "stack") +
      facet_grid(cols = vars(Expedition)) +
      # facet_wrap(~ Expedition, ncol = 1) +
      coord_flip() + #virando 90o o grafico
      ggtitle(label = "Proporção de ordens identificadas (2020)") +  ## alterando o titulo do grafico
      guides(fill = guide_legend("Ordens")) + ## alterar o titulo da legenda 
      theme(axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 23),
            axis.title.y = element_text(size = 23),
            plot.title = element_text(size = 22),
            legend.text = element_text(size = 19),
            legend.title = element_text(size = 19),
            strip.text = element_text(size = 19),
            legend.position = "bottom",
            legend.background = element_rect(fill = "lightgray", color = NA)) +
      scale_fill_manual(values = color1)
    
    ggsave(plot = alpha_plot_ordem_2020, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/alpha_plot_ordem_2020.pdf",
           device = "pdf", units = "cm", height = 12, width = 35, dpi = 600)
  }
    
  
# Numero de ASVs
## esse plot inclui apenas as spp de peixes ao nivel no minimo de genero

  ## 2021
  {
  asvs_ponto_21 <- ## barplot das ASVs encontradas em cada ponto
    alpha_all_ID_tbl %>% 
    filter(Class == "Actinopteri") %>%
    filter(!is.na(`Curated ID`)) %>%   # retirar os NA
    group_by(Sample,
             Ponto,
             Expedition,
             Mês) %>% 
    filter(Expedition %in% c(
      # "Nov_Dec/20",
      # "Nov/20",
      # "Dec/20"
      "out/21", # deixando apenas as amostras de 2021
      "Nov/21"
    )) %>%
    mutate(Ponto = factor(Ponto, levels = c("L1",
                                            "L2",
                                            "L3",
                                            "L4"
    ))) %>%
    summarise(`Num ASVs`= sum(`Num ASVs`)) %>%
    ungroup() %>% 
    ggplot(aes(x = Ponto,
               y = `Num ASVs`,
               fill = Ponto)) +
    geom_bar(stat = "identity") +
    ylim(0, 85) +
    guides(col = guide_legend(nrow = 6)) +
    xlab("Ponto") +
    ylab("Número de ASVs") +
    ggtitle(label = "Número de ASVs por ponto (2021)") +
    theme(legend.position = "none") +
    geom_text(stat='identity', ## codigo que exibe dentro das barras a contagem
              aes(label=`Num ASVs`),
              position = position_dodge(width = 0.9),
              vjust = 1.4,
              colour = "#ffffff", size = 6) +
    facet_grid2(cols = vars(Expedition), ##facetando por mes
                scales = "free_x",
                axes = "all",
                space = "free_x",
                independent ='x') +
    theme(axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 23),
          axis.title.y = element_text(size = 23),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          strip.text = element_text(size = 19)) +
    scale_fill_manual(values = viridis::viridis(n=10)[c(1,7,5,9)])
  
  ## Plotando
  ggsave(plot = asvs_ponto_21, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/asv_ponto_21.pdf",
         device = "pdf", units = "cm", height = 15, width = 20, dpi = 600) 
  }  
  ## 2020
  {
    asvs_ponto_20 <- ## barplot das ASVs encontradas em cada ponto
      alpha_all_ID_tbl %>% 
      filter(Class == "Actinopteri") %>%
      filter(!is.na(`Curated ID`)) %>%   # retirar os NA
      group_by(Sample,
               Ponto,
               Expedition,
               Mês) %>% 
      filter(Expedition %in% c(
        "Nov_Dec/20",
        "Nov/20",
        "Dec/20"
        # "out/21", # deixando apenas as amostras de 2021
        # "Nov/21"
      )) %>%
      mutate(Ponto = factor(Ponto, levels = c("L1",
                                              "L2",
                                              "L3",
                                              "L4"
      ))) %>%
      summarise(`Num ASVs`= sum(`Num ASVs`)) %>%
      ungroup() %>% 
      ggplot(aes(x = Ponto,
                 y = `Num ASVs`,
                 fill = Ponto)) +
      geom_bar(stat = "identity") +
      ylim(0, 35) +
      guides(col = guide_legend(nrow = 6)) +
      xlab("Ponto") +
      ylab("Número de ASVs") +
      ggtitle(label = "Número de ASVs por ponto (2020)") +
      theme(legend.position = "none") +
      geom_text(stat='identity', ## codigo que exibe dentro das barras a contagem
                aes(label=`Num ASVs`),
                position = position_dodge(width = 0.9),
                vjust = 1.4,
                colour = "#ffffff", size = 6) +
      facet_grid2(cols = vars(Expedition), ##facetando por mes
                  scales = "free_x",
                  axes = "all",
                  space = "free_x",
                  independent ='x') +
      theme(axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
            axis.text.y = element_text(size = 20),
            axis.title.x = element_text(size = 23),
            axis.title.y = element_text(size = 23),
            plot.title = element_text(size = 22),
            legend.text = element_text(size = 19),
            strip.text = element_text(size = 19)) +
      scale_fill_manual(values = viridis::viridis(n=10)[c(1,7,5,9)])
    
    ## Plotando
    ggsave(plot = asvs_ponto_20, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/asv_ponto_20.pdf",
           device = "pdf", units = "cm", height = 15, width = 20, dpi = 600) 
  } 

## Tabela ssp por ponto/amostra
## nova tabela agrupando Curated ID por ponto e por amostra
{
  alpha_spp <- 
    alpha_all_ID_tbl %>% 
    filter(Class == "Actinopteri") %>%
    filter(!is.na(`Curated ID`)) %>%
    mutate(Point = factor(Point, levels = c("L1",
                                            "L2",
                                            "L3",
                                            "L4"
    ))) %>%
    group_by(Point,
             Sample,
             Expedition) %>% 
    summarise(num_especies = n_distinct(`Curated ID`))
}

## Barplot spp por ponto para comparar com as ASVs
## 2021
{
  spp_ponto_21 <- ## barplot das SSP encontradas em cada ponto
    alpha_spp %>% 
    filter(Expedition %in% c(
      # "Nov_Dec/20",
      # "Nov/20",
      # "Dec/20"
      "out/21", # deixando apenas as amostras de 2021
      "Nov/21"
    )) %>%
    ggplot(aes(x = Point,
               y = num_especies,
               fill = Point)) +
    geom_bar(stat = "identity") + 
    ylim(0, 85) +
    guides(col = guide_legend(nrow = 6)) +
    xlab("Ponto") +
    ylab("Número de espécies") +
    ggtitle(label = "Número de espécies por ponto (2021)") +
    theme(legend.position = "none") +
    geom_text(stat='identity', ## codigo que exibe dentro das barras a contagem
              aes(label= num_especies),
              position = position_dodge(width = 0.9),
              vjust = 1.4,
              colour = "#ffffff", size = 6) +
    facet_grid2(cols = vars(Expedition), ##facetando por mes
                scales = "free_x",
                axes = "all",
                space = "free_x",
                independent ='x') +
    theme(axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 23),
          axis.title.y = element_text(size = 23),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          strip.text = element_text(size = 19)) +
    scale_fill_manual(values = viridis::viridis(n=10)[c(1,7,5,9)])
  
  ## Plotando
  ggsave(plot = spp_ponto_21, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/spp_ponto_21.pdf",
         device = "pdf", units = "cm", height = 15, width = 20, dpi = 600)
}
## 2020
{
  spp_ponto_20 <- ## barplot das SSP encontradas em cada ponto
    alpha_spp %>% 
    filter(Expedition %in% c(
      "Nov_Dec/20",
      "Nov/20",
      "Dec/20"
      # "out/21", # deixando apenas as amostras de 2021
      # "Nov/21"
    )) %>%
    ggplot(aes(x = Point,
               y = num_especies,
               fill = Point)) +
    geom_bar(stat = "identity") + 
    ylim(0, 35) +
    guides(col = guide_legend(nrow = 6)) +
    xlab("Ponto") +
    ylab("Número de espécies") +
    ggtitle(label = "Número de espécies por ponto (2020)") +
    theme(legend.position = "none") +
    geom_text(stat='identity', ## codigo que exibe dentro das barras a contagem
              aes(label= num_especies),
              position = position_dodge(width = 0.9),
              vjust = 1.4,
              colour = "#ffffff", size = 6) +
    facet_grid2(cols = vars(Expedition), ##facetando por mes
                scales = "free_x",
                axes = "all",
                space = "free_x",
                independent ='x') +
    theme(axis.text.x = element_text(size = 20), ## definindo o tamanho de cada texto
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 23),
          axis.title.y = element_text(size = 23),
          plot.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          strip.text = element_text(size = 19)) +
    scale_fill_manual(values = viridis::viridis(n=10)[c(1,7,5,9)])
  
  ## Plotando
  ggsave(plot = spp_ponto_20, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/qualificacao/spp_ponto_20.pdf",
         device = "pdf", units = "cm", height = 15, width = 20, dpi = 600)
}
  
  # Dados de exemplo
  dados <- data.frame(
    Campanha = c(rep("Campanha 1", 5), rep("Campanha 2", 5)),
    Ponto = c(1, 2, 3, 4, 5, 1, 2, 3, 4, 5),
    AlphaDiversidade = c(10, 8, 12, 9, 11, 7, 9, 8, 10, 6)
  )
  
  # Criar o gráfico
  ggplot(dados, aes(x = Ponto, y = AlphaDiversidade, color = Campanha, shape = Campanha)) +
    geom_point(size = 3) +
    labs(x = "Ponto", y = "Alpha Diversidade", color = "Campanha", shape = "Campanha") +
    scale_color_manual(values = c("blue", "red")) +
    scale_shape_manual(values = c(16, 17)) +
    theme_minimal()
  