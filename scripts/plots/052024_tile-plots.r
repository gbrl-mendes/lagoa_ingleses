---
title: "Tile plots Lagoa dos Ingleses"
author: "Gabriel Mendes"
date: "05/2024"
---

## Carregando bibliotecas ----
{
  required_packages <- c("Biostrings",
                         "DECIPHER",
                         "factoextra",
                         "future",
                         "ggh4x",
                         "ggpubr",
                         "Matrix",
                         "stringr",
                         "vegan",
                         "tidyverse",
                         "readxl",
                         "ggplot2")
  
  for (package in required_packages){
    if(!require(package,character.only = TRUE)) {
      install.packages(package, dependencies=TRUE)
    }
    library(package,character.only = TRUE)
  }
}
## Caminhos ----
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figuras")
  tbl_path <- paste0(prjct_path,"tabelas")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

# Obtencao dos dados ----

# Usando os objetos criados pelo codigo 052024_post_pipe.r:
# 1. grouped_filt

## Organizar as especies ----

grouped_filt$`Curated ID` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
grouped_filt$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  
# Organizar ordem

subset(grouped_filt, `Classe` == "Actinopteri")$`Order` %>%
  sort() %>% unique() # ver ordens
    
fish_ordens <- c(#"Acanthuriformes", # organizar as ordem das ordens de peixes
                 "Characiformes",
                 "Cichliformes",
                 "Clupeiformes",
                 "Cypriniformes",
                 "Cyprinodontiformes",
                 "Gymnotiformes",
                 "Salmoniformes",
                 "Siluriformes",
                 "Synbranchiformes")                       
    
subset(grouped_filt, `Classe` != "Actinopteri")$`Classe` %>%
  sort() %>% unique() # ver classes

nfish_classes <- c("Aves", # ordem das classes de nao peixes
                   "Mammalia")
  
  
# Organizar classes
  
subset(grouped_filt, `Classe` != "Actinopteri")$`Order` %>%
  sort() %>% unique() # ver ordens
    
nfish_ordens <- c(#"Anseriformes", 
                  "Artiodactyla",
                  #"Carnivora",
                  #"Cingulata",
                  #"Didelphimorphia",
                  "Galliformes",
                  "Lagomorpha",
                  #"Passeriformes",
                  "Primates",
                  "Rodentia"
                  #"Suliformes"
                  )
  
# Definir quais serao as especies e ordenar
  
grouped_filt$Sample %>% unique()

    spp_total <- grouped_filt %>%
      # filter(Sample %in% c("L1_nov_dec_20_mi",
      #                      "L2_nov20",
      #                      "L2_dez20")) %>% 
      # filter(`Classe` == "Actinopteri") %>%
      # filter(`Classe` == "Mammalia") %>%
      # mutate(words_count = str_count(`Curated ID`, "\\w+")) %>%
      # filter(words_count >= 2) %>% #so a nivel de spp
      # filter(str_detect(`Curated ID`, "sp.")) %>%
      select(`Curated ID`) %>%
      unique() %>%
      arrange(`Curated ID`) %>%
      pull()
    list(spp_total)
    
    {
    spp_fish <- grouped_filt %>%  # especies de peixes
      # filter(`Primer expected length` == "in range") %>%
      # filter(`Curated ID` != c("NA","")) %>%
      filter(`Classe` == "Actinopteri") %>%
      # mutate(words_count = str_count(`Curated ID`, "\\w+")) %>%
      # filter(words_count >= 2) %>% #so a nivel de spp
      select(`Curated ID`) %>%
      unique() %>% 
      arrange(`Curated ID`) %>%
      pull()
      
    list(spp_fish)
    }
  
  {
    no_fish <- grouped_filt %>%  # outras especies
      filter(`Classe` != "Actinopteri") %>% 
      # filter(`Curated ID` != "") %>% 
      select(`Curated ID`) %>% 
      unique() %>% 
      arrange(`Curated ID`) %>% 
      pull()
    
    list(no_fish)
  }  
  
  # Definir quais serao as amostras e ordenar
  {
    grouped_filt$Sample %>% unique()
    
    samples <- c("L1_nov_dec_20_mi",
                 "L2_dez20",
                 "L2_nov20",
                 "L1_out21",
                 "L2_out21",
                 "L3_out21",
                 "L4_out21",
                 "L1_nov21",
                 "L2_nov21",
                 "L3_nov21",
                 "L4_nov21",
                 "STX_L1_nov21",
                 "STX_L2_nov21",
                 "STX_L3_nov21",
                 "STX_L4_nov21",
                 "L1_jan22",
                 "L2_jan22",
                 "L3_jan22",
                 "L4_jan22"
                 )
  }
  
  # Definir quais serao as campanhas e ordenar
  {
    grouped_filt$expedition %>% unique()
    
    expeditions <- c("Novembro e Dezembro 2020",
                     "Novembro 2020",
                     "Dezembro 2020",
                     "Outubro 2021",
                     "Novembro 2021",
                     "Janeiro 2022")
  }
  
  # tabela so peixes
  {
    fish_ID_tbl <- grouped_filt %>% # ids apenas os peixes a nivel de spp
      filter(`Classe` == "Actinopteri") %>%
      # mutate(`Curated ID` = factor(`Curated ID`,
      #                              levels = rev(spp_fish))) %>%
      # filter(`Curated ID` %in% spp_fish) %>% 
      arrange(`Curated ID`) %>% 
      # correcao do RRA ja que as ASVs sem ID, de nao-peixes e <100 foram removidas alterando a abd total
      group_by(Sample) %>%
      mutate(RRA = RRA / sum(RRA)) %>%
      ungroup()
  
    View(fish_ID_tbl)  
  
    # conferir se a correcao do RRA funcionou
    fish_ID_tbl %>%
    # grouped_filt %>%
      group_by(Sample) %>% 
      summarize(total_RRA = sum(RRA))
    }
  
  # tabela nao peixes
  {
    nfish_ID_tbl <- grouped_filt %>% # ids nao peixes a nivel de spp
      filter(`Classe` == c("Aves", "Mammalia")) %>% 
      # mutate(`Order (DADA2)` = factor(`Order (DADA2)`,
      #                         levels = nfish_ordens)) %>% 
      # filter(`Curated ID` != "") %>% #tirei la em cima
      # filter(!is.na(`Curated ID`)) %>% #tirei la em cima
      arrange(`Curated ID`) %>% 
      # correcao do RRA ja que as ASVs sem ID e de peixes foram removidas alterando a abd total
      group_by(Sample) %>%
      mutate(RRA = RRA / sum(RRA)) %>%
      ungroup() %>% 
      select(all_of(names(fish_ID_tbl)))
    
    View(nfish_ID_tbl)
    
    # conferir se a correcao do RRA funcionou
    nfish_ID_tbl %>% 
      group_by(Sample) %>% 
      summarize(total_RRA = sum(RRA))
  }
  
    # tabela peixes & outros
  {
    all_ID_tbl <- fish_ID_tbl %>% 
      rbind(nfish_ID_tbl)
      
      View(all_ID_tbl)
  }


## Tile Plots ----

# Criacao do Tile Plot das amostras da Lagoa dos Ingleses 

  # Todos anos (apenas MCE)
    tile_all <- fish_ID_tbl %>%
      # filter(RRA >= 0.01) %>%
      mutate(`Curated ID` = factor(`Curated ID`, levels = rev(spp_fish))) %>%
      mutate(expedition = factor(expedition, levels = expeditions)) %>%
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(new_name = factor(new_name)) %>% 
      # mutate(sample = factor(sample)) %>% 
      mutate(`Order` = factor(`Order`)) %>%
      mutate(ano_nivel = case_when(
        year == "2020" ~ "2020 - Full",
        year == "2021" ~ "2021 - Empty",
        year == "2022" ~ "2022 - Full"
      )) %>%
      group_by(expedition, `Curated ID`, new_name, `Order`, Nivel) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(expedition %in% c("Novembro 2021")) %>%
      filter(filter %in% c("MCE")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`,
                 x = new_name,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # stat = "identity",
      # colour = "black", size = 4) +
      facet_grid(cols = vars(ano_nivel # Facet por ano
                             # ,filter
                             ), 
                 rows = vars(`Order`),
                 space = "free", 
                 scales = "free",
                 drop = TRUE) +
      scale_fill_continuous(
        trans = "log10",
        breaks = c(0.000001, 0.0087, 0.87),  # Defina os pontos de quebra desejados
        labels = c("0.001", "0.01", "100"),  # Rótulos correspondentes
        type = "viridis"
      ) +
      # theme_minimal() +
      theme(
        panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1),
        # axis.ticks = element_line(color = "grey"),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = rel(1.2)),
        # legend.background = element_blank(),
        # legend.key = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black", size = rel(1.2)),
        plot.title = element_text(color = "black", size = rel(1.5)),
        plot.subtitle = element_text(color = "black", size = rel(1.2)),
        strip.background = element_rect(fill = "#e4e4e4"),
        strip.text = element_text(color = "black", size = rel(1.2))) +
      theme(plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(size = 14, angle = 20, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 14, face = "italic"),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text = element_text(size = 14, face = "bold"),
            strip.text.y = element_text(angle=0),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14),
            legend.position = "bottom"
      ) +
      labs(x = "Sample sites",
           y = "Espécies",
           fill ='RRA (%)',
           title = "Detected species",
           subtitle = "Use of MCE filter")
    
    tile_all
    
ggsave(plot = tile_all, 
       filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2024/",
                        "tile_20-22_MCE", "-", Sys.Date(), ".pdf", sep = ""),
       device = "pdf", units = "cm", height = 40, width = 40, dpi = 600)
    
  # MCE versus Sterivex ----
  {
    mce_vs_sterivex <- fish_ID_tbl %>%
      # filter(RRA >= 0.01) %>%
      mutate(`Curated ID` = factor(`Curated ID`, levels = rev(spp_fish))) %>%
      mutate(expedition = factor(expedition, levels = expeditions)) %>%
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(new_name = factor(new_name)) %>% 
      # mutate(sample = factor(sample)) %>% 
      mutate(`Order` = factor(`Order`)) %>%
      mutate(ano_nivel = case_when(
        year == "2020" ~ "2020 - Cheio",
        year == "2021" ~ "2021 - Vazio",
        year == "2022" ~ "2022 - Cheio"
      )) %>%
      group_by(expedition, `Curated ID`, new_name, `Order`, Nivel) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(year %in% c(2021)) %>%
      filter(expedition %in% c("Novembro 2021")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`,
                 x = new_name,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # stat = "identity",
      # colour = "black", size = 4) +
      facet_grid(cols = vars(filter),
                 rows = vars(`Order`),
                 space = "free", 
                 scales = "free",
                 drop = TRUE) +
      scale_fill_continuous(
        trans = "log10",
        breaks = c(0.000001, 0.0087, 0.87),  # Defina os pontos de quebra desejados
        labels = c("0.001", "0.01", "100"),  # Rótulos correspondentes
        type = "viridis"
      ) +
      # theme_minimal() +
      theme(
        panel.background = element_blank(),
        # panel.grid.major = element_blank(),
        panel.grid.major = element_line(color = "grey",
                                        size = 0.2,
                                        linetype = 1),
        # axis.ticks = element_line(color = "grey"),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black", size = rel(1.2)),
        # legend.background = element_blank(),
        # legend.key = element_blank(),
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black", size = rel(1.2)),
        plot.title = element_text(color = "black", size = rel(1.5)),
        plot.subtitle = element_text(color = "black", size = rel(1.2)),
        strip.background = element_rect(fill = "#e4e4e4"),
        strip.text = element_text(color = "black", size = rel(1.2))) +
      theme(plot.title = element_text(size = 20, face = "bold"),
            axis.text.x = element_text(size = 14, angle = 20, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 14, face = "italic"),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            strip.text = element_text(size = 14, face = "bold"),
            strip.text.y = element_text(angle=0),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14),
            legend.position = "bottom"
            ) +
      labs(x = "Pontos amostrais",
           y = "Espécies",
           fill ='Abundância \nrelativa (%)',
           title = "Espécies detectadas",
           subtitle = "Filtro MCE versus Sterivex")
    
    mce_vs_sterivex
  }
  ggsave(plot = mce_vs_sterivex, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2024/",
                          "tile_MCE_Sterivex", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf", units = "cm", height = 40, width = 30, dpi = 600)

  # vazio versus cheio (versao SBG) ----
  
  {
    cheio_vs_vazio <- fish_ID_tbl %>%
      # filter(RRA >= 0.01) %>%
      mutate(expedition = factor(expedition, levels = expeditions)) %>% 
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(new_name = factor(new_name)) %>% 
      # mutate(sample = factor(sample)) %>% 
      mutate(`Order` = factor(`Order`)) %>%
      group_by(expedition, `Curated ID`, new_name, `Order`, Nivel) %>%   # Agrupa por expedição, Curated ID e Point
      # filter(new_name %in% c("Ponte", "Fundacao")) %>%
      ggplot(aes(y = `Curated ID`, 
                 group = `Curated ID`,
                 x = new_name,
                 # x = Sample,
                 fill = RRA)) +
      geom_tile() +
      # geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # stat = "identity",
      # colour = "black", size = 4) +
      facet_grid(cols = vars(Nivel),
                 labeller = labeller(Nivel = c("Cheio" = "Full pond", "Vazio" = "Empty pond")),
                 rows = vars(`Order`),
                 space = "free", 
                 scales = "free",
                 drop = TRUE) +
      labs(fill ='Relative \nabundance (%)',
           x = "Sample points",
           y= "Species") +
      # ggtitle(label = "Espécies detectadas: Lagoa cheia e vazia") +
      scale_fill_continuous(
        trans = "log10",
        # breaks = c( 0.001, 1, 75),
        # breaks = c(10, 90),
        # type = "viridis"
        breaks = c(0.01, 0.1, 0.5, 1),  # Defina os pontos de quebra desejados
        labels = c("1%", "10%", "50%", "100%"),  # Rótulos correspondentes
        type = "viridis"
      ) +
      theme(plot.title = element_text(size = 25, face = "bold", hjust=0.7),
            axis.text.x = element_text(size = 13, face = "bold", angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 13, face = "bold"),
            axis.title.x = element_text(size = 16, face = "bold"),
            axis.title.y = element_text(size = 16, face = "bold"),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 16),
            strip.text = element_text(size = 13, face = "bold")) +
      theme(strip.text.y = element_text(angle=0)) +
      theme(legend.position = "bottom") 
    
    cheio_vs_vazio
  }
  ggsave(plot = cheio_vs_vazio, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/reanalise/",
                          "tile_plot_cheio_vs_vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf", units = "cm", height = 20, width = 20, dpi = 600)
  