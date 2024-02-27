"Estatísticas"
"Mendes, GA; Hilário, OH"
"12/2023"

# Pacotes ----
{
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(base)
}

# Obtencao dos dados ----

# Usando os objetos criados pelo codigo 102023_tile-plots.r:
# 1. raw_results_tbl
# 2. grouped_by_ID_BLASTid

# Tabela suplementar 1 ----

## Tabela usada para avaliar os dados

# Tabela longer
dt_all_resume <- grouped_by_ID_BLASTid %>%
  group_by(Nivel, new_name) %>%
  mutate("Abd total periodo" = sum(Abundance)) %>%
  ungroup() %>%
  group_by(`Curated ID`,
           `Final ID (BLASTn)`,
           `new_name`, 
           # `Abundance`,
           `Nivel`
  ) %>%
  mutate("RRA no periodo" = (Abundance/`Abd total periodo`)*100) %>% 
  # ungroup() %>% 
  summarise(
    "RRA" = sum(`RRA no periodo`),
    "ASVs" = length(unique(`ASV header`)),
    "OTUs" = length(unique(OTU)),
    "Nivel" = Nivel,
    # "Ponto" = new_name,
    "Abundance" = sum(Abundance),
    "Class" = `Class (BLASTn)`,
    "Order" = `Order (BLASTn)`,
    "Family" = `Family (BLASTn)`,
    "Genus" = `Curated genus`
  ) %>% 
  mutate(
    "RRA_formatado" = format(RRA, scientific = TRUE)
  ) %>%
  unique()

View(dt_all_resume)

options(scipen = 99,
        digits = 3)

# Teste de RRA para ver se foi calculado corretamente
dt_all_resume %>%
  filter(Nivel == "Cheio") %>%
  # filter(Sample == "L2_nov20") %>%
  filter(new_name == "Fundação") %>%
  pull(RRA) %>%
  sum()  

# Tabela wider com raw data
wider_dt_all_resume <- dt_all_resume %>% 
  select(-c("RRA_formatado")) %>%
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  unite(new_name, Nivel,col= "ponto_nivel") %>% 
  pivot_wider(id_cols = c("Class","Curated ID", "Final ID (BLASTn)"),
              names_from = ponto_nivel,
              values_from =  c("RRA","ASVs","OTUs","Abundance"),
              names_glue = "{ponto_nivel}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate( "Class","Curated ID","Final ID (BLASTn)") 

View(wider_dt_all_resume)

write.csv(wider_dt_all_resume, "~/projetos/lagoa_ingleses/results/tabelas/wider_dt_all_resume.csv")

# Tabela definitiva 1 ----

## Apenas e' a Tabela suplementar 1 apenas com peixes
## e filtrada para Abundance > 100

# Tabela longer com os dados filtrados
dt_filt <- dt_all_resume %>% 
  filter(Class %in% "Actinopteri") %>%
  group_by(`Curated ID`, new_name, Nivel) %>% 
  summarise("Curated ID" = `Curated ID`,
            "Class" = Class,
            "Order" = Order,
            "Family" = Family,
            "Genus" = Genus,
            "Abundance" = sum(Abundance),
            "ASVs" = sum(ASVs),
            "OTUs" =sum(OTUs),
            "new_name" = new_name,
            "Nivel" = Nivel,
            "RRA" = sum(RRA)) %>% 
  ungroup() %>% 
  unique() %>%
  filter(Abundance >= 100)

View(dt_filt)

# Tabela wider com os dados filtrados
wider_dt_filt <- dt_filt %>% 
  mutate(RRA = round(RRA,digits = 4)) %>% 
  ungroup() %>% 
  unite(new_name, Nivel,col= "ponto_nivel") %>% 
  pivot_wider(id_cols = c("Curated ID"),
              names_from = ponto_nivel,
              values_from =  c("RRA","ASVs","OTUs","Abundance"),
              names_glue = "{ponto_nivel}_{.value}") %>%
  select(sort(colnames(.))) %>% 
  relocate("Curated ID") 

View(wider_dt_filt)

write.csv(wider_dt_filt, "~/projetos/lagoa_ingleses/results/tabelas/wider_dt_filt.csv")

# Analise das amostras ----

## Esta tabela sumariza, para cada classe, em cada ano, a quantidade de reads, o numero de ASVs, de OTUs, 
## de ids a nivel de spp e ids em outros niveis taxonomicos

analis_dt_filt <-
  dt_filt %>%
  # group_by(
    # `Curated ID`, 
           # new_name, 
           # Nivel
           # ) %>% 
  summarize(
    "Total reads" = sum(Abundance),
    "ASVs" = sum(ASVs),
    "OTUs" = sum(OTUs),
    "Orders" = length(unique(Order)),
    "Families" = length(unique(Family)),
    "Genera" = length(unique(Genus)),
    "Species" = length(unique(`Curated ID`[grepl("^[A-Za-z]+\\s[A-Za-z]+$", `Curated ID`) & !grepl("sp\\.", `Curated ID`)])),
    "Nspecies" = length(unique(`Curated ID`)) - Species
  ) %>%  
  ungroup()

# View(analis_pre_eco)
View(analis_dt_filt)

# Resumo 2

# Actinopteri
# Esta tabela mostra, por cada identificacao, a quantidade de reads atribuida, o numero de ASVs e de OTUs
stats_by_actinopteri <-
  grouped_by_ID_BLASTid %>%
  filter(`Class (BLASTn)` %in% "Actinopteri") %>%
  group_by(`Class (BLASTn)`,
           year,
           `Curated ID`
  ) %>%
  reframe(
    "Order" = `Curated Order`,
    "Especie" = `Curated ID`,
    "Reads" = sum(Abundance),
    "ASVs" = n_distinct(`ASV header`),
    "OTUs" = n_distinct(OTU)
  ) %>% 
  relocate(year, Order, Especie, `Reads`, ASVs, OTUs) %>% 
  filter(Reads >= 100) %>%
  unique()

  View(stats_by_actinopteri)
  
# 2020 Actinopteri ordens proporcoes

  fish_ord_2020 <- grouped_by_ID_BLASTid %>%
    filter(RRA >= 0.01) %>% 
    filter(Sample %in% c(
      "L1_nov_dec_20_mi",
      "L2_dez20",
      "L2_nov20"
      # "out/21", # deixando apenas as amostras de 2021
      # "Nov/21"
    )) %>%
    filter(`Class (BLASTn)` %in% "Actinopteri") %>%
    filter(`Primer expected length` %in% "in range") %>%
    ggplot(aes(x = new_name,
               y = length(`Curated ID`),
               fill = `Curated Order`)) +
    geom_col(position = "fill") +
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
    scale_y_continuous(labels = scales::percent) +
    xlab("Local") +
    ylab("Proporção")
  
  ggsave(plot = fish_ord_2020, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/defesa/alpha_plot_ordem_2020.pdf",
         device = "pdf", units = "cm", height = 12, width = 35, dpi = 600)
  
# 2020 Actinopteri ordens proporcoes
  
  fish_ord_2020_reads <- grouped_by_ID_BLASTid %>%
    mutate(Point = case_when(
      Sample == "L1_nov_dec_20_mi" ~ "L1",
      Sample == "L2_dez20" | Sample == "L2_nov20" ~ "L2"
    )) %>% 
    filter(Sample %in% c(
      "L1_nov_dec_20_mi",
      "L2_dez20",
      "L2_nov20"
      # "out/21", # deixando apenas as amostras de 2021
      # "Nov/21"
    )) %>%
    filter(`Class (BLASTn)` %in% "Actinopteri") %>%
    filter(`Primer expected length` %in% "in range") %>%
    ggplot(aes(x = Point,
               y = Abundance,
               fill = `Curated Order`)) +
    geom_col(position = "fill") +
    coord_flip() + #virando 90o o grafico
    ggtitle(label = "Proporção de reads atribuídas a cada ordem (2020)") +  ## alterando o titulo do grafico
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
    scale_y_continuous(labels = scales::percent) +
    xlab("Ponto") +
    ylab("Proporção")
  
  ggsave(plot = fish_ord_2020_reads, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/defesa/reads_plot_ordem_2020.pdf",
         device = "pdf", units = "cm", height = 12, width = 27, dpi = 600) 




