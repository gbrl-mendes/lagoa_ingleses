---
title: "Alpha diversity and NMDS plots Lagoa dos Ingleses "
author: "Gabriel Mendes"
date: "08/2023"
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
}

# Caminhos ----
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figuras")
  tbl_path <- paste0(prjct_path,"/tabelas/raw/run_2_4_5")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

# Obtencao de dados ----
{
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v062023.csv"), sep = ",", check.names = FALSE) %>% tibble()
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Hoplias brasiliensis")] <- "Hoplias sp." # Ver artigo do Heron 2023
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Hoplias intermedius")] <- "Hoplias sp." # Ver artigo do Heron 2023
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Hoplias spp.")] <- "Hoplias sp." # Corrigir o nome
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Astyanax fasciatus")] <- "Psalidodon fasciatus" # A. fasciatus. por Psalidodon fasciatus
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Pimelodus fur")] <- "Pimelodus sp." # Baixa resolucao pra Pimelodus
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Pimelodus maculatus")] <- "Pimelodus sp." # Baixa resolucao pra Pimelodus
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Pimelodus pohli")] <- "Pimelodus sp." # Baixa resolucao pra Pimelodus
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Pimelodus spp.")] <- "Pimelodus sp." # Corrigir o nome

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

# Agrupar por ID e Blast ID
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
      "OTU",
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
              `OTU`,
              `ASV (Sequence)`) %>%
    ungroup()
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

# tabela so peixes
{
  fish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids apenas os peixes a nivel de spp
    mutate(`Curated ID` = factor(`Curated ID`,
                                 levels = rev(spp_fish))) %>% 
    mutate(`Order` = factor(`Order`,
                            levels = fish_ordens)) %>%
    mutate('Nivel' = ifelse(`Year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                            "Cheio", "Vazio")) %>%
    filter(`Curated ID` != "") %>%
    filter(!is.na(`Curated ID`)) %>% 
    arrange(`Curated ID`)
}

# Tabela resumo runs 2,4 e 5

resume_tbl <- fish_ID_tbl %>%
  filter(RRA >= 1) %>%
  group_by(Sample, New_name, Year) %>% 
  summarise("Num. ASVs" = length(unique(`ASV (Sequence)`)),
            "Num. OTUs" = length(unique(as.character(`OTU`))),
            "Species" = length(unique(`Curated ID`)),
  ) %>% 
  ungroup()
  

# Longer table
{
  longer_tbl_IDS <- 
    fish_ID_tbl %>%
    filter(RRA >= 1) %>% # definicao de qual threshold de abundancia sera usado
    select(c("Sample", "Curated ID", "Point", "Year", "Nivel", "New_name", "Abundance","RRA")) %>%
    unique() %>% 
    mutate("New_name" = as.factor(New_name))
    
  # %>%
    # filter(Sample %in% c("LI1-neo-mi", "L2_nov20", "L2_dez20", "L2_nov21", "L3_nov21", "L4_nov21")) #com apenas as amostras que Daniel pediu
  # %>% 
    # rename("Nivel" = "Nivel") %>%
    # mutate("Nivel" = as.factor(Nivel))
}

## Rarefaction ----

# Ref: https://www.youtube.com/watch?v=_OEdFjc1D9I

my_rarefy <- function(x, sample){
  
  x <- x[x>0]
  sum(1-exp(lchoose(sum(x) - x, sample) - lchoose(sum(x), sample)))
  
}

# Formatando a longer_tbm_IDS para entrar no VEGAN

wider_tbl_IDS <-
  longer_tbl_IDS[,c(1,2,7)] %>% 
  group_by(Sample, `Curated ID`) %>%
  summarise(total_abundance = sum(Abundance)) %>% # tem que juntar as ASVs que sao da mesma especie
  pivot_wider(names_from = `Curated ID`, 
              # values_from = RRA, # Nao pode ser RRA, tem que ser valores inteiros
              values_from = total_abundance,
              values_fill = 0) %>% 
  as.data.frame()

rownames(wider_tbl_IDS) <- wider_tbl_IDS$Sample

wider_tbl_IDS <- wider_tbl_IDS[,-1]

# Apos formatar a tabela, ela esta pronta para entrar no VEGAN.
# O primeiro passo e' descobrir qual a menor amostra

min_n_seqs <- fish_ID_tbl %>% 
  group_by(Sample) %>% 
  summarise(sample_total = sum(Abundance)) %>% 
  summarise(min = min(sample_total)) %>% 
  pull(min)

# A menor amostra é L4_out21, com 114 reads.
# Vou rarefar por ela e depois excluir ela e entao realizar a 
# rarefacao pela entao amostra menor. Desta maneira terei dois 
# resultados distintos e poderei comparar.

# Rarefacao

vegan::rarefy(wider_tbl_IDS, min_n_seqs) %>% 
  as_tibble(rownames = "Sample") %>% 
  select(Sample, vegan=value) # numero de spp esperado observar
                              # se todas as amostras tivesse uma
                              # profundidade de reads igual a L4_out21

## Alpha diversity ----

## Alpha (within sample) diversity ##

# Common alpha diversity statistics include:
# Shannon: How difficult it is to predict the identity of a randomly chosen individual.
# Simpson: The probability that two randomly chosen individuals are the same species.
# Inverse Simpson: This is a bit confusing to think about. Assuming a theoretically 
# community where all species were equally abundant, this would be the number of 
# species needed to have the same Simpson index value for the community being analyzed.


# Fonte: https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html

## Calculo usando o tidyverse
{
  # Funcoes para calculo de riqueza, indice de shannon e simpson
  {
    # Funcao para calcular a riqueza de spp
    richness <- function(x){
      
      sum(x > 0)
    }
    
    # Funcao para calcular o indice de shannon
    shannon <- function(x){
      
      rabund <- x[x>0] /sum(x)
      -sum(rabund * log(rabund))
    }
    
    # Funcao para calcular o indice de simpson
    simpson <- function(x){
      
      n <- sum(x)
      
      sum(x*(x-1)/(n*(n-1)))  
    }
    
  }
  
# PLots from Riffomonas Project  
  longer_tbl_IDS %>% 
    group_by(Sample) %>% 
    summarise(sobs = richness(RRA),
              shannon = shannon(RRA)
              ,
              simpson = simpson(RRA),
              invsimpson = 1/simpson,
              n = sum(RRA)) %>% 
    pivot_longer(cols = c(sobs, shannon, invsimpson, simpson), 
                 names_to = "metric") %>% 
    ggplot(aes(x = n,
               y = value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow = 4, scales = "free_y")
  }

## Usando o Vegan
{
  longer_tbl_IDS %>% 
    group_by(Sample) %>% 
    summarise(sobs = specnumber(RRA),
              shannon = diversity(RRA, index = "shannon"),
              simpson = diversity(RRA, index = "simpson"),
              invsimpson = 1/simpson,
              n = sum(RRA)) %>% 
    pivot_longer(cols = c(sobs, shannon, invsimpson, simpson), 
                 names_to = "metric") %>% 
    ggplot(aes(x = n,
               y = value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow = 4, scales = "free_y")
}

## Outros plots
{
  index_tbl <- longer_tbl_IDS %>% 
    group_by(Sample) %>% 
    summarise(sobs = richness(RRA),
              shannon = shannon(RRA),
              simpson = simpson(RRA),
              invsimpson = 1/simpson,
              n = sum(RRA))
  
  index_merge_tbl <- longer_tbl_IDS %>% 
    select(Sample, Point, Year, Nivel, New_name) %>% 
    unique() %>% 
    merge(index_tbl, by = "Sample") %>% 
    mutate(Sample = factor(Sample, levels = c("LI1-neo-mi",
                                              "L2_nov20",
                                              "L2_dez20",
                                              "L1_out21",
                                              "L2_out21",
                                              "L3_out21",
                                              "L4_out21",
                                              "L1_nov21",
                                              "L2_nov21",
                                              "L3_nov21",
                                              "L4_nov21"))) %>%
    filter(Sample %in% c("LI1-neo-mi",
                         "L2_nov20",
                         "L2_dez20",
                         "L1_out21",
                         "L2_out21",
                         "L3_out21",
                         "L4_out21"
                         ,
                         "L1_nov21",
                         "L2_nov21",
                         "L3_nov21"
                         ,
                         "L4_nov21"
                         ))
  ## Scatterplot
  index_merge_tbl %>% 
    ggplot(aes(x = Sample,
               y = sobs,
               color = Nivel)) +
    geom_point() +
    theme_minimal()
  
  ## Boxplot
  
  library(ggpubr)
  
  # Calculo do P-value de comparacao entre os dois grupos Cheio e Vazio
  index_merge_tbl_df <- index_merge_tbl[,c(1,4,6)] %>% as.data.frame()
  compare_means(sobs ~ Nivel, data=index_merge_tbl_df, method = "t.test", paired = FALSE)
  
  boxplot_alfa <- index_merge_tbl %>%
    mutate(Nivel = ifelse(Nivel == "Cheio", "Full", "Empty")) %>%
    mutate(Nivel = factor(Nivel, levels = c("Full", "Empty"))) %>% 
    ggplot(aes(x = Nivel,
               y = sobs,
               fill = Nivel)) + 
    geom_boxplot(
      # fill = "blue",
      color = "black",
      size = 0.3) + 
    geom_jitter(height = 0.2, width = 0.1, color = "gray") +
    labs(x = "Pound level", y = "Richness", title = "Riqueza de Espécies:\nLagoa cheia vs. nível baixo") +
    theme_minimal() +
    theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.3),
          axis.text = element_text(size = 13, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "right") +
    scale_fill_manual(values = c("#00BFC4", "#F8766D")) +
    # P-value no grafico
    stat_compare_means()
  
}  
    # Salvar em pdf
  ggsave(plot =  boxplot_alfa, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/alfa_d/",
                          "alfa-d_boxplot_cheio-vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 15,
         width = 10,
         dpi = 600) 

  ## Versao do heron
  # diversidade na mao

    # Calculo do P-value de comparacao entre os dois grupos Cheio e Vazio
  index_df <- fish_ID_tbl %>%
    group_by(Sample,Nivel) %>% 
    summarise("Species" = length(unique(`Curated ID`))) %>%
    mutate(Nivel = ifelse(Nivel == "Cheio", "Full", "Empty")) %>%
    mutate(Nivel = factor(Nivel, levels = c("Full", "Empty"))) %>% 
    as.data.frame()
  
  compare_means(Species ~ Nivel, data=index_df, method = "t.test", paired = FALSE)

  boxplot_heron <- fish_ID_tbl %>% 
    group_by(Sample,Nivel) %>% 
    summarise("Species" = length(unique(`Curated ID`))) %>%
    mutate(Nivel = ifelse(Nivel == "Cheio", "Full", "Empty")) %>%
    mutate(Nivel = factor(Nivel, levels = c("Full", "Empty"))) %>%
    ggplot(aes(x=Nivel,
               fill=Nivel,
               y=Species,
               alpha )) +
    geom_boxplot(
      # fill = "blue",
      color = "black",
      size = 0.3) +
    geom_jitter(height = 0.2, width = 0.1, color = "gray")+ 
    theme_minimal() +
    theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.3),
          axis.text = element_text(size = 13, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "right") +
    labs(x = "Level") +
    scale_fill_manual(values = c("#00BFC4", "#F8766D"))
  
  # Salvar em pdf
  ggsave(plot =  boxplot_heron, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/alfa_d/",
                          "alfa-d_boxplot_cheio-vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 15,
         width = 10,
         dpi = 600) 

## Diagrama de Venn
{
  library("ggvenn", lib.loc="/opt/R/4.1.0/lib/R/library")
  
  
  list_nivel_id <- 
    longer_tbl_IDS %>%
    group_by(Nivel) %>% 
    summarise(
      `Curated ID` = unique(`Curated ID`)) %>% 
    with(split(`Curated ID`, 
               factor(Nivel, levels = unique(Nivel)))) %>% 
  
  # venn_cheio_vazio <-
    list_nivel_id %>% 
    ggvenn(c("Cheio", "Vazio"),
           fill_color = c("#f8766d", "#00bfc4"),
           text_size = 6,
           set_name_size = 18,
           show_elements = TRUE,
           label_sep = "\n"
           # ,
           # auto_scale = TRUE
           )
  
  ggsave(plot = venn_cheio_vazio, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/alfa_d/",
                          "venn_cheio-vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 30,
         width = 45,
         dpi = 600)
  }


## NMDS ----

# NMDS (Non-Metric Multidimensional Scaling), ou Escalonamento Multidimensional Não 
# Métrico, é uma técnica de análise multivariada utilizada para visualizar a semelhança 
# ou dissimilaridade entre objetos ou observações em um conjunto de dados. O NMDS é 
# particularmente útil quando os dados não possuem uma estrutura linear clara e quando 
# a distância entre os pontos não pode ser representada de forma métrica.

# Converter planilha de identificações por ponto para o formato Amostras X IDs
{
  # Função necessária para juntar as abds da tabela em IDs
  sum_uniq <- function(vec) {
    if (is.character(vec)) {
      suniq <- BiocGenerics::unique(vec)
    }
    if (is.numeric(vec)) {
      suniq <- sum(vec)
    }
    return(suniq)
  }
  
  # Número de linhas (IDs diferentes) por amostra
  fish_ID_tbl$File_name %>% table()
  
  fish_ID_tbl %>% colnames()
  
  FINAL_tbl_IDs <- fish_ID_tbl %>%
    # filter(Year %in% c("2021")) %>% # definir qual ano entrara na analise
    filter(RRA >= 1) %>% # definicao de qual threshold de abundancia sera usado
    # filter(New_name %in% c("Ponte", "Fundacao")) %>% 
    # filter(Sample %in% c("LI1-neo-mi", "L2_nov20", "L2_dez20", "L2_nov21", "L3_nov21", "L4_nov21")) %>% #com apenas as amostras que Daniel pediu
    select(c("Sample", "Curated ID", "Point", "Year", "Nivel", "New_name", "RRA")) %>% 
    pivot_wider(id_cols = c("Sample", "Point", "Year", "Nivel", "New_name"),
                names_from = `Curated ID`,
                values_from = `RRA`,
                values_fn = sum_uniq,
                names_sort = TRUE,
                names_prefix = "ID_") %>% 
    relocate(c("Sample", "Point", "Year", "Nivel", "New_name", starts_with("ID_"))) %>%  
    mutate(across(starts_with("ID_"), replace_na, replace = 0)) %>% 
    mutate("New_name" = as.factor(New_name)) %>% 
    mutate("Nivel" = as.factor(Nivel)) 
  
  # Verificando a soma das linhas
  FINAL_tbl_IDs %>% select(starts_with(match = "ID_")) %>% rowSums(na.rm = TRUE)
  
  # TODO: Verificando a soma das colunas
}

#  Rodando o NMDS
{
  # Preparar dados para entrada no pacote vegan
  colnames(FINAL_tbl_IDs)
  
  FINAL_tbl_IDs$`Sample` %>% unique() %>% sort()
  
  all_IDs_NMDS_tbl <- FINAL_tbl_IDs %>% 
    mutate("Sample number" = 0) %>% 
    relocate("Sample number" )
  
  # Associar números de amostra aos nomes de amostra
  for (sample in 1:nrow(all_IDs_NMDS_tbl)) {
    all_IDs_NMDS_tbl$`Sample number`[sample] <- sample
  }
  
  # Ordenar dataframe usado no NMDS
  all_IDs_NMDS_df <- all_IDs_NMDS_tbl %>% as.data.frame() %>% 
    mutate(Nivel = ifelse(Nivel == "Cheio", "Full", "Empty"))
  
  
  # Nomear linhas como números de amostra e remover coluna
  row.names(all_IDs_NMDS_df) <- all_IDs_NMDS_df$`Sample number`
  
  # Corrigir nomes das espécies para evitar problemas na plotagem
  colnames(all_IDs_NMDS_df)
  colnames(all_IDs_NMDS_df)[7:ncol(all_IDs_NMDS_df)] <- colnames(all_IDs_NMDS_df)[7:ncol(all_IDs_NMDS_df)] %>%
    str_replace_all(pattern = " ", replacement = "_") %>% 
    str_replace_all(pattern = "\\.", replacement = "") %>% 
    str_replace_all(pattern = "\\(", replacement = "") %>% 
    str_replace_all(pattern = "\\)", replacement = "")
  
  # Executando o NMDS
  # Esta e a funcao que faz o NMDS. Para ela fornece-se apenas
  # as colunas relativas as especies nas amostras
  all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)],
                                   comm = all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)],
                                   distance = "bray" # Leva em conta o RRA
                                   # distance = "jaccard" # Considera apenas presença/ausenciaa
                                   )
    
  dim(all_IDs_NMDS_df)

  # Fazer o fit das variáveis ambientais
  meta.envfit <- envfit(all_ps_vegan_ord_meta, all_IDs_NMDS_df[, c("Nivel", "Sample")],
                        permutations = 999,
                        na.rm=TRUE) # this fits environmental vectors
  
  # Espécies 
  
  # Fazer o fit das espécies para identificar a significância delas na explicação dos agrupamentos. 
  # Esse é o passo que mais demora quando com muitas amostras e espécies.
  meta.spp.fit <- envfit(all_ps_vegan_ord_meta, 
                         all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)], 
                         permutations = 999) # this fits species vectors
  
  # Obter valores de p para as espécies
  sps_pvals <- tibble("IDs" = names(meta.spp.fit$vectors$pvals),
                      "p-value" = meta.spp.fit$vectors$pvals)
  
  spp.scrs <- as.data.frame(scores(meta.spp.fit, display = "vectors")) %>%
    mutate("IDs" = rownames(.)) %>%
    left_join(y = sps_pvals, by = "IDs")
  
  # Selecionar espécies significativas
  sig.spp.scrs <- spp.scrs
  # %>%
  #   filter(`p-value` <=
  #            0.05) # definir p-value aqui!
  
  # Pontos amostrais
  
  # Definir os valores de NMDS1 e NMDS2, e os metadados de cada amostra/ponto amostral
  site.scrs <- as.data.frame(scores(all_ps_vegan_ord_meta, display = "sites")) %>%
    mutate("Sample number" = as.double(row.names(.))) %>%
    left_join(y = all_IDs_NMDS_df[, c("Sample number",
                                      "Sample",
                                      "Nivel",
                                      "New_name"
                                      )], 
              by = "Sample number")
  
  # Determinar centroides
  scrs <- scores(all_ps_vegan_ord_meta, display = "sites")
  
  cent <- aggregate(scrs ~ Nivel, data = site.scrs, FUN = "mean")
  
  # Calcular elipses
  NMDS <- data.frame("MDS1" = all_ps_vegan_ord_meta$points[, 1],
                     "MDS2" = all_ps_vegan_ord_meta$points[, 2],
                     "Nivel" = as.factor(all_IDs_NMDS_df$Nivel), check.names = FALSE)
  
  NMDS.mean <- aggregate(NMDS[, 1:2], list(group = NMDS$Nivel), "mean")
}

# Elipses

{
  plot(all_ps_vegan_ord_meta)
  
  # Sobrepor as elipses
  ord <- ordiellipse(ord = all_ps_vegan_ord_meta, 
                     groups = all_IDs_NMDS_df$Nivel,
                     display = "sites",
                     kind = "ehull", conf = 0.95, label = T)
  
  
  # Funcao do vegan de calcular elipses
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ell <- data.frame()
  
  for(g in levels(NMDS$Nivel)){
    print(g)
    df_ell <- 
      rbind(df_ell, 
            cbind(as.data.frame(with(NMDS[NMDS$Nivel==g,],
                                                     veganCovEllipse(
                                                       ord[[g]]$cov,
                                                       ord[[g]]$center,
                                                       ord[[g]]$scale))),
                  Nivel=g))}
}

# Configurar plot NMDS

# Definir paleta de cores
{
  paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")
  paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[1:3]
  paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[c(1:3, 2,7,9,10)]
  
  my_cols <- c("BAM" = "#0F6B99FF",
               "PARAC" = "#6B990FFF",
               "SAM" = "#99540FFF",
               "SFC" = "#A3CC51FF",
               "SFI" = "#7EC3E5FF",
               "SFM" = "#E5B17EFF") 
  }

# Plot
{
  NMDS_cheio_vazio <-
    ggplot(data = site.scrs,
           aes(x=NMDS1, 
               y=NMDS2)) +
    
    # Elipses 
    ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                               data = df_ell,
                               aes(x = NMDS1,
                                   y = NMDS2,
                                   group = Nivel,
                                   label = Nivel,
                                   col = Nivel,
                                   fill = Nivel
                               ), 
                               alpha=0.10,
                               # n = 200,
                               linetype=2,
                               expand = 0,
                               label.fontsize = 18,
                               con.cap = 0.1
    ) +
    # Niveis amostrais
    geom_point(aes(x=NMDS1,
                   y=NMDS2,
                   fill = Nivel,
                   # label = `Sampling unit`,  #descomentar se quiser exibir os nomes dos `Sampling sites`s
                   col=Nivel,
                   group = Nivel,
                   shape = New_name),
               stroke = 0.5,
               alpha = 0.75,
               size = 3) +
    
    # Nomes dos `Sampling sites`s amostrais
    ggrepel::geom_text_repel(aes(label = Sample),  #descomentar bloco se quiser exibir os nomes dos pontos
                             # hjust=0.5,
                             # vjust=2.75,
                             size=4.5,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha=0.1,
                             # min.segment.length = 1.5,
                             force = 3,
                             max.overlaps = 100,
                             fontface = "bold") +
    
    # Vetores das IDs
    geom_segment(data = sig.spp.scrs, aes(x = 0,
                                          xend = NMDS1,
                                          y = 0,
                                          yend = NMDS2),
                 arrow = arrow(length = unit(0.1, "cm")),
                 colour = "grey30",
                 alpha = 0.5,
                 lwd = 0.3) + #add vector arrows of significant species
    
    # Nomes das IDs
    ggrepel::geom_text_repel(data = sig.spp.scrs,
                             aes(x=NMDS1, y=NMDS2, label = IDs),
                             size = 3.5,
                             alpha = 0.75,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha = 0.1,
                             max.overlaps = 100) +
    
    # Centroides
    geom_point(data = cent,
               aes(x = NMDS1,
                   y = NMDS2, # colour = Nivel,
                   fill = Nivel
               ),
               size = 8,
               colour = "#222222",
               alpha = 0.75,
               shape  = 13
    ) + 
    coord_fixed(expand = c(0.5))+
    theme_light() +
    # +
    # scale_colour_manual(values = c("#fcca03","#6b0000","#02cc37"))+
    # scale_colour_manual(values = c("#fcca03","#6b0000","#02cc37"))+
    # scale_fill_manual(values = my_cols)+
    # scale_colour_manual(values = my_cols) +
    # +
    # scale_shape_manual(values = c("MG" = 24,                                    # se quiser formas customizadas para algum metadado
    #                               "RJ" = 21,
    #                               "SP" = 22,
    #                               "Minas Gerais" = 24,
    #                               "Rio de Janeiro" = 21,
  #                               "São Paulo" = 22
  #                               )) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
    labs(
      title  = paste0("Composição da ictiofauna: lagoa vazia e lagoa cheia"),
      subtitle = paste0("Stress: ",format(round(all_ps_vegan_ord_meta$stress,4)))
    ) +
    theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(plot.subtitle = element_text(size = 16, face = "bold")) +
    theme(legend.title = element_text(size = 16)) +
    theme(legend.text =  element_text(size = 12)) +
    theme(axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 16, face = "bold")) +
    theme(legend.position = "bottom") +
    guides(shape= "none")
  
  NMDS_cheio_vazio
}
  
  # Salvar em pdf
  ggsave(plot =  NMDS_cheio_vazio, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/nmds/",
                          "nmds_cheio_vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 25,
         width = 15,
         dpi = 600)
                    
                        


