---
  title: "Estatísticas"
author: 
  - "Mendes, GA; Hilário, OH; Carvalho, DC"
date: "07/12/2021"
output:
  html_document:
  code_download: yes
theme: flatly
toc: true
toc_depth: 4
toc_float: true
pdf_document: default
editor_options:
  chunk_output_type: console
---
  
# Pacotes ----
{
  library(dplyr)
  library(tidyr)
  library(tidyverse)
  library(moments)
  library(ggplot2)
  library(base)
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
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v2.csv"), sep = ",", check.names = FALSE) %>% tibble()
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Oreochromis niloticus")] <- "Tilapia rendalli" # Oreochromis niloticus é Tilapia
}

# Filtrando os dados ----

## Alpha ttbl original
{
  alpha_tbl <- raw_results_tbl %>%
    # filter(`Relative abundance on sample` >= 0.01) %>% # heron pediu para manter as ASVs espurias
    filter(!`Curated ID` %in% c("Actinopteri", ## tirando as ASVs que nao foram identificadas a nivel de especie, nao-peixes e NA
                                "Astyanax",
                                "Characidae",
                                "Characidium",
                                "Characiformes",
                                "Cichla",
                                "Cichlidae",
                                "Hoplias",
                                "Pimelodus",
                                "",
                                "Cavia magna",
                                "Cutibacterium acnes",
                                "Bos taurus",
                                "Canis familiaris",
                                "Didelphis albiventris (Gamba)",
                                "Homo sapiens",
                                "Hydrochaeris hydrochaeris (Capivara)",
                                "Nannopterum brasilianus",
                                "Oryctolagus cuniculus (Coelho-bravo)",
                                "Progne chalybea (Andorinha-grande)",
                                "Sus scrofa"
    )) %>%
    mutate("Mês" = str_split_fixed(string = .$Sample, # criando uma nova coluna quebrando as infos da coluna Sample
                                   pattern = "_",
                                   n = 2)[,2]) %>%
    mutate(Mês = factor(Mês, levels = c("out21",
                                        "nov21"))) %>%
    filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                             "Nov/21")) %>%
    group_by(Sample,
             Read_origin,
             Primer,
             `Curated ID`,
             Year,
             Point,
             Expedition,
             Mês
    ) %>%
    summarize(`Num ASVs` = length(unique(`ASV (Sequence)`)),
              `Num OTUs` = length(unique(`OTU`)),
              `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
              `Ponto` = Point) %>%
    mutate(Sample = factor(Sample, levels = c("L1_out21",
                                              "L1_nov21",
                                              "L2_out21",
                                              "L2_nov21",
                                              "L3_out21",
                                              "L3_nov21",
                                              "L4_out21",
                                              "L4_nov21"
    ))) %>%
    ungroup() %>%
    unique()
}

## Separando por amostra
{
  sample_tbl <- raw_results_tbl %>%
    filter(!`Curated ID` %in% c("Actinopteri", ## tirando as ASVs que nao foram identificadas a nivel de especie, nao-peixes e NA
                                "Astyanax",
                                "Characidae",
                                "Characidium",
                                "Characiformes",
                                "Cichla",
                                "Cichlidae",
                                "Hoplias",
                                "Pimelodus",
                                "",
                                "Cavia magna",
                                "Cutibacterium acnes",
                                "Bos taurus",
                                "Canis familiaris",
                                "Didelphis albiventris (Gamba)",
                                "Homo sapiens",
                                "Hydrochaeris hydrochaeris (Capivara)",
                                "Nannopterum brasilianus",
                                "Oryctolagus cuniculus (Coelho-bravo)",
                                "Progne chalybea (Andorinha-grande)",
                                "Sus scrofa"
    )) %>%
    mutate("Mês" = str_split_fixed(string = .$Sample, # criando uma nova coluna quebrando as infos da coluna Sample
                                   pattern = "_",
                                   n = 2)[,2]) %>%
    mutate(Mês = factor(Mês, levels = c("out21",
                                        "nov21"))) %>%
    filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                             "Nov/21")) %>%
    group_by(Sample,
             Read_origin,
             Primer,
             #`Curated ID`,
             Year,
             Point,
             Expedition,
             Mês
    ) %>%
    summarize(`Num ASVs` = length(unique(`ASV (Sequence)`)),
              `Num OTUs` = length(unique(`OTU`)),
              `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
              `Ponto` = Point) %>%
    mutate(Sample = factor(Sample, levels = c("L1_out21",
                                              "L1_nov21",
                                              "L2_out21",
                                              "L2_nov21",
                                              "L3_out21",
                                              "L3_nov21",
                                              "L4_out21",
                                              "L4_nov21"
    ))) %>%
    ungroup() %>%
    unique()
  
  sample_tbl <- raw_results_tbl %>%
    filter(!`Curated ID` %in% c("Actinopteri", ## tirando as ASVs que nao foram identificadas a nivel de especie, nao-peixes e NA
                                "Astyanax",
                                "Characidae",
                                "Characidium",
                                "Characiformes",
                                "Cichla",
                                "Cichlidae",
                                "Hoplias",
                                "Pimelodus",
                                "",
                                "Cavia magna",
                                "Cutibacterium acnes",
                                "Bos taurus",
                                "Canis familiaris",
                                "Didelphis albiventris (Gamba)",
                                "Homo sapiens",
                                "Hydrochaeris hydrochaeris (Capivara)",
                                "Nannopterum brasilianus",
                                "Oryctolagus cuniculus (Coelho-bravo)",
                                "Progne chalybea (Andorinha-grande)",
                                "Sus scrofa"
    )) %>%
    mutate("Mês" = str_split_fixed(string = .$Sample, # criando uma nova coluna quebrando as infos da coluna Sample
                                   pattern = "_",
                                   n = 2)[,2]) %>%
    mutate(Mês = factor(Mês, levels = c("out21",
                                        "nov21"))) %>%
    filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                             "Nov/21")) %>%
    group_by(Sample,
             Read_origin,
             Primer,
             #`Curated ID`,
             Year,
             Point,
             Expedition,
             Mês
    ) %>%
    summarize(`Num ASVs` = length((`ASV (Sequence)`)),
              `Num OTUs` = length((`OTU`)),
              `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
              `Ponto` = Point) %>%
    mutate(Sample = factor(Sample, levels = c("L1_out21",
                                              "L1_nov21",
                                              "L2_out21",
                                              "L2_nov21",
                                              "L3_out21",
                                              "L3_nov21",
                                              "L4_out21",
                                              "L4_nov21"
    ))) %>%
    ungroup() %>%
    unique()
  
}

## Separando por ponto
{
  point_tbl <- raw_results_tbl %>%
    filter(!`Curated ID` %in% c("Actinopteri", ## tirando as ASVs que nao foram identificadas a nivel de especie, nao-peixes e NA
                                "Astyanax",
                                "Characidae",
                                "Characidium",
                                "Characiformes",
                                "Cichla",
                                "Cichlidae",
                                "Hoplias",
                                "Pimelodus",
                                "",
                                "Cavia magna",
                                "Cutibacterium acnes",
                                "Bos taurus",
                                "Canis familiaris",
                                "Didelphis albiventris (Gamba)",
                                "Homo sapiens",
                                "Hydrochaeris hydrochaeris (Capivara)",
                                "Nannopterum brasilianus",
                                "Oryctolagus cuniculus (Coelho-bravo)",
                                "Progne chalybea (Andorinha-grande)",
                                "Sus scrofa"
                                )) %>%
    filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                             "Nov/21")) %>%
    group_by(#Sample,
      #`Curated ID`,
      #Year,
      Point
      #,Expedition,
      #Mês
      ) %>%
    summarize(`Num ASVs` = length((`ASV (Sequence)`)),
              `ID Abundance on sample (%)` = sum(`Relative abundance on sample`))
  }

## Separando por expedicao
{
  expedition_tbl <- raw_results_tbl %>%
    filter(!`Curated ID` %in% c("Actinopteri", ## tirando as ASVs que nao foram identificadas a nivel de especie, nao-peixes e NA
                                "Astyanax",
                                "Characidae",
                                "Characidium",
                                "Characiformes",
                                "Cichla",
                                "Cichlidae",
                                "Hoplias",
                                "Pimelodus",
                                "",
                                "Cavia magna",
                                "Cutibacterium acnes",
                                "Bos taurus",
                                "Canis familiaris",
                                "Didelphis albiventris (Gamba)",
                                "Homo sapiens",
                                "Hydrochaeris hydrochaeris (Capivara)",
                                "Nannopterum brasilianus",
                                "Oryctolagus cuniculus (Coelho-bravo)",
                                "Progne chalybea (Andorinha-grande)",
                                "Sus scrofa"
    )) %>%
    filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                             "Nov/21")) %>%
    group_by(#Sample,
      #`Curated ID`,
      #Year,
      # Point
      Expedition
      #,Mês
    ) %>%
    summarize(`Num ASVs` = length((`ASV (Sequence)`)),
              `ID Abundance on sample (%)` = sum(`Relative abundance on sample`))
}

# Calculo media/sd de ASVs por amostra
mean(sample_tbl$`Num ASVs`) #34.875
sd(sample_tbl$`Num ASVs`) #16.39136

# Calculo media/sd de ASVs por ponto
mean(point_tbl$`Num ASVs`) #69.75
sd(point_tbl$`Num ASVs`) #5.5

# Calculo media/sd de ASVs por expedicao
mean(expedition_tbl$`Num ASVs`) #139.5
sd(expedition_tbl$`Num ASVs`) #4.94

# Separando spp por ponto
{
  spp_tbl <- alpha_tbl %>% 
    select(Sample) %>%
    table() %>%
    as_tibble()
  spp_tbl <- rename(spp_tbl, Sample = .) 
  spp_tbl <- rename(spp_tbl, n_spp = n)
  view(spp_tbl)
}

# Calculo media/sd de spp por amostra
mean(spp_tbl$n_spp) #11.875
sd(spp_tbl$n_spp) #5.488625



