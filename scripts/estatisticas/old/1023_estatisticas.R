"Estatísticas"
"Mendes, GA; Hilário, OH"
"10/2023"

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

# Usando os objetos criados pelo codigo 102023_tile-plots.r

## Estatísticas das ASVs ----
{
  # Filtrando os dados ----
  
  ## Separando raw_results por ano
  
  raw_2020 <- raw_results_tbl %>% # Campanhas 2020
    filter(Year == 2020)
  raw_2021 <- raw_results_tbl %>% # Campanhas 2021
    filter(Year == 2021)
  
  ## Alpha tbl
  {
    # Para 2020
    alpha_tbl_2020 <- raw_2020 %>%
      filter(`Relative abundance on sample` >= 0.01) %>% # filtrar
      filter(!`Curated ID` %in% c("Actinopteri", ## tirando as ASVs que nao foram identificadas a nivel de especie, nao-peixes e NA
                                  "Astyanax",
                                  "Characidae",
                                  "Characidium",
                                  "Characiformes",
                                  "Cichla",
                                  "Cichlidae",
                                  "Hoplias",
                                  "Pimelodus",
                                  "Siluriformes",
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
      # mutate("Mês" = str_split_fixed(string = .$Sample, # criando uma nova coluna quebrando as infos da coluna Sample
      #                                pattern = "_",
      #                                n = 2)[,2]) %>%
      # mutate(Mês = factor(Mês, levels = c("dez20",
      #                                     "nov20",
      #                                     "out21",
      #                                     "nov21"))) %>%
      filter(Expedition %in% c("Nov_Dec/20",
                               "Nov/20",
                               "Dec/20",
                               "out/21",
                               "Nov/21"
                               )) %>%
      group_by(Sample,
               Read_origin,
               Primer,
               `Curated ID`,
               Year,
               Point,
               Expedition
      ) %>%
      summarize(`Num ASVs` = length(unique(`ASV (Sequence)`)),
                `Num OTUs` = length(unique(`OTU`)),
                `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
                `Ponto` = Point) %>%
      mutate(Sample = factor(Sample, levels = c("LI1-neo-mi",
                                                "L2_nov20",
                                                "L2_dez20",
                                                "L1_out21",
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
    
    ## Para 2021
    alpha_tbl_2021 <- raw_2021 %>%
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
                                  "Siluriformes",
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
      # mutate("Mês" = str_split_fixed(string = .$Sample, # criando uma nova coluna quebrando as infos da coluna Sample
      #                                pattern = "_",
      #                                n = 2)[,2]) %>%
      # mutate(Mês = factor(Mês, levels = c("dez20",
      #                                     "nov20",
      #                                     "out21",
      #                                     "nov21"))) %>%
      filter(Expedition %in% c("Nov_Dec/20",
                               "Nov/20",
                               "Dec/20",
                               "out/21",
                               "Nov/21"
      )) %>%
      group_by(Sample,
               Read_origin,
               Primer,
               `Curated ID`,
               Year,
               Point,
               Expedition
      ) %>%
      summarize(`Num ASVs` = length(unique(`ASV (Sequence)`)),
                `Num OTUs` = length(unique(`OTU`)),
                `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
                `Ponto` = Point) %>%
      mutate(Sample = factor(Sample, levels = c("LI1-neo-mi",
                                                "L2_nov20",
                                                "L2_dez20",
                                                "L1_out21",
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
                                  "Siluriformes",
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
                                  "Siluriformes",
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
                                  "Siluriformes",
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

  ## Separando por ano
  
  alpha_tbl_2020 <- alpha_tbl %>% # Campanhas 2020
    filter(Year == 2020)
  
  alpha_tbl_2021 <- alpha_tbl %>% # Campanhas 2021
    filter(Year == 2021)
  }  


## Estatísticas das Spp ----
{
  # Separando spp por amostra ----
  {
    spp_tbl <- alpha_tbl %>% 
      select(Sample) %>%
      table() %>%
      as_tibble()
    spp_tbl <- rename(spp_tbl, Sample = .) 
    spp_tbl <- rename(spp_tbl, n_spp = n)
    view(spp_tbl)
  }
  
  # Separando spp por ponto ----
  {
    spp_point_tbl <- alpha_tbl %>% 
      select(Point) %>%
      table() %>%
      as_tibble()
    spp_point_tbl <- rename(spp_point_tbl, Point = .) 
    spp_point_tbl <- rename(spp_point_tbl, n_spp = n)
    view(spp_point_tbl)
  }
  
  # Separando spp por expedicao ----
  {
    spp_exped_tbl <- alpha_tbl %>% 
      select(Expedition) %>%
      table() %>%
      as_tibble()
    spp_exped_tbl <- rename(spp_exped_tbl, Expedition = .) 
    spp_exped_tbl <- rename(spp_exped_tbl, n_spp = n)
    view(spp_exped_tbl)
  }
  
  # Calculo media/sd de spp por amostra
  mean(spp_tbl$n_spp) #11.875
  sd(spp_tbl$n_spp) #5.488625
  
  # Calculo media/sd de spp por ponto
  mean(spp_point_tbl$n_spp) #23.75
  sd(spp_point_tbl$n_spp) #6.130525
  
  # Calculo media/sd de spp por expedicao
  mean(spp_exped_tbl$n_spp) #47.5
  sd(spp_exped_tbl$n_spp) #7.778175
}
  
  
  ## Numeros para o relatorio do CAE ----
  
  # Obter numero de ASVs de cada campanha
  alpha_tbl_2020$Abundance %>% sum()
  alpha_tbl_2021$Abundance %>% sum()
  
  # Obter o numero de especies
  alpha_tbl_2020$`Curated ID` %>% unique() %>% sort() %>% length() # em 2020
  alpha_tbl_2021$`Curated ID` %>% unique() %>% sort() %>% length() # em 2021
  
  # Obter quais as especies de peixes que foram identificadas
  alpha_tbl_2020$`Curated ID` %>% unique() %>% sort() %>% cat(sep =", ")
  alpha_tbl_2021$`Curated ID` %>% unique() %>% sort() %>% cat(sep =", ")
    # ira apresentar o resultado como uma lista sem paragrafos, com os nomes 
    # separados por virgulas

  