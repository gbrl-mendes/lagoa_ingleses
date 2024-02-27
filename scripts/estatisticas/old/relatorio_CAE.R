"Relatório CAE"
"Mendes, GA; Hilário, OH"
"02/2023"

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
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v2023.csv"), sep = ",", check.names = FALSE) %>% tibble()
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Oreochromis niloticus")] <- "Tilapia rendalli" # Oreochromis niloticus é Tilapia
}

# Filtrando os dados ----

  ## Separando raw_results por ano
  
  raw_2020 <- raw_results_tbl %>% # Campanhas 2020
    filter(Year == 2020)
  raw_2021 <- raw_results_tbl %>% # Campanhas 2021
    filter(Year == 2021)
  
  ## Calculando total de ASVs apos o uso do DADA2 (sem filtros)
  
  raw_2020$Abundance %>% sum() #668.796 seqs em 2020
  raw_2021$Abundance %>% sum() #2.062.037 seqs em 2021
  
  ## Alpha tbl
  
    # Para 2020
    alpha_tbl_2020 <- raw_2020 %>%
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
      filter(Expedition %in% c("Nov_Dec/20",
                               "Nov/20",
                               "Dec/20"
      )) %>%
      group_by(Sample,
               Read_origin,
               Primer,
               `Curated ID`,
               Year,
               Point,
               Expedition,
               Class,
               Order,
               Family
      ) %>%
      summarize(`Num ASVs` = length(unique(`ASV (Sequence)`)),
                `Num OTUs` = length(unique(`OTU`)),
                `ID Abundance on sample (%)` = sum(`Relative abundance on sample`),
                `Ponto` = Point) %>%
      mutate(Sample = factor(Sample, levels = c("LI1-neo-mi",
                                                "L2_nov20",
                                                "L2_dez20"
                                                ))) %>%
      ungroup() %>%
      unique()

    ## Para 2021
    alpha_tbl_2021 <- raw_2021 %>%
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
      filter(Expedition %in% c("out/21",
                               "Nov/21"
      )) %>%
      group_by(Sample,
               Read_origin,
               Primer,
               `Curated ID`,
               Year,
               Point,
               Expedition,
               Class,
               Order,
               Family
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
    
  # Calculando o numero de ASVs geral
    #1 passo: comentar o trecho do codigo que filtra ASVs que nao sao de peixes
    #2 passo:
    alpha_tbl_2020$`Num ASVs` %>% sum()
    alpha_tbl_2021$`Num ASVs` %>% sum()
  
  # Calculando o numero de ASVs de peixes id a nivel de spp
    alpha_tbl_2020$`Num ASVs` %>% sum()
    alpha_tbl_2021$`Num ASVs` %>% sum()
  
  # Calculando o numero de spp
    alpha_tbl_2020$`Curated ID` %>% unique() %>% length()
    alpha_tbl_2021$`Curated ID` %>% unique() %>% length()
  
  # Recuperando o nome das spp
    alpha_tbl_2020$`Curated ID` %>% unique() %>% sort() %>% cat(sep = ", ")
    alpha_tbl_2021$`Curated ID` %>% unique() %>% sort() %>% cat(sep = ", ")
    
  # Spp por ponto em 2021
    L1_2021 <- alpha_tbl_2021 %>% filter(Point == "L1") %>% select("Curated ID")
    L2_2021 <- alpha_tbl_2021 %>% filter(Point == "L2") %>% select("Curated ID")
    L3_2021 <- alpha_tbl_2021 %>% filter(Point == "L3") %>% select("Curated ID")
    L4_2021 <- alpha_tbl_2021 %>% filter(Point == "L4") %>% select("Curated ID")
    
  # Spp exclusivas de cada ponto em 2021
    uniq_L1_2021 <- setdiff(as.character(unlist(L1_2021)), as.character(unlist((c(L2_2021, L3_2021, L4_2021))))) %>% sort() %>% cat(sep = ", ")
    uniq_L2_2021 <- setdiff(as.character(unlist(L2_2021)), as.character(unlist((c(L1_2021, L3_2021, L4_2021))))) %>% sort() %>% cat(sep = ", ")
    uniq_L3_2021 <- setdiff(as.character(unlist(L3_2021)), as.character(unlist((c(L1_2021, L2_2021, L4_2021))))) %>% sort() %>% cat(sep = ", ")
    uniq_L4_2021 <- setdiff(as.character(unlist(L4_2021)), as.character(unlist((c(L1_2021, L2_2021, L3_2021))))) %>% sort() %>% cat(sep = ", ")
    
  # Spp em comum com todos os pontos em 2021
    shared_spp_2021 <- intersect(as.character(unlist(L1_2021)), as.character(unlist((c(L2_2021, L3_2021, L4_2021)))))
    
  # Descobrir quais classes estao presentes
    # Em 2020
    alpha_tbl_2020$Class %>% unique() %>% sort() %>% cat(sep = ", ") # Actinopteri
    # Em 2021
    alpha_tbl_2021$Class %>% unique() %>% sort() %>% cat(sep = ", ") # Actinopteri
  # Descobrir quais ordens estao presentes
    # Em 2020
    alpha_tbl_2020$Order %>% unique() %>% sort() %>% cat(sep = ", ") # Characiformes, Cichliformes, Gymnotiformes, Siluriformes
    # Em 2021
    alpha_tbl_2021$Order %>% unique() %>% sort() %>% cat(sep = ", ") # Characiformes, Cichliformes, Cyprinodontiformes, Gymnotiformes, Salmoniformes, Siluriformes
  
      
    geom_tile() +
      # geom_text(aes(label= RRA),
      geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
                colour = "black", size = 3
      )
    
    
    