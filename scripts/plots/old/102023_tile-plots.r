---
title: "Tile plots Lagoa dos Ingleses"
author: "Gabriel Mendes"
date: "10/2023"
---

## Carregando bibliotecas ----
{ 
  # .libPaths(c(.libPaths(),"/opt/R/4.1.0/lib/R/library"))
  library("Biostrings")
  library("DECIPHER")
  library("factoextra")
  library("future")
  library("ggh4x")
  library("ggpubr")
  library("Matrix")
  library("phyloseq")
  library("ShortRead")
  library("stringr")
  library("vegan")
  library("tidyverse")
  library("dplyr")
  library("readxl")
}

## Caminhos ----
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figuras")
  tbl_path <- paste0(prjct_path,"tabelas/raw/run_2_4_5")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

## Obtencao dos dados ----

pre_raw_results_tbl <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/BLAST_12sdb_nt/lagoa_ingleses-Complete_analysis_results-2023-11-24.xlsx") %>% tibble()
curated_ids_tbl <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/BLAST_12sdb_nt/curated_lagoa_ingleses-ASVs_x_amostras-2023-12-08.xlsx") %>% tibble() 
blast_tax <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/BLAST_12sdb_nt/curated_tax_blast.xlsx")
# curated_ids_tbl <- read_excel(paste0(tbl_path,"/", "lagoa_ingleses-ASVs_x_amostras-2023-10-10.xlsx")) %>% tibble()
  
## BLAST vs DADA ----

#08/12: BLAST NT & 12sLGCdb + DADA2
curated_ids_tbl_08 <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/blast_vs_dada/curated_lagoa_ingleses-ASVs_x_amostras-2023-12-08.xlsx") %>% tibble()
#13/11: BLAST NT + DADA2
curated_ids_tbl_13 <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/blast_vs_dada/curated_lagoa_ingleses-ASVs_x_amostras-2023-11-13.xlsx") %>% tibble()
    
#Complete_analysis
pre_raw_results_tbl_24 <- read_excel("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/blast_vs_dada/lagoa_ingleses-Complete_analysis_results-2023-11-24.xlsx") %>% tibble()
# como eu nao gerei a complete_analysis da planilha do dia 13/11, irei usar a mesma para as duas planilhas. 
# A diferenca e que eu vou desconsiderar as colunas Final ID (BLASTn),	BLASTn pseudo-score,	Class (DADA2)	e Order (DADA2)
# Essas colunas so dizem respeito a planilha do dia 24/11

  

## Alteracoes pos-pipeline ----
 
## Adicionando a raw_results_tbl os nomes que foram curados manualmente
  
## 1o selecionando apenas as colunas que importam 
  
curated_ids_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  
curated_ids_tbl <-
  curated_ids_tbl %>% 
  select(c(
    "ASV header",
    "ASV (Sequence)",
    "Curated ID",
    "Final ID (BLASTn)",
    "BLASTn pseudo-score",
    "Class (BLASTn)",
    "Curated Order"
    )) %>% 
  unique()
    
View(curated_ids_tbl)

pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

pre_raw_results_tbl <- 
  pre_raw_results_tbl %>% 
  select(c("Sample",
           "Unique_File_name",
           "OTU",
           "Abundance",
           "Sample total abundance",
           "Relative abundance on sample",
           "Relative abundance to all samples",
           "Read origin",
           "Primer expected length",
           "ASV Size (pb)",
           "ASV header"))
  
## 2o renomeando segundo os nomes curados e reordenando 
  
raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  
raw_results_tbl <- pre_raw_results_tbl %>%
  left_join(curated_ids_tbl,
            by = c("ASV header")) %>%
  select(c("ASV header", #oriundas da curated_ids_tbl
           "ASV (Sequence)",
           "Curated ID",
           "Final ID (BLASTn)",
           "BLASTn pseudo-score",
           "Class (BLASTn)",
           "Curated Order",
           "Sample", #oriundas da pre_raw_samples_tbl
           "Unique_File_name",
           "OTU",
           "Abundance",
           "Sample total abundance",
           "Relative abundance on sample",
           "Relative abundance to all samples",
           "Read origin",
           "Primer expected length",
           "ASV Size (pb)",
           "ASV header"
           ))

View(raw_results_tbl)
  
# verificando o total de reads por amostra

raw_results_tbl %>% group_by(Sample) %>% 
  summarise(total_abd = sum(Abundance))
  
# comparar com o total que saiu do pipeline para ver se nao ha perdas:
  
  # > smp_abd_ID_Final %>% group_by(Sample) %>% 
  #   +  summarise(total_abd = sum(Abundance))
  # # A tibble: 11 × 2
  # Sample                total_abd
  # <chr>                     <int>
  # 1 L1_nov_dec_20_mi        41539
  # 2 L1_nov21                66633
  # 3 L1_out21               383952
  # 4 L2_dez20               357293
  # 5 L2_nov20               207421
  # 6 L2_nov21               184235
  # 7 L2_out21               530847
  # 8 L3_nov21                49611
  # 9 L3_out21               369880
  # 10 L4_nov21              514634
  # 11 L4_out21                 123
  
# Verificando se ha diferencas entre a tabela de ASVs e a tabela de IDs curadas
dif_raw_pre <- pre_raw_results_tbl %>%
  anti_join(raw_results_tbl)
  
View(dif_raw_pre) #vazio e` bom!

# Apos ver que tudo que esta na tabela eh o que saiu do pipeline, podemos retirar 
# tudo o que estiver fora do intervalo do amplicon e identificacoes NA

raw_results_tbl <-
  raw_results_tbl %>% 
  filter(!`Primer expected length` %in% "out of range") %>% 
  filter(!`Curated ID` == "NA") %>% 
  filter(!(`Class (BLASTn)` %in% c("Mammalia", "Aves") & `BLASTn pseudo-score` < 98))

View(raw_results_tbl)
  
## 3o adicionando metadados que faltaram ----
  
raw_results_tbl$Sample %>% unique() %>% sort() %>% paste0(collapse = '",\n"') %>% cat()
  
  # sample
  Sample <- c("L1_nov_dec_20_mi", 
              "L2_nov20",
              "L2_dez20",
              "L1_out21",
              "L2_out21",
              "L3_out21",
              "L4_out21",
              "L1_nov21",
              "L2_nov21",
              "L3_nov21",
              "L4_nov21"
              )
  
  # expedition
  expedition <- c("Novembro e Dezembro 2020",
                  "Novembro 2020",
                  "Dezembro 2020",
                  "Outubro 2021",
                  "Outubro 2021",
                  "Outubro 2021",
                  "Outubro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021",
                  "Novembro 2021"
                  )
  
  # year
  year <- c("2020",
            "2020",
            "2020",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021",
            "2021"
            )
  
  # new_name
  new_name <- c("Fundação",
                "Ponte",
                "Ponte",
                "Prainha",
                "Barragem",
                "Ponte",
                "Fundação",
                "Prainha",
                "Barragem",
                "Ponte",
                "Fundação"
                )
  # filter
  filter <- c("MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE",
              "MCE"
              )
  # run
  run <- c("run2_ago21",
           "run4_out21",
           "run4_out21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21",
           "run5_dez21"
           )

  # as tibble
  metadata_tbl <- tibble(run, 
                         Sample, 
                         new_name, 
                         expedition, 
                         year) 
  
  # mergindo os metadados com raw_results_tbl
  
  raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  
  raw_curated_tbl <- left_join(raw_results_tbl, metadata_tbl,
                      by = "Sample") %>% 
    relocate(c("Sample",
             "run",
             "new_name",
             "expedition",
             "year",
             "Curated ID",
             "Curated Order",
             # "Obs. Curadoria",
             "Final ID (BLASTn)",
             "BLASTn pseudo-score",
             "Class (BLASTn)",
             # "Order (DADA2)",
             # "Curated Max. taxonomy",
             # "Primer",
             "OTU",
             "Abundance",
             "Sample total abundance",
             "Relative abundance on sample",
             "Relative abundance to all samples",
             "Read origin",
             "Primer expected length",
             "ASV Size (pb)",
             "ASV header",
             "ASV (Sequence)"
             ))
  
  View(raw_curated_tbl)
  
  diff_cur_raw <- dplyr::anti_join(#curated_ids_tbl,
                                   #raw_results_tbl,
                                   raw_curated_tbl,
                                   # raw_results_tbl,
                                   #raw_curated_tbl,
                                   curated_ids_tbl,
                                   # ver se ocorreu sem problemas o left_join
                                    by = "ASV header")
  View(diff_cur_raw) #vazio e` bom!
  


## Tabelas ----

# Criacao da lista com os possiveis nomes atribuidos as ASVs

raw_curated_tbl %>% colnames()
raw_curated_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

# Agrupamento da ASVs que possuem os mesmos atributos abaixo

grouped_by_ID_BLASTid <- raw_curated_tbl %>%
  relocate(c("Sample",
             "expedition",
             "year",
             "Final ID (BLASTn)",
             "BLASTn pseudo-score",
             "Curated ID",
             "OTU",
             "Curated Order",
             "Relative abundance on sample",
             "new_name",
             "Class (BLASTn)",
             # "Order (DADA2)",
             "Curated Order",
             # "Curated Max. taxonomy",
             "Abundance",
             "ASV header",
             "ASV (Sequence)",
             "Primer expected length"
             )) %>%
  group_by(Sample, `Curated ID`,`Final ID (BLASTn)`,`BLASTn pseudo-score`, expedition, new_name, year, `ASV header`) %>%
  summarize("RRA" = sum(`Relative abundance on sample`),
            "Class (BLASTn)" = unique(`Class (BLASTn)`),
            "OTU" = unique(OTU),
            "Curated Order" = unique(`Curated Order`),
            # "Order (DADA2)" = unique(`Order (DADA2)`),
            # "Curated Max. taxonomy" = unique(`Curated Max. taxonomy`),
            "Abundance" = sum(Abundance),  # Agregando a coluna Abundance
            "ASV (Sequence)" = unique(`ASV (Sequence)`),
            "Primer expected length" = unique(`Primer expected length`)
            ) %>%
  mutate('Nivel' = ifelse(`year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                              "Cheio", "Vazio")) %>%
  ungroup()
    
View(grouped_by_ID_BLASTid)

# verificando o total de reads por amostra

grouped_by_ID_BLASTid %>% group_by(Sample) %>% 
  summarise(total_abd = sum(Abundance))

# comparar com o total que saiu do pipeline para ver se nao ha perdas:

# > smp_abd_ID_Final %>% group_by(Sample) %>% 
#   +  summarise(total_abd = sum(Abundance))
# # A tibble: 11 × 2
# Sample                total_abd
# <chr>                     <int>
# 1 L1_nov_dec_20_mi        41539
# 2 L1_nov21                66633
# 3 L1_out21               383952
# 4 L2_dez20               357293
# 5 L2_nov20               207421
# 6 L2_nov21               184235
# 7 L2_out21               530847
# 8 L3_nov21                49611
# 9 L3_out21               369880
# 10 L4_nov21              514634
# 11 L4_out21                 123

# adicionando a taxonomia do BLAST

blast_tax <- blast_tax %>% 
  rename(`ASV (Sequence)` = ASV)

blast_tax_less <- blast_tax %>% 
  select(c("ASV (Sequence)",
           "Order (BLASTn)",
           "Family (BLASTn)",
           "Genus (BLASTn)"
           ))


grouped_by_ID_BLASTid <- grouped_by_ID_BLASTid %>% 
  left_join(blast_tax_less,
            by = "ASV (Sequence)") %>% 
  mutate("Curated genus" = str_split_fixed(string = .$`Curated ID`, 
                                 pattern = " ",
                                 n = 2)[,1]) %>% 
  select(c("Sample",
           "Curated ID", 
           "Final ID (BLASTn)",
           "BLASTn pseudo-score",
           "ASV (Sequence)",
           "Class (BLASTn)",
           "Order (BLASTn)",
           "Curated Order",
           "Family (BLASTn)",
           "Genus (BLASTn)",
           "Curated genus",
           "expedition",
           "new_name",
           "year", 
           "ASV header",
           "RRA",
           "OTU", 
           "Abundance",
           "ASV (Sequence)",
           "Primer expected length",
           "Nivel" 
  ))

View(grouped_by_ID_BLASTid)

## Organizar as especies ----

grouped_by_ID_BLASTid$`Curated ID` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
grouped_by_ID_BLASTid$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  
# Organizar ordem

subset(grouped_by_ID_BLASTid, `Class (BLASTn)` == "Actinopteri")$`Curated Order` %>%
  sort() %>% unique() # ver ordens
    
fish_ordens <- c("Characiformes", # organizar as ordem das ordens de peixes
                 "Cichliformes",
                 "Gymnotiformes",
                 "Salmoniformes",
                 "Siluriformes")
    
subset(grouped_by_ID_BLASTid, `Class (BLASTn)` != "Actinopteri")$`Class (BLASTn)` %>%
  sort() %>% unique() # ver classes

nfish_classes <- c("Aves", # ordem das classes de nao peixes
                   "Mammalia")
  
  
# Organizar classes
  
subset(grouped_by_ID_BLASTid, `Class (BLASTn)` != "Actinopteri")$`Curated Order` %>%
  sort() %>% unique() # ver ordens
    
nfish_ordens <- c("Artiodactyla", # ordem das ordens de nao peixes
                  "Carnivora",
                  "Didelphimorphia",
                  "Lagomorpha",
                  "Passeriformes",
                  "Primates",
                  "Rodentia",
                  "Suliformes")
  
  
# Definir quais serao as especies e ordenar
  
    spp_total <- raw_curated_tbl %>%
      filter(Sample %in% c("L1_nov_dec_20_mi",
                           "L2_nov20",
                           "L2_dez20")) %>% 
      # filter(`Class (BLASTn)` == "Actinopteri") %>%
      # filter(`Class (BLASTn)` == "Mammalia") %>%
      # mutate(words_count = str_count(`Curated ID`, "\\w+")) %>%
      # filter(words_count >= 2) %>% #so a nivel de spp
      # filter(str_detect(`Curated ID`, "sp.")) %>%
      select(`Curated ID`) %>%
      unique() %>%
      arrange(`Curated ID`) %>%
      pull()
      
    list(spp_total)
    
    raw_curated_tbl %>% 
      
    
  
    {
    spp_fish <- raw_curated_tbl %>%  # especies de peixes
      filter(`Primer expected length` == "in range") %>% 
      filter(`Curated ID` != c("NA","")) %>% 
      filter(`Class (BLASTn)` == "Actinopteri") %>%
      mutate(words_count = str_count(`Curated ID`, "\\w+")) %>%
      filter(words_count >= 2) %>% #so a nivel de spp
      select(`Curated ID`) %>%
      unique() %>% 
      arrange(`Curated ID`) %>%
      pull()
    list(spp_fish)
    }
  
  {
    no_fish <- raw_curated_tbl %>%  # outras especies
      filter(`Class (BLASTn)` != "Actinopteri") %>% 
      # filter(`Curated ID` != "") %>% 
      select(`Curated ID`) %>% 
      unique() %>% 
      arrange(`Curated ID`) %>% 
      pull()
    list(no_fish)
  }  
  
  # Definir quais serao as amostras e ordenar
  {
    raw_curated_tbl$Sample %>% unique()
    
    samples <- c("L1_nov_dec_20_mi",
                 "L2_nov20",
                 "L2_dez20",
                 "L1_out21",
                 "L2_out21",
                 "L3_out21",
                 "L4_out21",
                 "L1_nov21",
                 "L2_nov21",
                 "L3_nov21",
                 "L4_nov21"
                 )
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
      filter(`Class (BLASTn)` == "Actinopteri") %>%
      # mutate(`Curated ID` = factor(`Curated ID`,
      #                              levels = rev(spp_fish))) %>%
      filter(`Curated ID` %in% spp_fish) %>% 
      # filter(`Curated ID` != "") %>% #tirei la em cima
      # filter(!is.na(`Curated ID`)) %>% #tirei la em cima
      arrange(`Curated ID`) %>% 
      mutate('Nivel' = ifelse(`year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                              "Cheio", "Vazio")) %>%
      # correcao do RRA ja que as ASVs sem ID e de nao-peixes foram removidas alterando a abd total
      group_by(Sample) %>%
      mutate(RRA = RRA / sum(RRA)) %>%
      ungroup()
  
    View(fish_ID_tbl)  
  
    # conferir se a correcao do RRA funcionou
    fish_ID_tbl %>% 
      group_by(Sample) %>% 
      summarize(total_RRA = sum(RRA))
    }
  
  # tabela nao peixes
  {
    nfish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids nao peixes a nivel de spp
      filter(`Class (BLASTn)` == c("Aves", "Mammalia")) %>% 
      mutate(`Order (DADA2)` = factor(`Order (DADA2)`,
                              levels = nfish_ordens)) %>% 
      mutate('Nivel' = ifelse(`year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                              "Cheio", "Vazio")) %>%
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

  
## Raw statistics ----
    
    raw_results_tbl %>% 
      filter(Sample %in% c("L1_nov_dec_20_mi",
                           "L2_nov20",
                           "L2_dez20")) %>% 
      group_by(`Class (BLASTn)`
               # ,Sample
               ) %>% 
      summarize(
        "Total reads" = sum(Abundance),
        "ASVs" = length(unique(`ASV header`)),
        "OTUs" = length(unique(OTU)),
        "Species" = length(unique(`Curated ID`))
        # ,
        # "Nspecies" = length(unique(`Curated ID`)) - Species
        ) %>%  
      ungroup()
    
    
    # grouped by sample
    raw_stat_tbl <- grouped_by_ID_BLASTid %>%
      
      # teste para provar que com os filtros abaixo aplicados, os valores de ASVs, OTUs e spp 
      # irao ficar iguais aos da analise
      
      # mutate(`Curated ID` = factor(`Curated ID`,
      #                              levels = rev(spp_fish))) %>% 
      # mutate(`Order (DADA2)` = factor(`Order (DADA2)`,
      #                                 levels = fish_ordens)) %>%
      # filter(`Curated ID` != "") %>%
      # filter(!is.na(`Curated ID`)) %>%
      # filter(RRA >= 0.01) %>% 
      group_by(Sample,
               year,
               expedition,
               # Nivel
               ) %>%
      summarize(
        "Total reads" = sum(Abundance),
        "Month" = expedition,
        "Pond level" = Nivel,
        "Year" = year,
        "ASVs" = length(unique(`ASV header`)),
        "OTUs" = length(unique(OTU)),
        "Species" = length(unique(`Curated ID`)), #nao sao species, sao identificacoes unicas, podendo ser ordem, familia, genero ou spp
        "Sampling points" = new_name
        ) %>% 
      arrange(match(Sample, samples)) %>% 
      ungroup() %>% 
      select(c("Sample",
               "Sampling points",
               "Year",
               "Month",
               "Pond level",
               "Total reads",
               "ASVs",
               "OTUs",
               "Species")) %>% 
      unique()
    View(raw_stat_tbl)

    
## Curated statistics ----
  
    # grouped by sample
  stat_tbl <- fish_ID_tbl %>%
    filter(RRA >= 0.01) %>%
    group_by(Sample,
             year,
             expedition
             # ,
             # Nivel
             ) %>%
    summarize(
              "Total reads" = sum(Abundance),
              "Month" = expedition,
              "Pond level" = Nivel,
              "Year" = year,
              "ASVs" = length(unique(`ASV header`)),
              "OTUs" = length(unique(OTU)),
              "Species" = length(unique(`Curated ID`)),
              "Sampling points" = new_name
              ) %>% 
    ungroup() %>% 
    select(c("Sample",
             "Sampling points",
             "Year",
             "Month",
             "Pond level",
             "Total reads",
             "ASVs",
             "OTUs",
             "Species")) %>% 
    unique()
  View(stat_tbl)
  
  # total
stat_tbl2 <- fish_ID_tbl %>%
  # stat_tbl2 <- grouped_by_ID_BLASTid %>%
    filter(RRA >= 0.01) %>%
    # group_by(Sample,
    #          year,
    #          expedition
    #          # ,
    #          # Nivel
    #          ) %>%
    summarize(
              "Total reads" = sum(Abundance),
              # "Month" = expedition,
              # "Pond level" = Nivel,
              # "Year" = year,
              "ASVs" = length(unique(`ASV header`)),
              "OTUs" = length(unique(OTU)),
              "Species" = length(unique(`Curated ID`)),
              # "Sampling points" = new_name
              ) %>%
    ungroup() %>%
    select(c(
      # "Sample",
             # "Sampling points",
             # "Year",
             # "Month",
             # "Pond level",
             "Total reads",
             "ASVs",
             "OTUs",
             "Species")) %>%
    unique()
  View(stat_tbl2)
  
  #grouped by ASV
stat_ASV <- fish_ID_tbl %>%
    filter(RRA >= 0.01) %>%
    group_by(`ASV header`
             ) %>%
    summarize(`Curated ID`,
              "Total reads" = sum(Abundance),
              "ASV (Sequence)" = `ASV (Sequence)`
              # ,
              # "ASVs" = length(unique(`ASV header`)),
              # "OTUs" = length(unique(OTU)),
              # "Species" = length(unique(`Curated ID`)),
              ) %>%
    ungroup() %>%
    select(c(`Curated ID`,
             `ASV header`,
             `ASV (Sequence)`,
             "Total reads"
             # ,
             # "ASVs",
             # "OTUs",
             # "Species"
             )) %>%
    unique()
  View(stat_ASV)

## Tile Plots ----

# Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5

## vazio versus cheio (versao SBG)
  
  {
    cheio_vs_vazio <- fish_ID_tbl %>%
      filter(RRA >= 0.01) %>%
      mutate(expedition = factor(expedition, levels = expeditions)) %>% 
      mutate(Sample = factor(Sample, levels = samples)) %>% 
      mutate(new_name = factor(new_name)) %>% 
      # mutate(sample = factor(sample)) %>% 
      mutate(`Order (DADA2)` = factor(`Order (DADA2)`)) %>%
      group_by(expedition, `Curated ID`, new_name, `Order (DADA2)`, Nivel) %>%   # Agrupa por expedição, Curated ID e Point
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
                 rows = vars(`Order (DADA2)`),
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
  
  
## Data resume paper ----

  # peixes
  dt_fish_resume <- fish_ID_tbl %>%
    filter(RRA >= 0.01) %>%
    group_by(Nivel) %>%
    mutate("Abd total ano" = sum(Abundance),
           "RRA no ano" = Abundance/`Abd total ano`) %>%
    ungroup() %>%
    group_by(`Curated ID`, Nivel) %>%
    mutate(
      "Samples" = paste(unique(Sample), collapse = ",")
    ) %>%
    reframe(
      "Specie" = unique(`Curated ID`),
      "ASVs" = length(unique(`ASV header`)),
      "OTUs" = length(unique(OTU)),
      "RRA" = sum(`RRA no ano`),
      "Samples" = Samples
    ) %>% 
    mutate(
      "RRA_formatado" = format(RRA, scientific = FALSE)
      ) %>% 
    unique()
  View(dt_fish_resume)
  
  
  #peixes e nao-peixes  
  dt_all_resume <- grouped_by_ID_BLASTid %>%
    group_by(Nivel, new_name) %>%
    mutate("Abd total periodo" = sum(Abundance)) %>%
    ungroup() %>%
    group_by(`Curated ID`, 
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
      "Class" = `Class (BLASTn)`
    ) %>% 
    # select("Specie",
    #        "RRA",
    #        "ASVs",
    #        "OTUs",
    #        "Class",
    #        "Nivel") %>%
    mutate(
      "RRA_formatado" = format(RRA, scientific = TRUE)
    ) %>%
    unique()
  View(dt_all_resume)
  
  options(scipen = 99,
          digits = 3)
  
  #teste de RRA para ver se foi calculado corretamente
  dt_all_resume %>%
    filter(Nivel == "Cheio") %>%
    # filter(Sample == "L2_nov20") %>%
    filter(new_name == "Fundação") %>%
    pull(RRA) %>%
    sum()  

  wider_dt_all_resume <- dt_all_resume %>% 
    select(-c("RRA_formatado")) %>%
    mutate(RRA = round(RRA,digits = 4)) %>% 
    ungroup() %>% 
    unite(new_name, Nivel,col= "ponto_nivel") %>% 
    pivot_wider(id_cols = c("Class","Curated ID"),
                names_from = ponto_nivel,
                values_from =  c("RRA","ASVs","OTUs","Abundance"),
                names_glue = "{ponto_nivel}_{.value}") %>%
    select(sort(colnames(.))) %>% 
    relocate( "Class","Curated ID" ) 
    View(wider_dt_all_resume)
  
  write.csv(wider_dt_all_resume, "~/projetos/lagoa_ingleses/results/figuras/2023/defesa/wider_dt_all_resume.csv")
  
    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  