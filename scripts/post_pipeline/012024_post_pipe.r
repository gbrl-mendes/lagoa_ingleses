---
title: "Post-processing of ASVs"
author: "Gabriel Mendes"
date: "01/2024"
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
  tbl_path <- paste0(prjct_path,"tabelas")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

# ## Alteracoes pos-pipeline ----
# 
# ## Obtencao dos dados
# {
#   pre_complete_results_tbl <- read_excel(paste0(tbl_path,"/", "pre_curated/pre_curated-Complete_analysis_results-2024-01-10.xlsx")) %>% tibble()
#   pre_curated_ids_tbl <- read_excel(paste0(tbl_path,"/", "pre_curated/pre_curated_lagoa_ingleses-ASVs_x_amostras-2024-01-09.xlsx")) %>% tibble() 
#   blast_tax <- read.csv("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/tax_blast.csv", sep = ",", check.names = FALSE) 
# }
#  
# ## Adicionando a raw_results_tbl os nomes que foram curados manualmente
#   
# ## 1o selecionando apenas as colunas que importam 
#   
# pre_curated_ids_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
#   
# # Pre-Curated IDs
# pre_curated_ids_tbl <-
#   pre_curated_ids_tbl %>% 
#   select(c(
#     "ASV header",
#     "ASV (Sequence)",
#     "Curated ID",
#     "Final ID (BLASTn)",
#     "BLASTn pseudo-score",
#     "Class (BLASTn)",
#     # "Obs. Curadoria",
#     # "Possible contamination",
#     # "Curated Order"
#     )) %>% 
#   unique()
#     
# View(pre_curated_ids_tbled_ids_tbl)
# 
# # Pre-Complete analysis
# pre_complete_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
# 
# pre_complete_results_tbl <- 
#   pre_complete_results_tbl %>% 
#   select(c("Sample",
#            "Unique_File_name",
#            "OTU",
#            "Abundance",
#            "Sample total abundance",
#            "Relative abundance on sample",
#            "Relative abundance to all samples",
#            "Obs. Curadoria",
#            "Possible contamination",
#            "Read origin",
#            "Primer expected length",
#            "ASV Size (pb)",
#            "ASV header")) %>% 
#   mutate(Sample = str_replace(Sample, "__", "_"))
#   
# ## 2o renomeando segundo os nomes curados e reordenando 
#   
# complete_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
#   
# complete_results_tbl <- pre_complete_results_tbl %>%
#   left_join(pre_curated_ids_tbl,
#             by = c("ASV header")) %>%
#   select(c("ASV header", #oriundas da pre_curated_ids_tbl
#            "ASV (Sequence)",
#            "Curated ID",
#            "Final ID (BLASTn)",
#            "BLASTn pseudo-score",
#            "Class (BLASTn)",
#            "Obs. Curadoria",
#            "Possible contamination",
#            # "Curated Order",
#            "Sample", #oriundas da pre_complete_results_tbl
#            "Unique_File_name",
#            "OTU",
#            "Abundance",
#            "Sample total abundance",
#            "Relative abundance on sample",
#            "Relative abundance to all samples",
#            "Read origin",
#            "Primer expected length",
#            "ASV Size (pb)",
#            "ASV header"
#            ))
# 
# View(complete_results_tbl)
#   
# # verificando o total de reads por amostra
# 
# raw_reads <- complete_results_tbl %>% group_by(Sample) %>% 
#   summarise(total_abd = sum(Abundance))
#   
# View(raw_reads)
# 
# # comparar com o total que saiu do pipeline para ver se nao ha perdas:
#   
#   # > raw_reads <- smp_abd_ID_Final %>% group_by(Sample) %>% summarise(total_abd = sum(Abundance))
#   # > print(raw_reads, n= 28)
#   # A tibble: 28 × 2
#   # Sample              total_abd
#   # <chr>                   <int>
#   # 1 EM113_NEGPCR1             1
#   # 2 EM135c4c5_NEGPCR1        87
#   # 3 EM135c4c5_NEGPCR2        35
#   # 4 EM149_NEGPCR1          3110
#   # 5 EM149_NEGPCR2           255
#   # 6 EM156_NegPCR1          1786
#   # 7 L1_jan22             188961
#   # 8 L1_nov21              66411
#   # 9 L1_nov_dec_20_mi      41539
#   # 10 L1_out21            384243
#   # 11 L2_dez20            357293
#   # 12 L2_jan22            249943
#   # 13 L2_nov20            207421
#   # 14 L2_nov21            184202
#   # 15 L2_out21            530846
#   # 16 L3_jan22                37
#   # 17 L3_nov21             49612
#   # 18 L3_out21            371205
#   # 19 L4_jan22                32
#   # 20 L4_nov21            515845
#   # 21 L4_out21               123
#   # 22 S724_NEGPCR1            31
#   # 23 STX__L1_nov21       304589
#   # 24 STX__L2_nov21       313557
#   # 25 STX__L3_nov21       278340
#   # 26 STX__L4_nov21         1265
#   # 27 br_jan_22               44
#   # 28 br_nov21                39
#   
# # Verificando se ha diferencas entre a tabela de ASVs e a tabela de IDs curadas
# dif_raw_pre <- pre_complete_results_tbl %>%
#   anti_join(complete_results_tbl)
#   
# View(dif_raw_pre) #vazio e` bom!
# 
# # Apos ver que tudo que esta na tabela eh o que saiu do pipeline, podemos retirar 
# # tudo o que estiver fora do intervalo do amplicon, identificacoes NA, brancos e
# # possiveis contaminacoes
# 
# complete_results_tbl_filt <-
#   complete_results_tbl %>% 
#   filter(!`Primer expected length` %in% "out of range") %>% 
#   filter(!`Curated ID` == "NA") %>% 
#   filter(!(`Class (BLASTn)` %in% c("Mammalia", "Aves") & `BLASTn pseudo-score` < 98)) %>% 
#   filter(!(Sample %in% c("EM113_NEGPCR1",
#                          "EM135c4c5_NEGPCR1",
#                          "EM135c4c5_NEGPCR2",
#                          "EM149_NEGPCR1",
#                          "EM149_NEGPCR2",
#                          "EM156_NegPCR1",
#                          "S724_NEGPCR1",
#                          "br_jan_22",
#                          "br_nov21"))) %>% 
#   filter(!(`Possible contamination` %in% c("Possible contamination"))) %>% 
#   filter(!(`Obs. Curadoria` %in% c("Possible contamination")))
# 
# View(complete_results_tbl_filt)
# 
# # verificando o total de reads por amostra apos filtragem
# 
# filt_raw_reads <- complete_results_tbl_filt %>% group_by(Sample) %>% 
#   summarise(total_abd_filt = sum(Abundance))
# 
# View(filt_raw_reads)
# 
# stats_filt_reads <- raw_reads %>% 
#   left_join(filt_raw_reads, 
#             by = c("Sample")) %>% 
#   mutate("Proportion left" = total_abd_filt / total_abd * 100) 
# %>%
# filter(!is.na(total_abd_filt))
# 
# View(stats_filt_reads)
# 
# ## 3o adicionando metadados que faltaram ----
#   
# complete_results_tbl_filt$Sample %>% unique() %>% sort() %>% paste0(collapse = '",\n"') %>% cat()
# {
#   # sample
#   Sample <- c("L1_nov_dec_20_mi", 
#               "L2_nov20",
#               "L2_dez20",
#               "L1_out21",
#               "L2_out21",
#               "L3_out21",
#               "L4_out21",
#               "L1_nov21",
#               "L2_nov21",
#               "L3_nov21",
#               "L4_nov21",
#               "STX_L1_nov21",
#               "STX_L2_nov21",
#               "STX_L3_nov21",
#               "STX_L4_nov21",
#               "L1_jan22",
#               "L2_jan22",
#               "L3_jan22",
#               "L4_jan22"
#   )
#   
#   # expedition
#   expedition <- c("Novembro e Dezembro 2020",
#                   "Novembro 2020",
#                   "Dezembro 2020",
#                   "Outubro 2021",
#                   "Outubro 2021",
#                   "Outubro 2021",
#                   "Outubro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Novembro 2021",
#                   "Janeiro 2022",
#                   "Janeiro 2022",
#                   "Janeiro 2022",
#                   "Janeiro 2022"
#   )
#   
#   # year
#   year <- c("2020",
#             "2020",
#             "2020",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2021",
#             "2022",
#             "2022",
#             "2022",
#             "2022"
#   )
#   
#   # new_name
#   new_name <- c("Fundação",
#                 "Ponte",
#                 "Ponte",
#                 "Prainha",
#                 "Barragem",
#                 "Ponte",
#                 "Fundação",
#                 "Prainha",
#                 "Barragem",
#                 "Ponte",
#                 "Fundação",
#                 "Prainha",
#                 "Barragem",
#                 "Ponte",
#                 "Fundação",
#                 "Prainha",
#                 "Barragem",
#                 "Ponte",
#                 "Fundação"
#   )
#   # filter
#   filter <- c("MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE",
#               "Sterivex",
#               "Sterivex",
#               "Sterivex",
#               "Sterivex",
#               "MCE",
#               "MCE",
#               "MCE",
#               "MCE"
#   )
#   # run
#   run <- c("run2_ago21",
#            "run4_out21",
#            "run4_out21",
#            "run5_dez21",
#            "run5_dez21",
#            "run5_dez21",
#            "run5_dez21",
#            "run5_dez21",
#            "run5_dez21",
#            "run5_dez21",
#            "run5_dez21",
#            "EM156",
#            "EM156",
#            "EM156",
#            "EM156",
#            "EM156",
#            "EM156",
#            "EM156",
#            "EM156"
#   )
# }
#   
#   # as tibble
#   metadata_tbl <- tibble(run, 
#                          Sample,
#                          filter,
#                          new_name, 
#                          expedition, 
#                          year) 
#   
#   # mergindo os metadados com complete_results_tbl
#   
#   complete_results_tbl_filt %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
#   
#   complete_results_tbl_filt$Sample %>% unique()
#   
#   pre_raw_curated_tbl <- left_join(complete_results_tbl_filt, metadata_tbl,
#                       by = "Sample") %>% 
#     select(c("Sample",
#              "run",
#              "new_name",
#              "expedition",
#              "year",
#              "filter",
#              "Curated ID",
#              # "Curated Order",
#              # "Obs. Curadoria",
#              "Final ID (BLASTn)",
#              "BLASTn pseudo-score",
#              "Class (BLASTn)",
#              # "Order (DADA2)",
#              # "Curated Max. taxonomy",
#              # "Primer",
#              "OTU",
#              "Abundance",
#              "Sample total abundance",
#              "Relative abundance on sample",
#              "Relative abundance to all samples",
#              "Read origin",
#              "Primer expected length",
#              "ASV Size (pb)",
#              "ASV header",
#              "ASV (Sequence)"
#              ))
#   
#   View(pre_raw_curated_tbl)
#   
#   diff_cur_raw <- dplyr::anti_join(#curated_ids_tbl,
#                                    #complete_results_tbl,
#                                    pre_raw_curated_tbl,
#                                    # complete_results_tbl,
#                                    #raw_curated_tbl,
#                                    complete_results_tbl_filt,
#                                    # ver se ocorreu sem problemas o left_join
#                                     by = "ASV header")
#   View(diff_cur_raw) #vazio e` bom!
# 
# ## Tabelas ----
# 
# # Criacao da lista com os possiveis nomes atribuidos as ASVs
#   
#   pre_raw_curated_tbl %>% colnames()
#   pre_raw_curated_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
# 
# # Agrupamento da ASVs que possuem os mesmos atributos abaixo
# 
# pre_grouped_by_ID_BLASTid <- pre_raw_curated_tbl %>%
#   relocate(c("Sample",
#              "expedition",
#              "year",
#              "Final ID (BLASTn)",
#              "BLASTn pseudo-score",
#              "Curated ID",
#              "OTU",
#              # "Curated Order",
#              "Relative abundance on sample",
#              "new_name",
#              "Class (BLASTn)",
#              # "Order (DADA2)",
#              # "Curated Order",
#              # "Curated Max. taxonomy",
#              "Abundance",
#              "ASV header",
#              "ASV (Sequence)",
#              "Primer expected length"
#              )) %>%
#   group_by(Sample, `Curated ID`,`Final ID (BLASTn)`,`BLASTn pseudo-score`, expedition, new_name, year, `ASV header`) %>%
#   summarize("RRA" = sum(`Relative abundance on sample`),
#             "Class (BLASTn)" = unique(`Class (BLASTn)`),
#             "OTU" = unique(OTU),
#             # "Curated Order" = unique(`Curated Order`),
#             # "Order (DADA2)" = unique(`Order (DADA2)`),
#             # "Curated Max. taxonomy" = unique(`Curated Max. taxonomy`),
#             "Abundance" = sum(Abundance),  # Agregando a coluna Abundance
#             "ASV (Sequence)" = unique(`ASV (Sequence)`),
#             "Primer expected length" = unique(`Primer expected length`)
#             ) %>%
#   mutate(Nivel = case_when(`year` == 2020 ~ "Cheio", 
#                            `year` == 2022 ~ "Cheio", 
#                            TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o nivel da lagoa
#   ungroup()
#     
# View(pre_grouped_by_ID_BLASTid)
# 
# # verificando o total de reads por amostra
# 
# grouped_reads <- pre_grouped_by_ID_BLASTid %>% group_by(Sample) %>% 
#   summarise(total_abd = sum(Abundance))
# 
# View(grouped_reads)
# 
# # comparar com o total que saiu do pipeline para ver se nao ha perdas:
# 
# # > raw_reads < smp_abd_ID_Final %>% group_by(Sample) %>% 
# #   +  summarise(total_abd = sum(Abundance))
# # > print(raw_reads, n= 21)
# 
# # A tibble: 21 × 2
# #   Sample             total_abd
# #   <chr>                  <dbl>
# # 1 L1_jan22              188658
# # 2 L1_nov21               66426
# # 3 L1_nov_dec_20_mi       41539
# # 4 L1_out21              384154
# # 5 L2_dez20              357293
# # 6 L2_jan22              249320
# # 7 L2_nov20              207421
# # 8 L2_nov21              184214
# # 9 L2_out21              530824
# # 10 L3_jan22                 37
# # 11 L3_nov21              49612
# # 12 L3_out21             370515
# # 13 L4_jan22                 32
# # 14 L4_nov21             514706
# # 15 L4_out21                123
# # 16 STX__L1_nov21        303077
# # 17 STX__L2_nov21        312437
# # 18 STX__L3_nov21        278561
# # 19 STX__L4_nov21          1264
# # 20 br_jan_22                44
# # 21 br_nov21                 37
# 
# # adicionando a taxonomia do BLAST
# 
# blast_tax <- blast_tax %>% 
#   rename(`ASV (Sequence)` = "ASV")
# 
# blast_tax_less <- blast_tax[2:86] %>%
#   rename("ASV (Sequence)" = ASV) %>% 
#   select(c("ASV (Sequence)",
#            "Order (BLASTn)",
#            "Family (BLASTn)",
#            "Genus (BLASTn)"
#            )) 
#   
# pre_grouped_by_ID_BLASTid <- pre_grouped_by_ID_BLASTid %>% 
#   left_join(blast_tax_less,
#             by = "ASV (Sequence)") %>% 
#   mutate("Curated genus" = str_split_fixed(string = .$`Curated ID`, 
#                                  pattern = " ",
#                                  n = 2)[,1]) %>% 
#   select(c("Sample",
#            "Curated ID", 
#            "Final ID (BLASTn)",
#            "BLASTn pseudo-score",
#            "ASV (Sequence)",
#            "Class (BLASTn)",
#            "Order (BLASTn)",
#            # "Curated Order",
#            "Family (BLASTn)",
#            "Genus (BLASTn)",
#            "Curated genus",
#            "expedition",
#            "new_name",
#            "year", 
#            "ASV header",
#            "RRA",
#            "OTU", 
#            "Abundance",
#            "ASV (Sequence)",
#            "Primer expected length",
#            "Nivel" 
#   ))
# 
# View(pre_grouped_by_ID_BLASTid)
# 


## Alteracoes post-curated IDs ----

## Obtencao dos dados ----

pre_raw_results_tbl <- read_excel(paste0(tbl_path,"/", "curated/curated-Complete_analysis_results-2024-01-10.xlsx")) %>% tibble()
curated_ids_tbl <- read_excel(paste0(tbl_path,"/", "curated/curated_lagoa_ingleses-ASVs_x_amostras-2024-01-09.xlsx")) %>% tibble()
blast_tax <- read.csv("/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/tax_blast.csv", sep = ",", check.names = FALSE)

## Adicionando o RRA correto a pre_raw_results ----
{
  pre_raw_results_tbl <- pre_raw_results_tbl %>%
    mutate("Curated Relative abundance to all samples" = 0,
           "Curated Relative abundance on sample" = 0,
           "Curated Sample total abundance" = 0)
  
  abd_total <- sum(pre_raw_results_tbl$Abundance)
  
  pre_raw_results_tbl <- pre_raw_results_tbl %>%
    group_by(Unique_File_name,`Read origin`) %>%   
    mutate("Curated Sample total abundance" = sum(Abundance),
           "Curated Relative abundance to all samples" = round((Abundance/abd_total),digits = 10),
           "Curated Relative abundance on sample" =  round((Abundance/`Curated Sample total abundance`),digits = 10)) %>%
    relocate(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`,
             `Curated Sample total abundance`,`Curated Relative abundance to all samples`,`Curated Relative abundance on sample`) %>%
    ungroup() %>% 
    select(-c(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`)) %>% 
    rename("Sample total abundance" = "Curated Sample total abundance",
           "Relative abundance to all samples" = "Curated Relative abundance to all samples",
           "Relative abundance on sample"= "Curated Relative abundance on sample")
}


## Adicionando a pre_raw_results_tbl os nomes que foram curados manualmente ----

  ## 1o selecionando apenas as colunas que importam
  
  # Curated IDs
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
        "Curated Order (BLASTn)"
        # ,
        # "Obs. Curadoria",
        # "Possible contamination",
        # "Curated Order"
        )) %>%
      unique()

    View(curated_ids_tbl)

    # Verificar se existem IDs diferentes para a mesma ASV
    curated_ids_tbl$`ASV header`[which(curated_ids_tbl$`ASV header` %>% duplicated())]
  
  # Complete analysis
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
               "Obs. Curadoria",
               "Possible contamination",
               "Read origin",
               "Primer expected length",
               "ASV Size (pb)",
               "ASV header")) %>%
      mutate(Sample = str_replace(Sample, "__", "_"))
    
    View(pre_raw_results_tbl)

  ## 2o renomeando segundo os nomes curados e reordenando

  pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()

  raw_results_tbl <- pre_raw_results_tbl %>%
    left_join(curated_ids_tbl,
              by = c("ASV header")) %>%
    select(c("ASV header", #oriundas da curated_ids_tbl
             "ASV (Sequence)",
             "Curated ID",
             "Final ID (BLASTn)",
             "BLASTn pseudo-score",
             "Class (BLASTn)",
             "Curated Order (BLASTn)",
             # "Curated Order",
             "Sample", #oriundas da pre_raw_samples_tbl
             "Unique_File_name",
             "OTU",
             "Abundance",
             "Sample total abundance",
             "Relative abundance on sample",
             "Relative abundance to all samples",
             "Obs. Curadoria",
             # "Possible contamination",
             "Read origin",
             "Primer expected length",
             "ASV Size (pb)",
             "ASV header"
             ))

  View(raw_results_tbl)

  # Verificar se existem IDs diferentes para a mesma ASV
  raw_results_tbl$`ASV header`[which(curated_ids_tbl$`ASV header` %>% duplicated())]

  # Verificando o total de reads por amostra
  raw_reads <- raw_results_tbl %>% group_by(Sample) %>%
    summarise(total_abd = sum(Abundance))

  View(raw_reads)

  # Comparar com o total que saiu do pipeline para ver se nao ha perdas:

    # > raw_reads < smp_abd_ID_Final %>% group_by(Sample) %>%
    #   +  summarise(total_abd = sum(Abundance))
    # > print(raw_reads, n= 21)

    # A tibble: 21 × 2
    #   Sample             total_abd
    #   <chr>                  <dbl>
    # 1 L1_jan22              188658
    # 2 L1_nov21               66426
    # 3 L1_nov_dec_20_mi       41539
    # 4 L1_out21              384154
    # 5 L2_dez20              357293
    # 6 L2_jan22              249320
    # 7 L2_nov20              207421
    # 8 L2_nov21              184214
    # 9 L2_out21              530824
    # 10 L3_jan22                 37
    # 11 L3_nov21              49612
    # 12 L3_out21             370515
    # 13 L4_jan22                 32
    # 14 L4_nov21             514706
    # 15 L4_out21                123
    # 16 STX__L1_nov21        303077
    # 17 STX__L2_nov21        312437
    # 18 STX__L3_nov21        278561
    # 19 STX__L4_nov21          1264
    # 20 br_jan_22                44
    # 21 br_nov21                 37

    # Verificando se ha diferencas entre a tabela de ASVs e a tabela de IDs curadas
    dif_raw_pre <- pre_raw_results_tbl %>%
      anti_join(raw_results_tbl)

    View(dif_raw_pre) #vazio e' bom!

    # Apos ver que tudo que esta na tabela eh o que saiu do pipeline, podemos retirar 
    # tudo o que estiver fora do intervalo do amplicon, identificacoes NA, brancos e
    # possiveis contaminacoes
    
    filt_results_tbl <-
      raw_results_tbl %>% 
      filter(!`Primer expected length` %in% "out of range") %>% 
      filter(!`Curated ID` == "NA") %>% 
      filter(!(`Class (BLASTn)` %in% c("Mammalia", "Aves") & `BLASTn pseudo-score` < 98)) %>% 
      filter(!(Sample %in% c("EM113_NEGPCR1",
                             "EM135c4c5_NEGPCR1",
                             "EM135c4c5_NEGPCR2",
                             "EM149_NEGPCR1",
                             "EM149_NEGPCR2",
                             "EM156_NegPCR1",
                             "S724_NEGPCR1",
                             "br_jan_22",
                             "br_nov21"))) %>% 
      # filter(!(`Possible contamination` %in% c("Possible contamination"))) %>% 
      filter(!(`Obs. Curadoria` %in% c("Possible contamination")))
    
    View(filt_results_tbl)
    
    # ASVs que foram retiradas
    out_results_tbl <-
      setdiff(raw_results_tbl,filt_results_tbl)
    
    # verificando o total de reads por amostra apos filtragem
    
    filt_raw_reads <- filt_results_tbl %>% group_by(Sample) %>% 
      summarise(total_abd_filt = sum(Abundance))
    
    View(filt_raw_reads)
    
    stats_filt_reads <- raw_reads %>% 
      left_join(filt_raw_reads, 
                by = c("Sample")) %>% 
      mutate("Proportion left" = total_abd_filt / total_abd * 100)
    
    View(stats_filt_reads)
    
    # corrigindo o RRA
    
    # Verificando o RRA
    filt_results_tbl %>% 
    # grouped_by_ID_BLASTid %>%
    group_by(Sample) %>% 
      summarize(total_RRA = sum(`Relative abundance on sample`)) # Veja como que algumas amostras perderam muito!
    
    # >     filt_results_tbl %>% 
    #   +     # grouped_by_ID_BLASTid %>%
    #   +     group_by(Sample) %>% 
    #   +       summarize(total_RRA = sum(`Relative abundance on sample`))
    # # A tibble: 19 × 2
    # Sample           total_RRA
    # <chr>                <dbl>
    #   1 L1_jan22             0.307
    # 2 L1_nov21             1.00 
    # 3 L1_nov_dec_20_mi     0.997
    # 4 L1_out21             0.999
    # 5 L2_dez20             1.00 
    # 6 L2_jan22             0.489
    # 7 L2_nov20             1.00 
    # 8 L2_nov21             1.00 
    # 9 L2_out21             1.00 
    # 10 L3_jan22             0.622
    # 11 L3_nov21             1.00 
    # 12 L3_out21             1.00 
    # 13 L4_jan22             0.656
    # 14 L4_nov21             1.00 
    # 15 L4_out21             1.00 
    # 16 STX_L1_nov21         0.949
    # 17 STX_L2_nov21         0.994
    # 18 STX_L3_nov21         0.910
    # 19 STX_L4_nov21         0.739
    
    filt_results_tbl <- filt_results_tbl %>%
      mutate("Curated Relative abundance to all samples" = 0,
             "Curated Relative abundance on sample" = 0,
             "Curated Sample total abundance" = 0)
    
    abd_total <- sum(filt_results_tbl$Abundance)
    
    filt_results_tbl <- filt_results_tbl %>%
      group_by(Unique_File_name,`Read origin`) %>%   
      mutate("Curated Sample total abundance" = sum(Abundance),
             "Curated Relative abundance to all samples" = round((Abundance/abd_total),digits = 10),
             "Curated Relative abundance on sample" =  round((Abundance/`Curated Sample total abundance`),digits = 10)) %>%
      relocate(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`,
               `Curated Sample total abundance`,`Curated Relative abundance to all samples`,`Curated Relative abundance on sample`) %>%
      ungroup() %>% 
      select(-c(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`)) %>% 
      rename("Sample total abundance" = "Curated Sample total abundance",
             "Relative abundance to all samples" = "Curated Relative abundance to all samples",
             "Relative abundance on sample"= "Curated Relative abundance on sample")

    # Verificando se funcionou
    filt_results_tbl %>% 
      # grouped_by_ID_BLASTid %>%
      group_by(Sample) %>% 
      summarize(total_RRA = sum(`Relative abundance on sample`)) # Compare com o resultado anterior!
    
    ## 3o adicionando metadados que faltaram 

    filt_results_tbl$Sample %>% unique() %>% sort() %>% paste0(collapse = '",\n"') %>% cat()
    {
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
                      "Novembro 2021",
                      "Novembro 2021",
                      "Novembro 2021",
                      "Novembro 2021",
                      "Novembro 2021",
                      "Janeiro 2022",
                      "Janeiro 2022",
                      "Janeiro 2022",
                      "Janeiro 2022"
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
                "2021",
                "2021",
                "2021",
                "2021",
                "2021",
                "2022",
                "2022",
                "2022",
                "2022"
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
                    "Fundação",
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
                  "MCE",
                  "Sterivex",
                  "Sterivex",
                  "Sterivex",
                  "Sterivex",
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
               "run5_dez21",
               "EM156",
               "EM156",
               "EM156",
               "EM156",
               "EM156",
               "EM156",
               "EM156",
               "EM156"
      )
    }
    
    # As tibble
    metadata_tbl <- tibble(run, 
                           Sample,
                           filter,
                           new_name, 
                           expedition, 
                           year)

    # Mergindo os metadados com raw_results_tbl

    filt_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
    
    filt_results_tbl$Sample %>% unique()
    
    curated_full_tbl <- left_join(filt_results_tbl, metadata_tbl,
                                     by = "Sample") %>%
      select(c("Sample",
               "run",
               "new_name",
               "expedition",
               "year",
               "filter",
               "Curated ID",
               # "Curated Order",
               # "Obs. Curadoria",
               "Final ID (BLASTn)",
               "BLASTn pseudo-score",
               "Class (BLASTn)",
               "Curated Order (BLASTn)",
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
    
    View(curated_full_tbl)
    
    diff_cur_raw <- dplyr::anti_join(#curated_ids_tbl,
      filt_results_tbl,
      curated_full_tbl,
      # ver se ocorreu sem problemas o left_join
      by = "ASV header")
    
    View(diff_cur_raw) #vazio e' bom!

    # Mergindo os metadados com out_results_tbl
    
    out_results_full <- left_join(out_results_tbl, metadata_tbl,
                                  by = "Sample") %>%
      select(c("Sample",
               "run",
               "new_name",
               "expedition",
               "year",
               "filter",
               "Curated ID",
               "Curated Order (BLASTn)",
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
    
    View(out_results_full)
    
## Dados de hits por DB ----
    # Tive que rodar novamente o BLASTr por que as infos sobre a origem dos hits nao foram salvos adequadamente.
    # As tabelas do pipeline geradas em 23-02 possuem esses dados completos. 
    
    ## Obter a nova tabela raw_results
    {
      NW_pre_raw_results_tbl <- read_excel("/home/gabriel/projetos/peixes-eDNA/analises/dez_23/runs_2_4_5_EM156/results/lagoa_ingleses-Complete_analysis_results-2024-02-23.xlsx") %>% tibble()
      
      ## Adicionando o RRA correto a NW_pre_raw_results 
      {
        NW_pre_raw_results_tbl <- NW_pre_raw_results_tbl %>%
          mutate("Curated Relative abundance to all samples" = 0,
                 "Curated Relative abundance on sample" = 0,
                 "Curated Sample total abundance" = 0)
        
        abd_total <- sum(NW_pre_raw_results_tbl$Abundance)
        
        NW_pre_raw_results_tbl <- NW_pre_raw_results_tbl %>%
          group_by(Unique_File_name,`Read origin`) %>%   
          mutate("Curated Sample total abundance" = sum(Abundance),
                 "Curated Relative abundance to all samples" = round((Abundance/abd_total),digits = 10),
                 "Curated Relative abundance on sample" =  round((Abundance/`Curated Sample total abundance`),digits = 10)) %>%
          relocate(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`,
                   `Curated Sample total abundance`,`Curated Relative abundance to all samples`,`Curated Relative abundance on sample`) %>%
          ungroup() %>% 
          select(-c(`Sample total abundance`,`Relative abundance to all samples`,`Relative abundance on sample`)) %>% 
          rename("Sample total abundance" = "Curated Sample total abundance",
                 "Relative abundance to all samples" = "Curated Relative abundance to all samples",
                 "Relative abundance on sample"= "Curated Relative abundance on sample")
      }
      
      ## Selecionando colunas
      NW_pre_raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
      
      NW_pre_raw_results_tbl <-
        NW_pre_raw_results_tbl %>%
        select(c("Sample",
                 "Unique_File_name",
                 "OTU",
                 "Abundance",
                 "Sample total abundance",
                 "Relative abundance on sample",
                 "Relative abundance to all samples",
                 "Obs. Curadoria",
                 "Possible contamination",
                 "Read origin",
                 "Primer expected length",
                 "ASV Size (pb)",
                 "ASV header")) %>%
        mutate(Sample = str_replace(Sample, "__", "_"))
      
      View(NW_pre_raw_results_tbl)
      
      ## Adicionando a info de contaminacao
      
      NW_contam <- pre_raw_results_tbl %>% 
        select(`Obs. Curadoria`,
               `Possible contamination`,
               `ASV header`)
      
      View(NW_contam)
      
      NW_pre_raw_results_tbl <- NW_pre_raw_results_tbl %>% 
        left_join(NW_contam,
                  by = "ASV header")
      
      View(NW_pre_raw_results_tbl)
      
      # Verificar se os resultados sao os mesmos com pre_raw_results_tbl

      dif_NW_pre <- NW_pre_raw_results_tbl %>%
        anti_join(pre_raw_results_tbl)
      
      View(dif_NW_pre)
    }
    
    ## Obter a nova tabela curated_IDs
    NW_curated_ids_tbl <- read_excel("/home/gabriel/projetos/peixes-eDNA/analises/dez_23/runs_2_4_5_EM156/results/lagoa_ingleses-ASVs_x_amostras-2024-02-23.xlsx") %>% tibble()
    
    ## Selecionando apenas as colunas que importam para comparar
    
    # Curated IDs
  
    TT_curated_ids_tbl <-
      curated_ids_tbl %>%
      select(c(
        "ASV header",
        "ASV (Sequence)",
        # "Curated ID",
        "Final ID (BLASTn)",
        "BLASTn pseudo-score",
        "Class (BLASTn)",
        # "Curated Order (BLASTn)"
        # ,
        # "Obs. Curadoria",
        # "Possible contamination",
        # "Curated Order"
      )) %>%
      unique()
    
    View(TT_curated_ids_tbl)
    
    NWt_curated_ids_tbl <-
    NWt_curated_ids_tbl %>%
      select(c(
        "ASV header",
        "ASV (Sequence)",
        # "Curated ID",
        "Final ID (BLASTn)",
        "BLASTn pseudo-score",
        "Class (BLASTn)",
        # "Curated Order (BLASTn)"
        # ,
        # "Obs. Curadoria",
        # "Possible contamination",
        # "Curated Order"
      )) %>%
      unique()
    
    View(NWt_curated_ids_tbl)
    
    # Verificar se os resultados sao os mesmos com pre_raw_results_tbl
    
    dif_NW_pre <- NW_curated_ids_tbl %>%
      anti_join(curated_ids_tbl)    
    
    View(dif_NW_pre)
    
    # Podemos ver como a diferenca entre as tabelas e` unicamente a correcao do Final ID, mas aparentemente
    # os resultados do BLASTn foram identicos. Posso entao usar a ASV_header para pegar os resultados do DB
    
    hits_DB <- NW_curated_ids_tbl %>% 
      select(`ASV header`,
             `Final ID (BLASTn)`,
             `1_DB`,
             `2_DB`,
             `3_DB`)
    
    View(hits_DB)
    
    ## Tabelas ----
    
    # Criacao da lista com os possiveis nomes atribuidos as ASVs
    
    curated_full_tbl %>% colnames()
    curated_full_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
    
    # Agrupamento da ASVs que possuem os mesmos atributos abaixo
    
    grouped_by_ID_BLASTid <- curated_full_tbl %>%
      relocate(c("Sample",
                 "expedition",
                 "year",
                 "Final ID (BLASTn)",
                 "BLASTn pseudo-score",
                 "Curated ID",
                 "OTU",
                 "filter",
                 "Curated Order (BLASTn)",
                 "Relative abundance on sample",
                 "new_name",
                 "Class (BLASTn)",
                 # "Order (DADA2)",
                 # "Curated Order",
                 # "Curated Max. taxonomy",
                 "Abundance",
                 "ASV header",
                 "ASV (Sequence)",
                 "Primer expected length"
      )) %>%
      group_by(Sample, `Curated ID`,`Final ID (BLASTn)`,`BLASTn pseudo-score`, expedition, new_name, year, filter, `ASV header`, `Curated Order (BLASTn)`) %>%
      summarize("RRA" = sum(`Relative abundance on sample`),
                "Class (BLASTn)" = unique(`Class (BLASTn)`),
                "OTU" = unique(OTU),
                # "Curated Order" = unique(`Curated Order`),
                # "Order (DADA2)" = unique(`Order (DADA2)`),
                # "Curated Max. taxonomy" = unique(`Curated Max. taxonomy`),
                "Abundance" = sum(Abundance),  # Agregando a coluna Abundance
                "ASV (Sequence)" = unique(`ASV (Sequence)`),
                "Primer expected length" = unique(`Primer expected length`)
      ) %>%
      mutate(Nivel = case_when(`year` == 2020 ~ "Cheio", 
                               `year` == 2022 ~ "Cheio", 
                               TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o nivel da lagoa
      ungroup()
    
    View(grouped_by_ID_BLASTid)
    
    grouped_by_ID_BLASTid[which(grouped_by_ID_BLASTid %>% duplicated())]
      
    # verificando o total de reads por amostra
    
    grouped_reads <- grouped_by_ID_BLASTid %>% group_by(Sample) %>% 
      summarise(total_abd = sum(Abundance))
    
    View(grouped_reads)
    
    # comparar com o total que saiu do pipeline para ver se nao ha perdas:
    
    # > raw_reads < smp_abd_ID_Final %>% group_by(Sample) %>% 
    #   +  summarise(total_abd = sum(Abundance))
    # > print(raw_reads, n= 21)
    
    # A tibble: 21 × 2
    #   Sample             total_abd
    #   <chr>                  <dbl>
    # 1 L1_jan22              188658
    # 2 L1_nov21               66426
    # 3 L1_nov_dec_20_mi       41539
    # 4 L1_out21              384154
    # 5 L2_dez20              357293
    # 6 L2_jan22              249320
    # 7 L2_nov20              207421
    # 8 L2_nov21              184214
    # 9 L2_out21              530824
    # 10 L3_jan22                 37
    # 11 L3_nov21              49612
    # 12 L3_out21             370515
    # 13 L4_jan22                 32
    # 14 L4_nov21             514706
    # 15 L4_out21                123
    # 16 STX__L1_nov21        303077
    # 17 STX__L2_nov21        312437
    # 18 STX__L3_nov21        278561
    # 19 STX__L4_nov21          1264
    # 20 br_jan_22                44
    # 21 br_nov21                 37
    
    # adicionando a taxonomia do BLAST
    
    blast_tax <- blast_tax[2:86] %>% 
      rename(`ASV (Sequence)` = "ASV")
    
    blast_tax_less <- blast_tax %>%
      select(c("ASV (Sequence)",
               "Order (BLASTn)",
               "Family (BLASTn)",
               "Genus (BLASTn)"
      )) 
    grouped_by_ID_BLASTid %>% colnames()
    
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
               "Curated Order (BLASTn)",
               # "Order (BLASTn)",
               # "Curated Order",
               "Family (BLASTn)",
               "Genus (BLASTn)",
               "Curated genus",
               "expedition",
               "new_name",
               "year", 
               "filter",
               "ASV header",
               "RRA",
               "OTU", 
               "Abundance",
               "ASV (Sequence)",
               "Primer expected length",
               "Nivel" 
      )) %>% 
      unique()
    
    View(grouped_by_ID_BLASTid)
    
    # Agrupamento da ASVs que possuem os mesmos atributos abaixo VERSAO OUT ----
    {
      out_grouped_by_ID_BLASTid <- out_results_full %>%
        relocate(c("Sample",
                   "expedition",
                   "year",
                   "Final ID (BLASTn)",
                   "BLASTn pseudo-score",
                   "Curated ID",
                   "OTU",
                   "filter",
                   # "Curated Order",
                   "Relative abundance on sample",
                   "new_name",
                   "Class (BLASTn)",
                   # "Order (DADA2)",
                   # "Curated Order",
                   # "Curated Max. taxonomy",
                   "Abundance",
                   "ASV header",
                   "ASV (Sequence)",
                   "Primer expected length"
        )) %>%
        group_by(Sample, `Curated ID`,`Final ID (BLASTn)`,`BLASTn pseudo-score`, expedition, new_name, year, filter, `ASV header`) %>%
        summarize("RRA" = sum(`Relative abundance on sample`),
                  "Class (BLASTn)" = unique(`Class (BLASTn)`),
                  "OTU" = unique(OTU),
                  # "Curated Order" = unique(`Curated Order`),
                  # "Order (DADA2)" = unique(`Order (DADA2)`),
                  # "Curated Max. taxonomy" = unique(`Curated Max. taxonomy`),
                  "Abundance" = sum(Abundance),  # Agregando a coluna Abundance
                  "ASV (Sequence)" = unique(`ASV (Sequence)`),
                  "Primer expected length" = unique(`Primer expected length`)
        ) %>%
        mutate(Nivel = case_when(`year` == 2020 ~ "Cheio", 
                                 `year` == 2022 ~ "Cheio", 
                                 TRUE ~ "Vazio")) %>%  ## com essa linha a gente inclui o nivel da lagoa
        ungroup()
      
      View(out_grouped_by_ID_BLASTid)
      
      # adicionando a taxonomia do BLAST
      
      out_grouped_by_ID_BLASTid <- out_grouped_by_ID_BLASTid %>% 
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
                 # "Curated Order",
                 "Family (BLASTn)",
                 "Genus (BLASTn)",
                 "Curated genus",
                 "expedition",
                 "new_name",
                 "year", 
                 "filter",
                 "ASV header",
                 "RRA",
                 "OTU", 
                 "Abundance",
                 "ASV (Sequence)",
                 "Primer expected length",
                 "Nivel" 
        )) %>% 
        unique()
      
      View(out_grouped_by_ID_BLASTid)
      
      write.csv(out_grouped_by_ID_BLASTid, "~/projetos/lagoa_ingleses/results/tabelas/out_grouped_2024.csv")
    }
    

    
