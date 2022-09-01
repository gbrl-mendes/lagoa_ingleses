---
  title: "Análises de composição de espécies por amostras coletadas em época chuvosa e época seca na Lagoa dos Ingleses"
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
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v2.csv"), sep = ",", check.names = FALSE) %>% tibble()
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Oreochromis niloticus")] <- "Tilapia rendalli" # Oreochromis niloticus é Tilapia
}

# Tabela wider para NMDS, Beta Diversidade e PCA ----

# Essa tabela e uma matriz onde as linhas sao as amostras e as colunas sao as especies
# Cada valor representa a abundancia da especie na amostra. Essa matriz pode ser usada
# em diversas analises

# Para realizar as analises que envolvam apenas as abundancias das especies,
# usar as colunas de 4 a 35 c[,4:35]

## Tabela só peixes
{
  spp_sample_tbl_f <- raw_results_tbl %>%
    select(c(
      "Sample",
      # "Run",
      # "Group",
      "Expedition",
      # "Coleta",
      "Year",
      #"Sample.Name",
      #"File_name",
      # "OTU",
      # "Read_origin",
      "Curated ID",
      # "final ID",
      # "Abundance",
      # "Relative abundance to all samples",
      "Relative abundance on sample",
      # "Sample.total.abundance",
      # "Type",
        "Point")) %>%
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
      filter(!Sample %in% c("LI1-neo-mi",
                            "L2_dez20",
                            "L2_nov20")
      ) %>% 
      group_by(Sample,
               Expedition,
               Point,
               Year,
               `Curated ID`
      ) %>%
      summarize(`ID Abundance on sample (%)` = sum(`Relative abundance on sample`)
      ) %>%
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
      unique()  %>% 
      pivot_wider(names_from = `Curated ID`, 
                  values_from = `ID Abundance on sample (%)`)
  
  ## incluindo a coluna de forma para fazer o plot
  spp_sample_tbl_f <- spp_sample_tbl_f %>%  mutate("Shape" = if_else(Expedition %in% c("out/21"), 1, 2)) %>%
    relocate("Shape") # passar a coluna para o começo
  
  ## transformando a tabela em um data frame
  spp_sample_tbl_f <- as.data.frame(spp_sample_tbl_f) 
  
  ## nomes das amostras como row names
  row.names(spp_sample_tbl_f) <- spp_sample_tbl_f$Sample 
  
  ## retirando a coluna Samples
  spp_sample_tbl_f <- spp_sample_tbl_f %>% select(-c(Sample)) 
  
  ## substituindo NA por 0
  spp_sample_tbl_f[is.na(spp_sample_tbl_f)] = 0 
  
  ## Substituir o espaco no nome das especies por _
  colnames(spp_sample_tbl_f) <- colnames(spp_sample_tbl_f) %>%  
    str_replace_all(pattern = " ",replacement = "_")
}

## Tabela todas Ids a nível de especie
{
  spp_sample_tbl_all <- raw_results_tbl %>%
    select(c(
      "Sample",
      # "Run",
      # "Group",
      "Expedition",
      # "Coleta",
      "Year",
      #"Sample.Name",
      #"File_name",
      # "OTU",
      # "Read_origin",
      "Curated ID",
      # "final ID",
      # "Abundance",
      # "Relative abundance to all samples",
      "Relative abundance on sample",
      # "Sample.total.abundance",
      # "Type",
      "Point")) %>%
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
                                "Cutibacterium acnes"
    )) %>%
    filter(!Sample %in% c("LI1-neo-mi",
                          "L2_dez20",
                          "L2_nov20")
    ) %>% 
    group_by(Sample,
             Expedition,
             Point,
             Year,
             `Curated ID`
    ) %>%
    summarize(`ID Abundance on sample (%)` = sum(`Relative abundance on sample`)
    ) %>%
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
    unique()  %>% 
    pivot_wider(names_from = `Curated ID`, 
                values_from = `ID Abundance on sample (%)`)
  
  ## incluindo a coluna de forma para fazer o plot
  spp_sample_tbl_all <- spp_sample_tbl_all %>%  mutate("Shape" = if_else(Expedition %in% c("out/21"), 1, 2)) %>%
    relocate("Shape") # passar a coluna para o começo
  
  ## transformando a tabela em um data frame
  spp_sample_tbl_all <- as.data.frame(spp_sample_tbl_all) 
  
  ## nomes das amostras como row names
  row.names(spp_sample_tbl_all) <- spp_sample_tbl_all$Sample 
  
  ## retirando a coluna Samples
  spp_sample_tbl_all <- spp_sample_tbl_all %>% select(-c(Sample)) 
  
  ## substituindo NA por 0
  spp_sample_tbl_all[is.na(spp_sample_tbl_all)] = 0
  
  ## Substituir o espaco no nome das especies por _
  colnames(spp_sample_tbl_all) <- colnames(spp_sample_tbl_all) %>%  
    str_replace_all(pattern = " ",replacement = "_")
}

## Tabela todas amostras todas ids e samples
{
  spp_sample_tbl_tudo <- raw_results_tbl %>%
    select(c(
      "Sample",
      # "Run",
      # "Group",
      "Expedition",
      # "Coleta",
      "Year",
      #"Sample.Name",
      #"File_name",
      # "OTU",
      # "Read_origin",
      "Curated ID",
      # "final ID",
      # "Abundance",
      # "Relative abundance to all samples",
      "Relative abundance on sample",
      # "Sample.total.abundance",
      # "Type",
      "Point")) %>%
    filter(!`Curated ID` %in% c("Actinopteri",
                                "Astyanax",
                                "Characidae",
                                "Characidium",
                                "Characiformes",
                                "Cichla",
                                "Cichlidae",
                                "Hoplias",
                                "Pimelodus",
                                "",
                                "Cutibacterium acnes"
                                )) %>% 
    group_by(Sample,
             Expedition,
             Point,
             Year,
             `Curated ID`
    ) %>%
    summarize(`ID Abundance on sample (%)` = sum(`Relative abundance on sample`)
    ) %>%
    mutate(Sample = factor(Sample, levels = c("L1_out21",
                                              "L1_nov21",
                                              "L2_out21",
                                              "L2_nov21",
                                              "L3_out21",
                                              "L3_nov21",
                                              "L4_out21",
                                              "L4_nov21",
                                              "L2_dez20",
                                              "L2_nov20",
                                              "LI1-neo-mi"
    ))) %>%
    ungroup() %>%
    unique()  %>% 
    pivot_wider(names_from = `Curated ID`, 
                values_from = `ID Abundance on sample (%)`)
  
  ## incluindo a coluna de forma para fazer o plot
  spp_sample_tbl_tudo <- spp_sample_tbl_tudo %>%  
    mutate("Shape" = case_when(
      Expedition %in% c("Nov_Dec/20") ~"0", 
      Expedition %in% c("Dec/20") ~"1",
      Expedition %in% c("Nov/20") ~"2",
      Expedition %in% c("Nov/21") ~"5",
      Expedition %in% c("out/21") ~"6")) %>% 
    relocate("Shape") # passar a coluna para o começo
  
  ## transformando a tabela em um data frame
  spp_sample_tbl_tudo <- as.data.frame(spp_sample_tbl_tudo) 
  
  ## nomes das amostras como row names
  row.names(spp_sample_tbl_tudo) <- spp_sample_tbl_tudo$Sample 
  
  ## retirando a coluna Samples
  spp_sample_tbl_tudo <- spp_sample_tbl_tudo %>% select(-c(Sample)) 
  
  ## substituindo NA por 0
  spp_sample_tbl_tudo[is.na(spp_sample_tbl_tudo)] = 0 
  
  ## Substituir o espaco no nome das especies por _
  colnames(spp_sample_tbl_tudo) <- colnames(spp_sample_tbl_tudo) %>%  
    str_replace_all(pattern = " ",replacement = "_")
}

# Bar Plots Proporcao Especies/Ids ----

## Barplots com os resultados das proporcoes

## Utilizando a abundancia relativa nas amostras para calcular a proporcao de peixes
{
  grouped_by_ID_especie_tbl <- raw_results_tbl %>%
    select(c(
      "Sample",
      "Expedition",
      "Year",
      "Sample.Name",
      "File_name",
      "Curated ID",
      "Relative abundance to all samples",
      "Relative abundance on sample",
      "Point",
      "Abundance",
      "Sample total abundance"
    )) %>% 
    mutate("Peixe" = if_else (`Curated ID` %in% c("",
                                                  "Cutibacterium acnes",
                                                  "Bos taurus",
                                                  "Canis familiaris",
                                                  "Didelphis albiventris (Gamba)",
                                                  "Homo sapiens",
                                                  "Hydrochaeris hydrochaeris (Capivara)",
                                                  "Nannopterum brasilianus",
                                                  "Oryctolagus cuniculus (Coelho-bravo)",
                                                  "Progne chalybea (Andorinha-grande)",
                                                  "Sus scrofa"),"FALSE", "TRUE")) %>%
    mutate("Especie" = if_else (`Curated ID` %in% c("Hoplias",
                                                    "Cichla",
                                                    "Actinopteri",
                                                    "Characiformes",
                                                    "Siluriformes",
                                                    "Astyanax",
                                                    "Characidium",
                                                    "Cichlidae",
                                                    "Characidae",
                                                    "Pimelodus"), "FALSE", "TRUE")) %>% 
    group_by(Peixe, Especie) %>% 
    summarize("Peixe" = unique(Peixe), #identificacao se as ASVs sao de peixes ou nao-peixes
              "Especie" = unique(Especie), #identificacao se as ASVs foram identificadas ao nivel de especie
              # `RRA` = sum(`Relative abundance on sample`)/11, #proporcao calculada com a RRA
              `ASVs` = sum(`Abundance`), #numero de ASVs
              `Proporcao` = sum(`ASVs` / 2730833 * 100) #proporcao calculada com o total de ASVs
    ) %>%
    ungroup()
}

## Proporcoes barplots
{
  ### Proporcao de peixes
  {peixe <- c("Não-peixes", "Peixes")
  ASVs <- c(25130, (160274 + 2545429))
  proporcao <- c(0.92, (5.87 + 93.21))
  peixes <- data_frame(peixe, ASVs, proporcao)
  }
  ## Barplots da proporcao de peixes
  {
    print(plot_peixe <- peixes %>% 
            ggplot(aes(x=peixe, y=proporcao, fill = peixe)) + 
            geom_bar(stat = "identity") +
            geom_text(aes(label = proporcao), ## colocar o valor da proporcao acima das barras
                      position = position_dodge(width = 0.9),
                      vjust = -0.1,
                      colour = "black", size = 6) +
            guides(col = guide_legend(nrow = 6)) +
            xlab("Espécies") + ## alterar o nome do eixo x
            ylab("Proporção de ASVs %") + ## alterar o nome do eixo y
            ggtitle(label = "Proporção de identificações associadas a peixes") + ## alterar o titulo do plot
            theme_bw(base_size = 16) +
            scale_fill_manual(values = viridis::viridis(n=4)[c(2,3)]) +
            guides(fill = guide_legend("")) + ## alterar o titulo da legenda
            theme(plot.title = element_text(hjust=0.3))
    )
    ## Plotando
    ggsave(plot = plot_peixe, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/barplots/plot_peixes.pdf",
           device = "pdf", units = "cm", height = 15, width = 15, dpi = 600)
  }
  
  ### Proporcao de ids
  especie <- c("Espécie", "Outros")
  e_ASVs <- c(2545429, 160274)
  e_proporcao <- c(0.94, 0.06)
  especies <- data_frame(especie, e_ASVs, e_proporcao)
  
  ## Barplots da proporcao de ids
  
  print(plot_ids <- especies %>% 
          ggplot(aes(x = especie, y = e_proporcao, fill = especie)) + 
          geom_bar(stat = "identity") +
          geom_text(aes(label = e_proporcao), ## colocar o valor da proporcao acima das barras
                    position = position_dodge(width = 0.9),
                    vjust = -0.1,
                    colour = "black", size = 6) +
          guides(col = guide_legend(nrow = 6)) +
          xlab("Nível taxonômico") + ## alterar o nome do eixo x
          ylab("Proporção de ASVs %") + ## alterar o nome do eixo y
          ggtitle(label = "Proporção de identificações e \n respectivos níveis taxonômicos") + ## alterar o titulo do plot
          theme_bw(base_size = 16) +
          scale_fill_manual(values = viridis::viridis(n=4)[c(2,3)]) +
          guides(fill = guide_legend("Nível")) + ## alterar o titulo da legenda
          theme(plot.title = element_text(hjust=0.5))
  )
  ## Plotando
  ggsave(plot = plot_ids, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/barplots/plot_ids.pdf",
         device = "pdf", units = "cm", height = 15, width = 15, dpi = 600)
}



# Tile Plots ----

## Criacao das tabelas all_ids_tbl e few_ids_tbl
{
  ## Criacao da lista com os possiveis nomes atribuidos as ASVs
  {
    raw_results_tbl %>% colnames()
    raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
  }
  
  ## Agrupamento da ASVs que possuem os mesmos atributos abaixo
  {
    grouped_by_ID_tbl <- raw_results_tbl %>%
      select(c(
        "Sample",
        "Expedition",
        "Year",
        "Sample.Name",
        "File_name",
        "Curated ID",
        "Relative abundance on sample",
        "Point"
      )) %>% 
      group_by(Sample, `Curated ID`, Expedition, Point, Sample.Name, File_name) %>% 
      summarize(`Sample` = unique(Sample),
                `Curated ID` = unique(`Curated ID`),
                `Expedition` = unique(Expedition),
                `Year` = unique(Year),
                `Point` = unique(Point),
                `Sample.Name` = unique(Sample.Name),
                `File_name` = unique(File_name),
                `RRA` = sum(`Relative abundance on sample`)) %>%
      ungroup()
  }
  
  ## Organizar as especies
  {
    grouped_by_ID_tbl$`Curated ID` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
    grouped_by_ID_tbl$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  }
  
  ## Organizar as ordem das especies usando fatores
  
  ### com todas as identificacoes
  {
    all_ID_tbl <- grouped_by_ID_tbl %>%
      mutate(`Curated ID` = factor(`Curated ID`,
                                   levels = rev(c(
                                     "",
                                     "Actinopteri",
                                     "Acinocheirodon melanogramma",
                                     "Astyanax",
                                     "Astyanax fasciatus",
                                     "Astyanax lacustris",
                                     "Brycon orthotaenia",
                                     "Bryconamericus stramineus",
                                     "Characidae",
                                     "Characidium",
                                     "Characiformes",
                                     "Cichla",
                                     "Cichlidae",
                                     "Colossoma macropomum",
                                     "Coptodon zillii",
                                     "Eigenmannia virescens",
                                     "Gymnotus carapo",
                                     "Hemigrammus gracilis",
                                     "Hemigrammus marginatus",
                                     "Hoplias",
                                     "Hoplias intermedius",
                                     "Hoplias malabaricus",
                                     "Hypomasticus steindachneri",
                                     "Leporellus vittatus",
                                     "Leporinus piau",
                                     "Leporinus reinhardti",
                                     "Leporinus taeniatus",
                                     "Megaleporinus elongatus",
                                     "Megaleporinus garmani",
                                     "Moenkhausia costae",
                                     "Myleus micans",
                                     "Oreochromis niloticus",
                                     "Orthospinus franciscensis",
                                     "Pimelodus",
                                     "Pimelodus fur",
                                     "Pimelodus maculatus",
                                     "Pimelodus pohli",
                                     "Planaltina myersi",
                                     "Poecilia reticulata",
                                     "Prochilodus costatus",
                                     "Pseudoplatystoma corruscans",
                                     "Pygocentrus piraya",
                                     "Rhamdia quelen",
                                     "Salmo salar",
                                     "Serrasalmus brandtii",
                                     "Tilapia rendalli",
                                     "Wertheimeria maculata",
                                     #nao-peixes
                                     "Bos taurus",
                                     "Canis familiaris",
                                     "Cavia magna",
                                     "Cutibacterium acnes",
                                     "Didelphis albiventris (Gamba)",
                                     "Homo sapiens",
                                     "Hydrochaeris hydrochaeris (Capivara)",
                                     "Nannopterum brasilianus",
                                     "Oryctolagus cuniculus (Coelho-bravo)",
                                     "Progne chalybea (Andorinha-grande)",
                                     "Sus scrofa"
                                   ))))
    }
  
  ### Sem os grupos
  {
    few_ID_tbl <- grouped_by_ID_tbl %>% #ids apenas ao nivel de especie e sem a bacteria
      mutate(`Curated ID` = factor(`Curated ID`,
                                   levels = rev(c(
                                     #"",
                                     #"Actinopteri",
                                     "Acinocheirodon melanogramma",
                                     #"Astyanax",
                                     "Astyanax fasciatus",
                                     "Astyanax lacustris",
                                     "Brycon orthotaenia",
                                     "Bryconamericus stramineus",
                                     #"Characidae",
                                     #"Characidium",
                                     #"Characiformes",
                                     #"Cichla",
                                     #"Cichlidae",
                                     "Colossoma macropomum",
                                     "Coptodon zillii",
                                     "Eigenmannia virescens",
                                     "Gymnotus carapo",
                                     "Hemigrammus gracilis",
                                     "Hemigrammus marginatus",
                                     #"Hoplias",
                                     "Hoplias intermedius",
                                     "Hoplias malabaricus",
                                     "Hypomasticus steindachneri",
                                     "Leporellus vittatus",
                                     "Leporinus piau",
                                     "Leporinus reinhardti",
                                     "Leporinus taeniatus",
                                     "Megaleporinus elongatus",
                                     "Megaleporinus garmani",
                                     "Moenkhausia costae",
                                     "Myleus micans",
                                     "Oreochromis niloticus",
                                     "Orthospinus franciscensis",
                                     #"Pimelodus",
                                     "Pimelodus fur",
                                     "Pimelodus maculatus",
                                     "Pimelodus pohli",
                                     "Planaltina myersi",
                                     "Poecilia reticulata",
                                     "Prochilodus costatus",
                                     "Pseudoplatystoma corruscans",
                                     "Pygocentrus piraya",
                                     "Rhamdia quelen",
                                     "Salmo salar",
                                     "Serrasalmus brandtii",
                                     "Tilapia rendalli",
                                     "Wertheimeria maculata",
                                     #nao-peixes
                                     "Cavia magna",
                                     #"Cutibacterium acnes",
                                     "Bos taurus",
                                     "Canis familiaris",
                                     "Didelphis albiventris (Gamba)",
                                     "Homo sapiens",
                                     "Hydrochaeris hydrochaeris (Capivara)",
                                     "Nannopterum brasilianus",
                                     "Oryctolagus cuniculus (Coelho-bravo)",
                                     "Progne chalybea (Andorinha-grande)",
                                     "Sus scrofa"
                                   ))))
  }
  }

## Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5

### Por ponto

  ## Tile plot filtrando as ids
  {
    tile_ponto_few <- few_ID_tbl %>%
      mutate(Expedition = factor(Expedition)) %>% # transformando variaveis categoricas em fatores com niveis
      mutate(Sample = factor(Sample,
                             levels = c("L1_nov21",
                                        "L1_out21",
                                        "L2_nov21",
                                        "L2_out21",
                                        "L3_nov21",
                                        "L3_out21",
                                        "L4_nov21",
                                        "L4_out21",
                                        "L1-neo-mi",
                                        "L2 Dez/20",
                                        "L2 Nov/20"
                                        ))) %>%
      mutate(Expedition = factor(Expedition,
                                 levels = c("Nov_Dec/20",
                                            "Nov/20",
                                            "Dec/20",
                                            "out/21",
                                            "Nov/21"
                                            ))) %>%
      mutate(Point = factor(Point)) %>%
      mutate(File_name = factor(File_name)) %>%
      mutate(Expedition = factor(Expedition)) %>%
      # filter(RRA >= 0.01) %>% # retirando as ASVs "espurias", com abundancia menor que 0.01 ATUAL: Heron pediu para voltar com elas
      filter(!is.na(`Curated ID`)) %>%   # retirar os NA
      filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                               "Nov/21")) %>%
      ### Tile plot
      ggplot(aes(y = `Curated ID`, 
                 # x = Point,
                 x = Expedition,
                 fill = RRA
                 )) +
      geom_tile() +
      # geom_text(aes(label= RRA),
      geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2))), # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
                    colour = "black", size = 3
                    ) +
        # facet_grid(~Expedition, # facetando por mes
      facet_grid(~Point, # facetando por ponto
                   scales = "free_x",
                   space = "free_x") +
      labs(fill ='Abundância \nrelativa (%)',
               x = "Amostras",
               y= "Espécies") +
      ggtitle(label = "Espécies identificadas em cada ponto") +
      geom_hline(yintercept = c(8.5)) + # linha que separa as especies de peixes e nao-peixes
      scale_fill_continuous(
        trans = "log10", # exibir a abundancia em escala logaritmica para favorecer a exibicao de baixo RRA
        breaks = c(0.001, 0.01, 0.1, 1, 10, 75), # definindo os valores que aparecem na escala
        type = "viridis") +
      theme(text=element_text(size = 10, face = "bold")) +
      guides(col = guide_legend(nrow = 6)) +
      theme_bw(base_size = 16) +
      theme(plot.title = element_text(hjust = 0.5))
      
      ## Plotando
      ggsave(plot = tile_ponto_few, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/tile_plots/tile_ponto_few.pdf",
           device = "pdf", units = "cm", height = 20, width = 30, dpi = 600)
}

  ## Title plot exibindo todas as ids
  {
    tile_ponto_all <- all_ID_tbl %>%
      mutate(Expedition = factor(Expedition)) %>% # transformando variaveis categoricas em fatores com niveis
      mutate(Sample = factor(Sample,
                             levels = c("L1_nov21",
                                        "L1_out21",
                                        "L2_nov21",
                                        "L2_out21",
                                        "L3_nov21",
                                        "L3_out21",
                                        "L4_nov21",
                                        "L4_out21",
                                        "L1-neo-mi",
                                        "L2 Dez/20",
                                        "L2 Nov/20"
                                        ))) %>%
      mutate(Expedition = factor(Expedition,
                                 levels = c("Nov_Dec/20",
                                            "Nov/20",
                                            "Dec/20",
                                            "out/21",
                                            "Nov/21"
                                            ))) %>%
      mutate(Point = factor(Point)) %>%
      mutate(File_name = factor(File_name)) %>%
      mutate(Expedition = factor(Expedition)) %>%
      # filter(RRA >= 0.01) %>% # retirando as ASVs "espurias", com abundancia menor que 0.01 #exibindo todas as Ids, mesmo com RRA abaixo de 0.01
      # filter(!is.na(`Curated ID`)) %>%   # retirar os NA
      # filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
      #                          "Nov/21")) %>%
      ### Tile plot
      ggplot(aes(y = `Curated ID`, 
                 # x = Point,
                 x = Expedition,
                 fill = RRA
                 )) +
      geom_tile() +
      geom_text(aes(label= sprintf("%0.2f", round(RRA, digits = 2)))) + # exibindo os valores de RRA dentro dos tiles para facilitar a discussao (opcional)
      # facet_grid(~Expedition, # facetando por mes
      facet_grid(~Point, # facetando por ponto
                 scales = "free_x",
                 space = "free_x") +
      labs(fill ='Abundância \nrelativa (%)',
           x = "Amostras",
           y = "Espécies") +
      ggtitle(label = "Espécies identificadas em cada ponto e mês (todas ids)") +
      geom_hline(yintercept = c(11.5)) +
      scale_fill_continuous(
        trans = "log10", # exibir a abundancia em escala logaritmica para favorecer a exibicao de baixo RRA
          breaks = c(0.001, 0.01, 0.1, 1, 10, 75), # definindo os valores que aparecem na escala
        type = "viridis") +
      theme(text=element_text(size = 14)) +
      guides(col = guide_legend(nrow = 6)) +
      theme(plot.title = element_text(hjust=0.5))
    
    ## Plotando
    ggsave(plot = tile_ponto_all, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/tile_plots/tile_ponto_all.pdf",
           device = "pdf", units = "cm", height = 25, width = 30, dpi = 600)
}

### Por Ano

  ## Exibindo por ano filtrando ids
  {
    tile_ano_few <- few_ID_tbl %>%
      mutate(Expedition = factor(Expedition)) %>% #transformando variaveis categoricas em fatores com niveis
      mutate(Sample = factor(Sample,
                           levels = c("L1_nov21",
                                      "L1_out21",
                                      "L2_nov21",
                                      "L2_out21",
                                      "L3_nov21",
                                      "L3_out21",
                                      "L4_nov21",
                                      "L4_out21",
                                      "L1-neo-mi",
                                      "L2 Dez/20",
                                      "L2 Nov/20"
                           ))) %>%
      mutate(Expedition = factor(Expedition,
                               levels = c("Nov_Dec/20",
                                          "Nov/20",
                                          "Dec/20",
                                          "out/21",
                                          "Nov/21"
                               ))) %>%
      mutate(Point = factor(Point)) %>%
      mutate(File_name = factor(File_name)) %>%
      mutate(Expedition = factor(Expedition)) %>%
      filter(!`Curated ID` %in% c("Astyanax", #retirando as ids a nivel de genero
                              "Characidae",
                              "Cichlidae",
                              "Hoplias",
                              "Pimelodus")) %>%
      # filter(RRA >= 0.01) %>% #retirando as ASVs "espurias", com abundancia menor que 0.01 ATUAL. Heron pediu para voltar com essas 
      filter(!is.na(`Curated ID`)) %>% #retirar os NA
      ### Tile plot
      ggplot(aes(y = `Curated ID`, 
               x = Year,
               fill = RRA,
               )) +
      scale_x_continuous(breaks = 0:2100) +
      geom_tile() + 
      labs(fill='Abundância \nrelativa (%)',
           x = "Ano",
           y = "Espécies") +
      geom_hline(yintercept = c(10.5)) +
      scale_fill_continuous(
        trans = "log10", # exibir a abundancia em escala logaritmica para favorecer a exibicao de baixo RRA
        breaks = c(0.001, 0.01, 0.1, 1, 10, 75), # definindo os valores que aparecem na escala
        type = "viridis") +
      ggtitle(label = "Espécies identificadas por ano") +
      theme(text = element_text(size = 14)) +
      guides(col = guide_legend(nrow = 6)) +
      theme(plot.title = element_text(hjust = 0.5))
    
    ## Plotando
    ggsave(plot = tile_ano_few, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/tile_plots/tile_ano_few.pdf",
           device = "pdf", units = "cm", height = 25, width = 20, dpi = 600)
    }

  ## Exibindo  por ano com todas ids
  {
  tile_ano_all <- all_ID_tbl %>%
      mutate(Expedition = factor(Expedition)) %>% #transformando variaveis categoricas em fatores com niveis
      mutate(Sample = factor(Sample,
                             levels = c("L1_nov21",
                                        "L1_out21",
                                        "L2_nov21",
                                        "L2_out21",
                                        "L3_nov21",
                                        "L3_out21",
                                        "L4_nov21",
                                        "L4_out21",
                                        "L1-neo-mi",
                                        "L2 Dez/20",
                                        "L2 Nov/20"
                                        ))) %>%
      mutate(Expedition = factor(Expedition,
                                 levels = c("Nov_Dec/20",
                                            "Nov/20",
                                            "Dec/20",
                                            "out/21",
                                            "Nov/21"
                                            ))) %>%
      mutate(Point = factor(Point)) %>%
      mutate(File_name = factor(File_name)) %>%
      mutate(Expedition = factor(Expedition)) %>%
      # filter(!`Curated ID` %in% c("Astyanax", #retirando as ids a nivel de genero
      #                           "Characidae",
      #                           "Cichlidae",
      #                           "Hoplias",
      #                           "Pimelodus")) %>%
      # filter(RRA >=0.01) %>% #retirando as ASVs "espurias", com abundancia menor que 0.01 ATUAL. Heron pediu para manter as reads com baixa abundancia
      # filter(!is.na(`Curated ID`)) %>% #retirar os NA
      ### Tile plot
      ggplot(aes(y = `Curated ID`,
                 x = Year,
                 fill = RRA,
                 )) +
      scale_x_continuous(breaks = 0:2100) +
      geom_tile() +
      labs(fill='Abundância \nrelativa (%)',
           x = "Ano",
           y = "Espécies") +
      geom_hline(yintercept = c(11.5)) +
      scale_fill_continuous(
        trans = "log10", # exibir a abundancia em escala logaritmica para favorecer a exibicao de baixo RRA
        breaks = c(0.001, 0.01, 0.1, 1, 10, 75), # definindo os valores que aparecem na escala
        type = "viridis") +
        ggtitle(label = "Espécies identificadas por ano (todas ids)") +
        theme(text=element_text(size = 14)) +
        guides(col = guide_legend(nrow = 6)) +
        theme(plot.title = element_text(hjust=0.5))
    
    ## Plotando
    ggsave(plot = tile_ano_all, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/tile_plots/tile_ano_all.pdf",
           device = "pdf", units = "cm", height = 25, width = 20, dpi = 600)
    }




# Alfa diversidade Ponto/Mes ----
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
    mutate(Mês = factor(Mês, levels = c("out21", "nov21"))) %>%
    filter(!Sample %in% c("LI1-neo-mi",
                          "L2_dez20")) %>% 
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
  
  ## Facetar por ponto 
  {
    print(alpha_plot_ponto <- alpha_tbl %>% 
            mutate(`Curated ID` = factor(`Curated ID`)) %>%
            ggplot(aes(x = Sample,
                       fill = Mês)) + ## preenchendo as cores por mes de coleta
            geom_bar(stat = "count", position = "stack") +
            guides(col = guide_legend(nrow = 6)) +
            xlab("Amostra") +
            ylab("Riqueza de espécies") +
            ggtitle(label = "Riqueza de espécies por amostra") +
            theme_bw(base_size = 16) +
            theme(legend.position = "bottom") +
            theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
            geom_text(stat='count', ## codigo que exibe dentro das barras a contagem
                      aes(label=..count..),
                      position = position_dodge(width = 0.9),
                      vjust = 1.4,
                      colour = "#ffffff", size = 6) +
            facet_grid2(cols = vars(Point), ##facetando por ponto
                        scales = "free_x",
                        axes = "all",
                        space = "free_x",
                        independent ='x') +
            scale_fill_manual(values = viridis::viridis(n=6)[c(2,5)]))
    
    ## Plotando
    ggsave(plot = alpha_plot_ponto, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/barplots/alpha_plot_ponto.pdf",
           device = "pdf", units = "cm", height = 15, width = 20, dpi = 600)
    }
  
  ## Facetar por Mes/ponto
  {
    print(alpha_plot_mes_ponto <- alpha_tbl %>% 
            mutate(`Curated ID` = factor(`Curated ID`)) %>% 
            ggplot(aes(x = Sample,
                       fill = Ponto)) + ## preenchendo as barras por mes de coleta
            geom_bar(stat = "count", position = "stack") +
            guides(col = guide_legend(nrow = 6)) +
            xlab("Amostra") +
            ylab("Riqueza de espécies") +
            ggtitle(label = "Riqueza de espécies por amostra") +
            theme_bw(base_size = 16) +
            theme(legend.position = "bottom") +
            theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
            geom_text(stat='count', ## codigo que exibe dentro das barras a contagem
                      aes(label=..count..),
                      position = position_dodge(width = 0.9),
                      vjust = 1.4,
                      colour = "#ffffff", size = 6) +
            facet_grid2(cols = vars(Mês), ##facetando por mes
                        scales = "free_x",
                        axes = "all",
                        space = "free_x",
                        independent ='x') +
            scale_fill_manual(values = viridis::viridis(n=10)[c(1,7,5,9)]))
    
    ## Plotando
    ggsave(plot = alpha_plot_mes_ponto, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/barplots/alpha_plot_mes_ponto.pdf",
           device = "pdf", units = "cm", height = 15, width = 20, dpi = 600)
  }
  
  ## Facetar por Mes
  {
    print(alpha_plot_mes_join <- alpha_tbl %>% 
            mutate(`Curated ID` = factor(`Curated ID`)) %>% 
            ggplot(aes(x = Mês,
                       fill = Ponto)) + ## preenchendo as barras por mes de coleta
            geom_bar(stat = "count", position = "stack") +
            guides(col = guide_legend(nrow = 6)) +
            xlab("Mês") +
            ylab("Riqueza de espécies") +
            ggtitle(label = "Riqueza de espécies por mês") +
            theme_bw(base_size = 16) +
            geom_text(stat = 'count', aes(label = ..count..), position = position_stack(vjust = 0.5), size = 4))
    
    ## Plotando
    ggsave(plot = alpha_plot_mes_join, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/barplots/alpha_plot_mes_join.pdf",
           device = "pdf", units = "cm", height = 15, width = 14, dpi = 600)
  }
  
  ## Ids por amostra
  
  {
    print(alpha_plot_id_sample <- alpha_tbl %>% 
      mutate(`Curated ID` = factor(`Curated ID`)) %>% 
      ggplot(aes(x = Sample,
                 y = `ID Abundance on sample (%)`,
                 group = Sample,
                 fill = `Curated ID`)) +
        ggtitle(label = "Espécies por amostra") +
        xlab("Amostra") +
        ylab("Proporção %") +
        guides(fill = guide_legend("Espécie")) + ## alterar o titulo da legenda
        # theme(legend.position = "bottom") +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_discrete(viridis::turbo(n = 32)))
    
    
    ## Plotando
    ggsave(plot = alpha_plot_id_sample, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/barplots/alpha_plot_id_sample.pdf",
           device = "pdf", units = "cm", height = 15, width = 28, dpi = 600)
  }
  
}

# Beta diversidade ----
{
  ## Calculo dos componentes de Beta diversidade de Jaccard
  {
    bd_j <- beta.div.comp(spp_sample_tbl_f[,c(5:36)], coef = "J", quant = T)
    bd_j$part
    bd_j$rich
    
   }
  
  ## Calculo dos componentes de Beta diversidade de Sorensen
  {
    bd_s <- beta.div.comp(spp_sample_tbl_f[,c(5:36)], coef = "S", quant = T)
    bd_s$part
    }
}

# NMDS Plots ----
{
  # NMDS Versao Heron ----
  
  # Criacao da tabela
  {
    # Obtencao dos dados
    {
      raw_results_NMDS <- read.csv(file = paste0(tbl_path,"/", "run_2_4_5_lagoa_ingleses_v2.csv"), sep = ",",
                                   header = TRUE,
                                   check.names = FALSE,
                                   na.strings=c("NA","NaN", "")) %>% 
        as_tibble()
      }
    
    # Criacao da lista com os possiveis nomes atribuidos as ASVs
    {
      raw_results_NMDS %>% colnames()
      raw_results_NMDS %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
    }
    
    # Agrupamento da ASVs que possuem os mesmos atributos abaixo
    {
      grouped_by_ID_NMDS <- raw_results_NMDS %>%
        select(c(
          "Sample",
          # "Run",
          # "Group",
          "Expedition",
          # "Coleta",
          "Year",
          "Sample.Name",
          "File_name",
          # "OTU",
          # "Read_origin",
          "Curated ID",
          # "final ID",
          # "Abundance",
          # "Relative abundance to all samples",
          "Relative abundance on sample",
          # "Sample.total.abundance",
          # "Type",
          "Point", # "Sub.point", # "Depth", # "Num.replicates", # "Obs", # "Primer", # "Quantidade.de.ovos.ou.larvas", 
          # "Kingdom", # "Phylum", # "Class", # "Order", # "Family", # "Genus", # "Species", # "Specimen", # "Basin", 
          # "exact.Genus", # "exact.Species", # "exact.GenSp", # "X1_subject.header", # "X1_subject", # "X1_indentity", 
          # "X1_length", # "X1_mismatches", # "X1_gaps", # "X1_query.start", # "X1_query.end", # "X1_subject.start", 
          # "X1_subject.end", # "X1_e.value", # "X1_bitscore", # "X1_qcovhsp", # "X2_subject.header", # "X2_subject", 
          # "X2_indentity", # "X2_length", # "X2_mismatches", # "X2_gaps", # "X2_query.start", # "X2_query.end", 
          # "X2_subject.start", # "X2_subject.end", # "X2_e.value", # "X2_bitscore", # "X2_qcovhsp", # "X3_subject.header", 
          # "X3_subject", # "X3_indentity", # "X3_length", # "X3_mismatches", # "X3_gaps", # "X3_query.start", 
          # "X3_query.end", # "X3_subject.start", # "X3_subject.end", # "X3_e.value", # "X3_bitscore", # "X3_qcovhsp", 
          # "Tag.pairs", # "Tag.FWD", # "Tag.REV", # "Control", # "Size..pb.", # "ASV.header", # "ASV..Sequence.", 
          # "Remove", # "Probable.bacteria", # "Abd..higher.than.in.control"
        )) %>% 
        group_by(Sample, `Curated ID`,  Expedition, Point, Sample.Name, File_name) %>% 
        summarize(`Sample` = unique(Sample),
                  `Curated ID` = unique(`Curated ID`),
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
      grouped_by_ID_NMDS$Sample %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat() ## valores possiveis de Amostra
      grouped_by_ID_NMDS$Expedition %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat() ## valores possiveis de Coleta
      grouped_by_ID_NMDS$Year %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat() ## valores possiveis de Ano
      grouped_by_ID_NMDS$Point %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat() ## valores possiveis de ponto
    }
    
    # Organizar as ordem das especies usando fatores
    
    # Sem os seguintes grupos
    {
      grouped_by_ID_NMDS <- grouped_by_ID_NMDS %>%
        mutate(`Curated ID` = factor(`Curated ID`,
                                     levels = rev(c(
                                       #"Actinopteri",
                                       "Acinocheirodon melanogramma",
                                       #"Astyanax",
                                       "Astyanax fasciatus",
                                       "Astyanax lacustris",
                                       "Brycon orthotaenia",
                                       "Bryconamericus stramineus",
                                       #"Characidae",
                                       #"Characidium",
                                       #"Characiformes",
                                       #"Cichla",
                                       #"Cichlidae",
                                       "Colossoma macropomum",
                                       "Coptodon zillii",
                                       "Eigenmannia virescens",
                                       "Gymnotus carapo",
                                       "Hemigrammus gracilis",
                                       "Hemigrammus marginatus",
                                       #"Hoplias",
                                       "Hoplias intermedius",
                                       "Hoplias malabaricus",
                                       "Hypomasticus steindachneri",
                                       "Leporellus vittatus",
                                       "Leporinus piau",
                                       "Leporinus reinhardti",
                                       "Leporinus taeniatus",
                                       "Megaleporinus elongatus",
                                       "Megaleporinus garmani",
                                       "Moenkhausia costae",
                                       "Myleus micans",
                                       "Oreochromis niloticus",
                                       "Orthospinus franciscensis",
                                       #"Pimelodus",
                                       "Pimelodus fur",
                                       "Pimelodus maculatus",
                                       "Pimelodus pohli",
                                       "Planaltina myersi",
                                       "Poecilia reticulata",
                                       "Prochilodus costatus",
                                       "Pseudoplatystoma corruscans",
                                       "Pygocentrus piraya",
                                       "Rhamdia quelen",
                                       "Salmo salar",
                                       "Serrasalmus brandtii",
                                       "Tilapia rendalli",
                                       "Wertheimeria maculata",
                                       #nao-peixes
                                       #"NA"
                                       "Cavia magna",
                                       #"Cutibacterium acnes",
                                       "Bos taurus",
                                       "Canis familiaris",
                                       "Didelphis albiventris (Gamba)",
                                       "Homo sapiens",
                                       "Hydrochaeris hydrochaeris (Capivara)",
                                       "Nannopterum brasilianus",
                                       "Oryctolagus cuniculus (Coelho-bravo)",
                                       "Progne chalybea (Andorinha-grande)",
                                       "Sus scrofa"
                                     ))))
    }
    
    # Fatorizar as variaveis para melhor plotagem
    
    {
      fact_NMDS <- grouped_by_ID_NMDS %>%
        mutate(Sample = factor(Sample, levels = c("L1_nov21",
                                                  "L1_out21",
                                                  "L2_dez20",
                                                  "L2_nov21",
                                                  "L2_out21",
                                                  "L3_nov21",
                                                  "L3_out21",
                                                  "L4_nov21",
                                                  "L4_out21",
                                                  "LI1-neo-mi")), 
               Expedition = factor(Expedition, levels = c("Dec/20",
                                                          "Nov_Dec/20",
                                                          "Nov/20",
                                                          "Nov/21",
                                                          "out/21")), 
               Year = factor(Year, levels = c("2020",
                                              "2021")),
               Point = factor(Point, levels = c("L1",
                                                "L2",
                                                "L3",
                                                "L4")))
    }  
  }
  
  # Criacao do plot NMDS
  {
    fact_NMDS$Sample %>% unique()
    fact_NMDS %>% colnames()
    
  }
  # 1- Preparar os dados para entrar no vegan
  {
    fact_NMDS_tbl <- fact_NMDS %>% 
      select(c(Sample, Point, `Curated ID`, Expedition, Year, RRA )) %>%
      group_by(Sample,`Curated ID`, Expedition, Year, Point) %>%
      summarise(RRA = sum(RRA)) %>%
      pivot_wider(c(Sample, Expedition, Year, Point), names_from = `Curated ID` ,values_from = RRA) %>%
      mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
      mutate("Sample number" = 0)  %>% 
      ungroup() %>%
      select(`Sample number`, 1:(ncol(.)-1))
  }
  #2- Associar um numero a cada amostra
  {
    # Etapa necessaria para os dados entrarem no pacote vegan
    for (sample in 1:nrow(fact_NMDS_tbl)) {
      fact_NMDS_tbl$`Sample number`[sample] <- sample
    }
    
    colnames(fact_NMDS_tbl)
    hist(colSums(fact_NMDS_tbl[,-c(1:5)]))
    hist(rowSums(fact_NMDS_tbl[,-c(1:5)]))
    fact_NMDS_tbl[,-c(1:5)]
    
    fact_NMDS_tbl %>% select(Sample, `Sample number`) %>% unique()
    
  }
  #3- Criar data.frame de contagem de especies: rownames sao Sample numbers
  {
    
    fact_NMDS_df <- fact_NMDS_tbl  %>% 
      select(sort(colnames(.))) %>%  ## reorganiza as colunas por ordem alfabetica
      select(-c( "Sample", ## retira as colunas Sample, Expedition, Year e Point nas primeiras posicoes
                 "Expedition",
                 "Year",
                 "Point")) %>% 
      as.data.frame() #salva a tabela como data.frame
    
  } 
  #4- name rows como Sample numbers e remover coluna
  {
    row.names(fact_NMDS_df) <- fact_NMDS_df$`Sample number` ## row.names como sample numbers
    
    fact_NMDS_df <- fact_NMDS_df %>% ## remover a coluna Sample numbers
      select(-c(`Sample number`))
    
    # Substituir o espaco no nome das especies por _
    colnames(fact_NMDS_df) <- colnames(fact_NMDS_df) %>%  
      str_replace_all(pattern = " ",replacement = "_") 
  }
  
  #5- Detrended Correspondence Analysis Plot 
  {
    
    ## Técnica de estatística Multivariada que permite clusterizar grupos de variaveis
    ## em conjuntos de dados muito grandes e com gradiente, permitindo observar tendencias
    
  }
  ## tentar deixar esse plot mais visivel
  {
    fact_NMDS_ps_ord <- decorana(veg = fact_NMDS_df)
    fact_NMDS_ps_ord %>% summary()
    fact_NMDS_ps_ord %>% str()
    fact_NMDS_ps_ord$cproj
    
    plot(fact_NMDS_ps_ord)
    plot(fact_NMDS_ps_ord,type = "p")
    plot(fact_NMDS_ps_ord,type = "c") 
    
    points(fact_NMDS_ps_ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
    text(fact_NMDS_ps_ord, display = "sites", cex=0.7, col="blue")
    text(fact_NMDS_ps_ord, display = "spec", cex=0.7, col="blue")
  }
  
  #6- NMDS analisys 
  {
    
    #6a- Calculate distances
    
    fact_NMDS_vg_dist <- vegdist(fact_NMDS_df, method="bray")
    scores(fact_NMDS_vg_dist)
    
    fact_NMDS_ps_ord %>% ncol()
    fact_NMDS_ps_ord <- fact_NMDS_df[,(colnames(fact_NMDS_df) %in% expected_sps)]
    
    fact_NMDS_ps_vegan_ord_meta <- metaMDS(veg = fact_NMDS_ps_ord, comm = fact_NMDS_df)
    # actually autotransform = FALSE doesn't seem to change the results
    plot(fact_NMDS_ps_vegan_ord_meta, type = "t")
    
    fact_NMDS_ps_vegan_ord_meta %>% str()
    fact_NMDS_ps_vegan_ord_meta
    plot(all_ps_vegan_ord_meta, type = "t")
    plot(all_ps_vegan_ord_meta, type = "p")
    
    fact_NMDS_ps_vegan_ord_meta$stress
    
    #6b- extract NMDS scores from results
    
    all_vegan_meta <- (vegan::scores(fact_NMDS_ps_vegan_ord_meta) %>% 
                         tidyr::as_tibble(rownames = "Sample number")) %>% 
      mutate(`Sample number` = as.numeric(`Sample number`))
    
    #7- bring NMDS scores to complete table
    
    all_vegan_meta_tbl <- left_join(x = unique(fact_NMDS_tbl[,c(1:4)]),
                                    y = all_vegan_meta, by = "Sample number") %>% 
      # mutate(Primer=factor(Primer,levels = c("NeoFish", "MiFish", "COI")),
      # Sample = as.factor(Sample)) %>% 
      select(-c("Sample number"))
    
    nmds_PLOT_ord <- ggord(fact_NMDS_ps_vegan_ord_meta, 
                           grp_in = fact_NMDS_tbl$Year, 
                           ellipse = F,
                           size = 10,
                           arrow = 0.5, veccol = "dark grey",
                           txt = 3,
                           repel = T,
                           max.overlaps = 55
                           # ,
                           # # facet=T,
                           # cols = viridis::viridis(option = "turbo",n = nrow(all_IDs_NMDS_df), alpha = 1)
    )+
      annotate(geom = "text",
               x=c(0.2),
               y=c(0.20),
               label=c(paste0("Stress: ",format(round(fact_NMDS_ps_vegan_ord_meta$stress,4)))),
               size=5) +
      scale_colour_manual(name = "Samples", values = viridis::viridis(option = "turbo",n = 4, alpha = 1))
    nmds_PLOT_ord
    
    
    ### Tirar as especies com ids muito genericas e fazer esse grafico apenas
    ### para 2021 separando os clusters por pontos coletados (L1, L2, L3 e L4)
    
    
    nmds_PLOT <- all_vegan_meta_tbl %>% 
      # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>% 
      ggplot(aes(x = NMDS1,y = NMDS2, col = Sample,
                 shape = Year,
                 label = Sample,
                 Group = Year))+
      geom_point(size = 11)+
      theme_linedraw(base_size = 18) +
      theme(legend.position="bottom") +
      coord_fixed(ratio = 1) +
      scale_color_manual(values = viridis::viridis(option = "turbo",n = 30, alpha = 1, begin = 0, end = 1, direction = 1)
      ) +
      annotate(geom = "text",
               x=c(0.30),y=c(-0.3),label=c(paste0("Stress: ",fact_NMDS_ps_vegan_ord_meta$stress)),size=5)  
    
    nmds_PLOT
    
    
    ggsave(file = paste0(figs_path,"/",prjct_radical,"_NMDS.pdf"),
           plot = nmds_PLOT,
           device = "png",
           width = 31,
           height =20,
           units = "cm",
           dpi = 300)
    
    
  }
  
  # NMDS Versao Ecological Applications in R ----  
  # Esta e uma abordagem heuristica, nao estamos testando nenhuma hipotese ainda!
  # O unico objetivo e tentar observar padroes em nossos dados!
  {
    # Transformar os dados da comunidade na matriz spp_sample_tbl para entrar no NMDS
    spp_sample_hel <- decostand(spp_sample_tbl[c,4:35], method = "hellinger")
    # Criando o NMDS
    nmds1 <- metaMDS(spp_sample_hel, autotransform = FALSE) # Esta dando uma Warning Message. 
    # stress is (nearly) zero: you may have insufficient data
    # Isso seria um problema?
    # Plotando no vegan
    nmds1_vplot <- ordiplot(nmds1, type = "t") # type = "t" permite ver o nome dos pontos e das especies
    
    # Plotando no ggvegan
    nmds1_gvplot <- autoplot(nmds1) # ainda nao e um plot ideal. Os proximos passos vao deixar o plot mais publicavel
    
    # Plotando no ggplot
    fort <- fortify(nmds1) # "fortificando" os dados para eles entrarem no ggplot
    nmds1_ggplot <- ggplot() +
      geom_point(data = subset(fort, Score == 'sites'), # criacao do scatterplot com os pontos (como visto acima)
                 mapping = aes(x = NMDS1,
                               y = NMDS2),
                 colour = "black",
                 alpha = 0.5) +
      geom_segment(data = subset(fort, Score == 'species'), # criacao dos vetores (setas)
                   mapping = aes(x = 0,
                                 y = 0,
                                 xend = NMDS1,
                                 yend = NMDS2),
                   arrow = arrow(length = unit(0.015, "npc"),
                                 type = "closed"),
                   colour = "darkgrey",
                   size = 0.8) +
      geom_text(data = subset(fort, Score == 'species'), # colocando o nome dos pontos e especies 
                mapping = aes(label = Label, x = NMDS1 * 1.1,
                              y = NMDS2 * 1.1)) + # multiplicando as coordenadas para afastar os nomes do fim das setas
      geom_abline(intercept = 0, slope =  0, linetype = "dashed", size = 0.8, colour = "gray") +
      geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.8, colour = "gray") +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"))
    
    ## Plotando
    ggsave(plot = nmds1_ggplot, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/NDMS/nmds1_ggplot.pdf",
           device = "pdf", units = "cm", height = 15, width = 20, dpi = 600)
    
    # Plotando no ggord
    nmds1_ggord <- ggord(nmds1,
                         grp_in = sort(row.names(spp_sample_tbl)),
                         vectyp = "dotted",
                         parse = F,
                         ellipse = F,
                         size = 8,
                         arrow = 1, 
                         veccol = "dark grey",
                         txt = 3,
                         repel = T,
                         veclsz = 0.75,
    ) +
      annotate(geom = "text",
               x = c(-1.25),
               y = c(1),
               label = c(paste0("Stress: ",
                                format(round(nmds1$stress,4)))),
               size = 5)
    
    ggsave(plot = nmds1_ggord, filename = "/home/gabriel/projetos/lagoa_ingleses/results/figuras/agosto/NDMS/nmds1_ggord.pdf",
           device = "pdf", units = "cm", height = 15, width = 20, dpi = 600)
  }
  # Minha Versao ----
  {
    
  } 
}
# PCA Plots ----
{
  pca1 <- rda(spp_sample_tbl)
  pca1
  
  # o BiodiversityR nao roda em servidores devido a problemas com o X11 
  # por isso irei rodar na minha maquina domestica, e para isso vou baixar a matriz que possui as amostras e abundancias
  write.csv(spp_sample_tbl, "/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5/ssp_sample_tbl", row.names = TRUE)
  
  }
# Mapa da Lagoa dos Ingleses e os respectivos pontos ----
  
  # Localizacao da Lagoa dos Ingleses: 
  lagoa_ingleses <- c(-20.17874819739011, -43.96130043170407)
  
  # Localizacao dos pontos:
  {
    fundacao_4 <- c(-20.163437545269186, -43.955100550084005)
    ponte_3 <- c(-20.173468237193315, -43.950422777741935)
    barragem_2 <- c(-20.177738135335836, -43.94299842263384)
    prainha_1 <- c(-20.17910770033768, -43.95947791566818)
    pontos <- data.frame(prainha_1, barragem_2, ponte_3, fundacao_4)
    }
      
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

  