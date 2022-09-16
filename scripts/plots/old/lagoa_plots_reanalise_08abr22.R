library(dplyr)
library(ggplot2)

#Obtencao dos dados
{
  raw_results_tbl <- read.csv("~/projetos/lagoa_ingleses/tabelas/raw/08mar22/v3_ed_li_22mar22-todas_info_da_analise_2022-04-08.csv",
                              sep = ";", header = T, dec = ",")
}

#Criacao da lista com os possiveis nomes atribuidos as ASVs
{
  raw_results_tbl %>% colnames()
  raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
}

#Agrupamento da ASVs que possuem os mesmos atributos abaixo
{
  grouped_by_ID_tbl <- raw_results_tbl %>%
    select(c(
      # ASV  Sequence ",
      "Sample",
      # "Abundance",
      # "Run",
      "Coleta",
      "Ano",
      "Sample.Name",
      "File_name",
      # "Type",
      "Point",
      # "Filter", # "Num replicates", # "Obs", # "Primer", # "Tag pairs", # "Tag FWD", # "Tag REV", # "Control", # "Kingdom", # "Phylum", # "Class", # "Order", # "Family", # "Genus", # "Species", # "Specimen", # "Basin", # "Read_origin", # "exact Genus",
      # "exact Species", # "X1_res", # "X1_subject header", # "X1_query", # "X1_subject", # "X1_indentity", # "X1_length", # "X1_mismatches", # "X1_gaps", # "X1_query start", # "X1_query end", # "X1_subject start", # "X1_subject end", # "X1_e value", # "X1_bitscore", # "X1_qcovhsp",
      # "X2_res", # "X2_subject header", # "X2_query", # "X2_subject", # "X2_indentity", # "X2_length", # "X2_mismatches", # "X2_gaps", # "X2_query start", # "X2_query end", # "X2_subject.start", # "X2_subject.end",# "X2_e.value", # "X2_bitscore", # "X2_qcovhsp", # "X3_res", # "X3_subject.header",
      # "X3_query", # "X3_subject", # "X3_indentity", # "X3_length", # "X3_mismatches", # "X3_gaps", # "X3_query.start", # "X3_query.end", # "X3_subject.start", # "X3_subject.end", # "X3_e.value", # "X3_bitscore", # "X3_qcovhsp", # "Relative.abundance.to.all.samples", 
      "Relative.abundance.on.sample", # "Sample total abundance", #"exact GenSp", #"final ID",
      "Curated.ID", # "Size  pb ", # "ASV header", # "Remove
    )) %>% group_by(Sample, `Curated.ID`, Coleta, Point , `Sample.Name`, File_name) %>%
    summarize(`Sample` = unique(Sample),
              `Curated.ID` = unique(`Curated.ID`),
              `Coleta` = unique(Coleta),
              `Ano` = unique(Ano),
              `Point` = unique(Point),
              `Sample.Name` = unique(`Sample.Name`),
              `File_name` = unique(File_name),
              `RRA` = sum(`Relative.abundance.on.sample`)) %>%
    ungroup()
}

#Trocando os nomes dos animais não-peixes
{
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Bos taurus"] <- "Bos taurus (Boi)"
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Canis familiaris"] <- "Canis familiaris (Cachorro doméstico)"
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Cavia magna"] <- "Cavia magna (Preá)"
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Hydrochaeris hydrochaeris"] <- "Hydrochaeris hydrochaeris (Capivara)"
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Nannopterum brasilianus"] <- "Nannopterum brasilianus (Biguá)"
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Progne chalybea"] <- "Progne chalybea (Andorinha-grande)"
  grouped_by_ID_tbl$`Curated.ID`[grouped_by_ID_tbl$`Curated.ID` == "Sus scrofa"] <- "Sus scrofa (Javali)"
}

#organizar as especies
{
  grouped_by_ID_tbl$`Curated.ID` %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  grouped_by_ID_tbl$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
}

#organizar as ordem das especies usando fatores

#com todas as identificacoes
{
  grouped_by_ID_tbl <- grouped_by_ID_tbl %>%
    mutate(Curated.ID = factor(Curated.ID,
                               levels = rev(c( #"",
                                 "Astyanax",
                                 "Astyanax fasciatus",
                                 "Astyanax lacustris",
                                 "Brycon orthotaenia",
                                 "Bryconamericus stramineus",
                                 "Characidae",
                                 "Characidium",
                                 "Cichlidae",
                                 "Colossoma macropomum",
                                 "Eigenmannia virescens",
                                 "Gymnotus carapo",
                                 "Hemigrammus gracilis",
                                 "Hemigrammus marginatus",
                                 "Hoplias",
                                 "Hoplias intermedius",
                                 "Hoplias malabaricus",
                                 "Hoplias sp",
                                 "Hypomasticus steindachneri",
                                 "Leporellus vittatus",
                                 "Leporinus piau",
                                 "Leporinus reinhardti",
                                 "Leporinus taeniatus",
                                 "Loricariidae",
                                 "Megaleporinus elongatus",
                                 "Megaleporinus garmani",
                                 "Moenkhausia costae",
                                 "Myleus micans",
                                 "Oreochromis niloticus",
                                 "Oreochromis niloticus ",
                                 "Orthospinus franciscensis",
                                 "Phalloceros sp",
                                 "Pimelodus",
                                 "Pimelodus fur",
                                 "Pimelodus maculatus",
                                 "Pimelodus pohli",
                                 "Planaltina myersi",
                                 "Prochilodus costatus",
                                 "Pseudoplatystoma corruscans",
                                 "Pygocentrus piraya",
                                 "Rhamdia quelen",
                                 "Salmo salar",
                                 "Serrasalmus brandtii",
                                 "Tilapia rendalli",
                                 "Wertheimeria maculata",
                                 #não-peixes
                                 "Homo sapiens",
                                 "Bos taurus (Boi)",
                                 "Canis familiaris (Cachorro doméstico)",
                                 "Cavia magna (Preá)",
                                 "Hydrochaeris hydrochaeris (Capivara)",
                                 "Nannopterum brasilianus (Biguá)",
                                 "Progne chalybea (Andorinha-grande)",
                                 "Sus scrofa (Javali)"
                               ))))
}

#Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5

#exibindo  por ponto, por mês
{
  grouped_by_ID_tbl %>%
    #transformando variaveis categoricas em fatores com niveis
    mutate(Coleta = factor(Coleta)) %>%
    mutate(Sample = factor(Sample,
                           levels = c("L1_nov_dec_20_mi", 
                                      "L1_nov_dec_20_neo", 
                                      "L2_dez20", 
                                      "L2_nov20",
                                      "L1_nov21", 
                                      "L1_out21", 
                                      "L2_nov21", 
                                      "L2_out21", 
                                      "L3_nov21", 
                                      "L3_out21", 
                                      "L4_nov21", 
                                      "L4_out21" ))) %>% 
    mutate(Coleta = factor(Coleta, 
                           levels = c("Nov_Dec/20",
                                      "Nov/20",
                                      "Dec/20",
                                      "out/21",
                                      "Nov/21"
                                      ))) %>%
    mutate(Point = factor(Point)) %>%
    mutate(File_name = factor(File_name)) %>%
    mutate(Coleta = factor(Coleta)) %>%
    
    #filtrando apenas o que vc quer mostras
    filter(!Curated.ID %in% c("Astyanax",
                              "Characidae",
                              "Characidium",
                              "Cichlidae",
                              "Hoplias",
                              "Loricariidae",
                              "Pimelodus",
                              "Salmo salar")) %>%
    
    #retirar os NA
    filter(!is.na(Curated.ID)) %>%
    
    #retirando as ASVs "espurias", com abundancia menor que 0.01%
    filter(RRA >=0.01) %>%
    
    #Tile plot
    ggplot(aes(y = Curated.ID,
               x = Point,
               fill = RRA,
               # col = Coleta,
    )) +
    geom_tile() +
    # scale_color_manual()+
    facet_grid(~Coleta,
               scales = "free_x",
               space = "free_x") +
    # facet_grid(~Point,
    #            scales = "free_x",
    #            space = "free_x") +
    #facet_grid(~Year,
    #            scales = "free_x",
    #            space = "free_x") +
    labs(fill='Relative Read\nAbundance (%)',
         x = "Amostras",
         y= "Espécies") +
    geom_hline(yintercept = c(6.5)) +
      viridis::scale_fill_viridis(
      alpha = 1,
      begin = 0.75,
      end = 0,
      direction = 1,
      discrete = FALSE,
      option = "D"
    )
    theme(text=element_text(size = 10, face = "bold"))
  
}

#Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5
#exibindo  por ano apenas
{
  grouped_by_ID_tbl %>%
    
    #transformando variaveis categoricas em fatores com niveis
    mutate(Coleta = factor(Coleta)) %>%
    mutate(Sample = factor(Sample,
                           levels = c("L1_nov_dec_20_mi", 
                                      "L1_nov_dec_20_neo", 
                                      "L2_dez20", 
                                      "L2_nov20",
                                      "L1_nov21", 
                                      "L1_out21", 
                                      "L2_nov21", 
                                      "L2_out21", 
                                      "L3_nov21", 
                                      "L3_out21", 
                                      "L4_nov21", 
                                      "L4_out21" 
                                      
                           ))) %>% mutate(Coleta = factor(Coleta,
                                                              levels = c("Nov_Dec/20",
                                                                         "Nov/20",
                                                                         "Dec/20",
                                                                         "out/21",
                                                                         "Nov/21"
                                                              ))) %>%
    mutate(Point = factor(Point)) %>%
    mutate(File_name = factor(File_name)) %>%
    mutate(Coleta = factor(Coleta)) %>%
    
    #filtrando apenas o que vc quer mostras
    # filter(!Curated.ID %in% c("Astyanax",
    #                           "Characidae",
    #                           "Characidium",
    #                           "Cichlidae",
    #                           "Hoplias",
    #                           "Loricariidae",
    #                           "Pimelodus",
    #                           "Salmo salar")) %>%
    
    #retirando as ASVs "espurias", com abundancia menor que 0.01
    # filter(RRA >=0.01) %>%
    
    #retirar os NA
    filter(!is.na(Curated.ID)) %>%
    
    #Tile plot
    ggplot(aes(y = Curated.ID,
               x = Ano,
               fill = RRA,
               # col = Coleta,
    )) +
    scale_x_continuous(breaks = 0:2100) +
    
    # geom_tile() +
    geom_tile() +
    # scale_color_manual()+
    #facet_grid(~Coleta,
    #           scales = "free_x",
    #           space = "free_x") +
    # facet_grid(~Point,
    #            scales = "free_x",
    #            space = "free_x") +
    #facet_grid(~Year,
    #            scales = "free_x",
    #            space = "free_x") +
    labs(fill='Relative Read\nAbundance (%)',
         x = "Amostras",
         y= "Espécies") +
    geom_hline(yintercept = c(8.5)) +
    viridis::scale_fill_viridis(
      alpha = 1,
      begin = 0.75,
      end = 0,
      direction = 1,
      discrete = FALSE,
      option = "D"
    )
    theme(text=element_text(size = 10, face = "bold"))
}

