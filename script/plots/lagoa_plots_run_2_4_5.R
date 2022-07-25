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
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(phyloseq)
  library(Biostrings)
  library(Matrix)
  library(ShortRead)
  library(DECIPHER)
  library(future)
  library(vegan)
}

# Caminhos 
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  
  results_path <- paste0(prjct_path,"/results")
  
  figs_path <- paste0(results_path,"/figuras")
  
  tbl_path <- paste0(prjct_path,"/tabelas/raw/run_2_4_5")
  
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

# Obtencao dos dados
{
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v2.csv"), sep = ",")
}

# Tile Plots ----

# Criacao da lista com os possiveis nomes atribuidos as ASVs
{
raw_results_tbl %>% colnames()
raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
}

# Agrupamento da ASVs que possuem os mesmos atributos abaixo
{
grouped_by_ID_tbl <- raw_results_tbl %>%
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
    "Curated.ID",
    # "final.ID",
    # "Abundance",
    # "Relative.abundance.to.all.samples",
    "Relative.abundance.on.sample",
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
    group_by(Sample, Curated.ID, Expedition, Point, Sample.Name, File_name) %>% 
    summarize(`Sample` = unique(Sample),
              `Curated.ID` = unique(Curated.ID),
              `Expedition` = unique(Expedition),
              `Year` = unique(Year),
              `Point` = unique(Point),
              `Sample.Name` = unique(Sample.Name),
              `File_name` = unique(File_name),
              `RRA` = sum(Relative.abundance.on.sample)) %>%
    ungroup()
}

# Organizar as especies
{
grouped_by_ID_tbl$Curated.ID %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
grouped_by_ID_tbl$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
}

# Organizar as ordem das especies usando fatores

  #com todas as identificacoes
{
grouped_by_ID_tbl <- grouped_by_ID_tbl %>%
  mutate(Curated.ID = factor(Curated.ID,
                             levels = rev(c(
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
                                        "Serrasalmus brandtii",
                                        "Tilapia rendalli",
                                        "Wertheimeria maculata",
                                        #não-peixes
                                        "Cavia magna",
                                        "Salmo salar",
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

                                        ))))
}

# Sem os grupos
{
  grouped_by_ID_tbl <- grouped_by_ID_tbl %>%
    mutate(Curated.ID = factor(Curated.ID,
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
                                 #"Salmo salar",
                                 "Serrasalmus brandtii",
                                 "Tilapia rendalli",
                                 "Wertheimeria maculata",
                                 #não-peixes
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
                                 
                               ))))
}

# Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5
# exibindo  por ponto, por mês
{
grouped_by_ID_tbl %>%
  # transformando variaveis categoricas em fatores com niveis
  mutate(Expedition = factor(Expedition)) %>%
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
                                    
  ))) %>% mutate(Expedition = factor(Expedition,
                                 levels = c("Nov_Dec/20",
                                            "Nov/20",
                                            "Dec/20",
                                            "out/21",
                                            "Nov/21"
 ))) %>%
  mutate(Point = factor(Point)) %>%
  mutate(File_name = factor(File_name)) %>%
  mutate(Expedition = factor(Expedition)) %>%
 
  # filtrando apenas o que vc quer mostras
  filter(!Curated.ID %in% c("Astyanax",
                            "Characidae",
                            "Cichlidae",
                            "Hoplias",
                            "Pimelodus")) %>%
  
  # retirando as ASVs "espurias", com abundancia menor que 0.01
  filter(RRA >=0.01) %>%
  
  # retirar os NA
  filter(!is.na(Curated.ID)) %>%
    
    #Tile plot
   
     ggplot(aes(y = Curated.ID,
             x = Point,
             fill = RRA,
             # col = Expedition,
             )) +
  # geom_tile() +
  geom_tile() +
  # scale_color_manual()+
  facet_grid(~Expedition,
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
  geom_hline(yintercept = c(7.5)) +
  scale_fill_continuous(type = "viridis") +
       theme(text=element_text(size = 10, face = "bold"))
    
}

# Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5
# Exibindo  por ano apenas
{
  grouped_by_ID_tbl %>%
    #transformando variaveis categoricas em fatores com niveis
    mutate(Expedition = factor(Expedition)) %>%
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
                                      
                           ))) %>% mutate(Expedition = factor(Expedition,
                                                              levels = c("Nov_Dec/20",
                                                                         "Nov/20",
                                                                         "Dec/20",
                                                                         "out/21",
                                                                         "Nov/21"
                                                              ))) %>%
    mutate(Point = factor(Point)) %>%
    mutate(File_name = factor(File_name)) %>%
    mutate(Expedition = factor(Expedition)) %>%
    
    #filtrando apenas o que vc quer mostras
    filter(!Curated.ID %in% c("Astyanax",
                              "Characidae",
                              "Cichlidae",
                              "Hoplias",
                              "Pimelodus")) %>%
    
    #retirando as ASVs "espurias", com abundancia menor que 0.01
    filter(RRA >=0.01) %>%
    
    #retirar os NA
    filter(!is.na(Curated.ID)) %>%
    
    #Tile plot
    
    ggplot(aes(y = Curated.ID,
               x = Year,
               fill = RRA,
               # col = Expedition,
    )) +
    scale_x_continuous(breaks = 0:2100) +
    # geom_tile() +
    geom_tile() +
    # scale_color_manual()+
    #facet_grid(~Expedition,
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
    geom_hline(yintercept = c(7.5)) +
    scale_fill_continuous(type = "viridis") +
    theme(text=element_text(size = 10, face = "bold"))
}

# Criacao do NMDS ----

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
      # "final.ID",
      # "Abundance",
      # "Relative.abundance.to.all.samples",
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
    grouped_by_ID_NMDS$Sample %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  grouped_by_ID_NMDS$Expedition %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  grouped_by_ID_NMDS$Year %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  grouped_by_ID_NMDS$Point %>% unique() %>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
  }

# Fatorizar as variaveis para melhor plotagem

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

# Criacao do plot NMDS

fact_NMDS$Sample %>% unique()
fact_NMDS %>% colnames()

# 1- Preparar os dados para entrar no vegan ----

fact_NMDS_tbl <- fact_NMDS %>% 
  select(c(Sample, Point, `Curated ID`, Expedition, Year, RRA )) %>%
  group_by(Sample,`Curated ID`, Expedition, Year, Point) %>%
  summarise(RRA = sum(RRA)) %>%
  pivot_wider(c(Sample, Expedition, Year, Point), names_from = `Curated ID` ,values_from = RRA) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>%
  mutate("Sample number" = 0)  %>% 
  ungroup() %>%
  select(`Sample number`, 1:(ncol(.)-1))


#2- Asscoiar um numero a cada amostra ----

# Etapa necessaria para os dados entrarem no pacote vegan
for (sample in 1:nrow(fact_NMDS_tbl)) {
  fact_NMDS_tbl$`Sample number`[sample] <- sample
  }

colnames(fact_NMDS_tbl)
hist(colSums(fact_NMDS_tbl[,-c(1:5)]))
hist(rowSums(fact_NMDS_tbl[,-c(1:5)]))
fact_NMDS_tbl[,-c(1:5)]

fact_NMDS_tbl %>% select(Sample, `Sample number`) %>% unique()

#3- Criar data.frame de contagem de especies: rownames sao Sample numbers ----

fact_NMDS_df <- fact_NMDS_tbl  %>% 
  select(base::sort(colnames(.))) %>% 
  relocate(c("Sample",
             "Expedition",
             "Year",
             "Point")) %>% 
  as.data.frame() 

#4- name rows como Sample numbers e remover coluna ----
row.names(fact_NMDS_df) <- fact_NMDS_df$`Sample number`

fact_NMDS_df <- fact_NMDS_df %>% 
  select(-c(`Sample number`))


colnames(fact_NMDS_df)[5:64] <- colnames(fact_NMDS_df)[5:64] %>% 
  str_replace_all(pattern = " ",replacement = "_") %>% 
  str_replace_all(pattern = "Gamb\xe1", replacement = "Gamba") %>% 
  str_replace_all(pattern = "Sus_scrofa_\\*",replacement = "Sus_scrofa")
  


#5- ----

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



#6- NMDS analisys ----

#6a- Calculate distances

fact_NMDS_vg_dist <- vegdist(fact_NMDS_df[5:64], method="bray")

vegan::scores(fact_NMDS_vg_dist)
# 
fact_NMDS_ps_ord %>% ncol()
fact_NMDS_ps_ord <- fact_NMDS_df[,(colnames(fact_NMDS_df) %in% expected_sps)]

# 
# fact_NMDS_ps_ord %>% ncol()
# fact_NMDS_vg_dist <- vegdist(fact_NMDS_df, method="bray")
# 
# fact_NMDS_ps_ord <- decorana(veg = all_IDs_NMDS_df)
# 
# fact_NMDS_ps_ord %>% summary()
# 
# fact_NMDS_ps_ord %>% str()
# 
# fact_NMDS_ps_ord$cproj
# 
# 
# plot(fact_NMDS_ps_ord)
# plot(fact_NMDS_ps_ord,type = "p")
# plot(fact_NMDS_ps_ord,type = "c") 
# vegan::scores(fact_NMDS_ps_ord)

fact_NMDS_ps_vegan_ord_meta <- metaMDS(veg = fact_NMDS_ps_ord[5:64], comm = fact_NMDS_df[5:64])
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
# all_vegan_meta <- as.data.frame(vegan::scores(fact_NMDS_ps_vegan_ord_meta))

#Using the scores function from vegan to extract the site scores and convert to a data.frame

# all_vegan_meta$`Sample number` <- rownames(all_vegan_meta) %>% as.numeric()  

# all_vegan_meta %>% left_join()# create a column of site names, from the rownames of data.scores

# all_vegan_meta <- all_vegan_meta  %>% as_tibble() # create a column of site names, from the rownames of data.scores

#7- bring NMDS scores to complete table

all_vegan_meta_tbl <- left_join(x = unique(fact_NMDS_tbl[,c(1:4)]),
                                y = all_vegan_meta, by = "Sample number") %>% 
  # mutate(Primer=factor(Primer,levels = c("NeoFish", "MiFish", "COI")),
         # Sample = as.factor(Sample)) %>% 
  select(-c("Sample number"))





library(factoextra)
library(ggforce)
library(ggord)




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
# + 
#   # labs(title='NEW LEGEND TITLE') +
#   # labs(fill='NEW LEGEND TITLE') +
#   # labs(colour='NEW LEGEND TITLE') + 
#   labs(shape='NEW LEGEND TITLE') 
#   

nmds_PLOT_ord


### Tirar as espécies com ids muito genéricas e fazer esse gráfico apenas
### para 2021 separando os clusters por pontos coletados (L1, L2, L3 e L4)


nmds_PLOT <- all_vegan_meta_tbl %>% 
  # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>% 
  ggplot(aes(x = NMDS1,y = NMDS2, col = Sample,
             shape = Year,
             label = Sample,
             Group = Year))+
  # stat_ellipse()+ 
  geom_point(size = 11)+
  theme_linedraw(base_size = 18) +
  theme(legend.position="bottom") +
  coord_fixed(ratio = 1) +
  # ggrepel::geom_label_repel(label.size = 0.8,size = 3,min.segment.length = 2) +
  # ggrepel::geom_text_repel(col="black",size = 3,min.segment.length = 2) +
  # scale_shape_manual() %>% 
  scale_color_manual(values = viridis::viridis(option = "turbo",n = 30, alpha = 1, begin = 0, end = 1, direction = 1)
                     # labels = c("NeoFish", "MiFish", "COI"), values = alpha(colour = colors_norm )
  ) +
  annotate(geom = "text",
           x=c(0.30),y=c(-0.3),label=c(paste0("Stress: ",fact_NMDS_ps_vegan_ord_meta$stress)),size=5)  

# +
#   
#   
#   # ADD ggforce's ellipses
#   ggforce::geom_mark_ellipse(inherit.aes = FALSE,
#                              aes(x = NMDS1,y = NMDS2,
#                                  group=`Sample`,
#                                  label=`Sample`),
#                              n = 100,
#                              expand = 0.07,
#                              label.fontsize = 20,con.cap = 0.1) 

nmds_PLOT
# 
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_NMDS.pdf",
#      plot = nmds_PLOT,
#      device = "pdf",
#      width = 14,
#      height =10,
#      dpi = 600)
# 
# ggsave(file = "~/prjcts/fish_eDNA/sfjq/results/figs/SFJQ_NMDS.svg",
#      plot = nmds_PLOT,
#      device = "svg",
#      width = 14,
#      height =10,
#      dpi = 600)


ggsave(file = paste0(figs_path,"/",prjct_radical,"_NMDS.pdf"),
       plot = nmds_PLOT,
       device = "png",
       width = 31,
       height =20,
       units = "cm",
       dpi = 300)


