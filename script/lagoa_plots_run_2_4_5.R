library(dplyr)
library(ggplot2)

#Obtencao dos dados
{
  raw_results_tbl <- read.csv("~/projetos/lagoa_ingleses/tabelas/raw/run_2_4_5_lagoa_ingleses_csv.csv",sep = ",")
}

#Criacao da lista com os possiveis nomes atribuidos as ASVs
{
raw_results_tbl %>% colnames()
raw_results_tbl %>% colnames() %>% paste0(collapse = '",\n"') %>% cat()
}

#Agrapamento da ASVs que possuem os mesmos atributos abaixo
{
grouped_by_ID_tbl <- raw_results_tbl %>%
  select(c(
    "Sample",
    # "Run",
    # "Group",
    "Expedition",
    # "Coleta",
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
    "Point", # "Sub.point", # "Depth", # "Num.replicates", # "Obs", # "Primer", # "Quantidade.de.ovos.ou.larvas", # "Kingdom", # "Phylum", # "Class", # "Order", # "Family", # "Genus", # "Species", # "Specimen", # "Basin", # "exact.Genus", # "exact.Species", # "exact.GenSp", # "X1_subject.header", # "X1_subject", # "X1_indentity", # "X1_length", # "X1_mismatches", # "X1_gaps", # "X1_query.start", # "X1_query.end", # "X1_subject.start", # "X1_subject.end", # "X1_e.value", # "X1_bitscore", # "X1_qcovhsp", # "X2_subject.header", # "X2_subject", # "X2_indentity", # "X2_length", # "X2_mismatches", # "X2_gaps", # "X2_query.start", # "X2_query.end", # "X2_subject.start", # "X2_subject.end", # "X2_e.value", # "X2_bitscore", # "X2_qcovhsp", # "X3_subject.header", # "X3_subject", # "X3_indentity", # "X3_length", # "X3_mismatches", # "X3_gaps", # "X3_query.start", # "X3_query.end", # "X3_subject.start", # "X3_subject.end", # "X3_e.value", # "X3_bitscore", # "X3_qcovhsp", # "Tag.pairs", # "Tag.FWD", # "Tag.REV", # "Control", # "Size..pb.", # "ASV.header", # "ASV..Sequence.", # "Remove", # "Probable.bacteria", # "Abd..higher.than.in.control"
  )) %>% group_by(Sample,Curated.ID,Expedition, Point ,Sample.Name, File_name) %>%
  summarize(`Sample` = unique(Sample),
            `Curated.ID` = unique(Curated.ID),
            `Expedition` = unique(Expedition),
            `Point` = unique(Point),
            `Sample.Name` = unique(Sample.Name),
            `File_name` = unique(File_name),
            `RRA` = sum(Relative.abundance.on.sample)) %>%
  ungroup()
}


#organizar as especies
{
grouped_by_ID_tbl$Curated.ID %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
grouped_by_ID_tbl$Sample %>% unique()%>% sort() %>%  paste0(collapse = '",\n"') %>% cat()
}

#organizar as ordem das especies usando fatores
{
grouped_by_ID_tbl <- grouped_by_ID_tbl %>%
  mutate(Curated.ID = factor(Curated.ID,
                             levels = rev(c(
                                        "Acinocheirodon melanogramma",
                                        "Actinopteri",
                                        "Astyanax",
                                        "Astyanax fasciatus",
                                        "Astyanax lacustris",
                                        "Brycon orthotaenia",
                                        "Bryconamericus stramineus",
                                        "Cavia magna",
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
                                        #não-peixes
                                        "NA",
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

#Criacao do Tile Plot das amostras da Lagoa dos Ingleses sequenciadas nas corridas 2, 4 e 5
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
  labs(fill='Relative Read\nAbundance (%)',
       x = "Amostras",
       y= "Espécies") +
  geom_hline(yintercept = c(9.5)) +
  scale_fill_continuous(type = "viridis")
}