
---- ## Escolhendo as maiores seqs entre as ASVs de cada espécie ## ----

# Carregando bibliotecas ----
{
  library(base)
  library(stringr)
  library(tibble)
  library(tidyr)
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

# Escolher maiores ASVs de cada especie ----

bigger_seqs <- raw_results_tbl %>% 
  select(`Curated ID`,
         `ASV header`,
         `ASV (Sequence)`,
         Expedition) %>% 
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
                              "Sus scrofa")) %>%
  filter(Expedition %in% c("out/21", # deixando apenas as amostras de 2021
                           "Nov/21")) %>% 
  group_by(`Curated ID`) %>% 
  slice(which.max(nchar(`ASV (Sequence)`)))

  
  

