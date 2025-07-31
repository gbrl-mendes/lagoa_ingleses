---
  title: "Extraindo ASVs da raw_results_tbl para um FASTA"
author: "Gabriel Mendes"
date: "01/2024"
---
  
  ## Atencao! Esse script funciona a partir da extracao das seqs de ASVs contidas no script 102023_tile-plots, 
  ## onde estao as analises realizadas com os dados dos sequenciamoentos 2, 4, 5 e EM156.
  ## Para usar com outros dados sera necessario adaptar.
  
  ## Carregando bibliotecas ----
{
  library(Biostrings)
  library(DECIPHER)
  library(ShortRead)
  library(stringr)
  library(tidyverse)
  library(ggtree)
  library(dplyr)
}

## Criando a tabela com as seqs e seus respectivos cabecalhos ----

# fish_ID_tbl_fasta <- fish_ID_tbl %>% 
# fish_ID_tbl_fasta <- grouped_by_ID_BLASTid %>%
fish_ID_tbl_fasta <- raw_results_tbl %>%
  # filter(`Final ID (BLASTn)`%in% c("Poeciliidae sp.")) %>% #### corrigir
  # filter(`Final ID (BLASTn)`%in% c("Astyanax lacustris", "Astyanax paranae", "Astyanax bimaculatus")) %>% 
  # filter(`Final ID (BLASTn)`%in% c("Astyanax bimaculatus", "Piabina argentea")) %>% 
  # filter(`Final ID (BLASTn)`%in% c("Coptodon rendalli", "Coptodon zillii")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Eigenmannia virescens")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Eigenmannia virescens")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Gymnotus carapo", "Gymnotus sylvius")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Hemigrammus marginatus")) %>%
  filter(`Final ID (BLASTn)`%in% c("Hoplias intermedius", "Hoplias malabaricus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Iheringichthys labrosus", "Pimelodus fur", "Pimelodus mysteriosus", "Pimelodus albicans", "Pimelodus maculatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Megaleporinus garmani", "Leporinus reinhardt", "Leporinus crassilabris", "Leporinus steindachneri", "Leporinus octofasciatus", "Leporinus piau")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Megaleporinus garmani", "Leporinus reinhardt", "Leporinus crassilabris", "Leporinus steindachneri", "Leporinus octofasciatus", "Leporinus piau", "Leporellus vittatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Moenkhausia costae")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Myleus micans", "Piaractus brachypomus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Coptodon rendalli", "Coptodon zillii", "Oreochromis niloticus", "Oreochromis sp.")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Oreochromis niloticus", "Oreochromis sp.")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Iheringichthys labrosus","Pimelodus fur", "Pimelodus mysteriosus", "Pimelodus albicans", "Pimelodus maculatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Pygocentrus nattereri", "Serrasalmus brandtii")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Rhamdia quelen")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Salmo salar", "Brycon orthotaenia")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Trichomycterus sp.", "Astyanax lacustris")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Wertheimeria maculata")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Acestrorhynchus falcatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Aequidens pallidus", "Cichlasoma dimerus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Astyanax cf fasciatus", "Astyanax paranae", "Psalidodon rivularis")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Prochilodus harttii", "Prochilodus argenteus", " Prochilodus vimboides")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Astyanax bimaculatus", "Piabina argentea", "Planaltina myersi", "Knodus sp.")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Astyanax lacustris", "Astyanax bimaculatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Poecilia latipinna", "Poeciliidae sp.")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Eigenmannia virescens", "Eigenmannia sp.")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Pseudoplatystoma corruscans", "Pseudoplatystoma reticulatum")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Odontostilbe sp", "Serrapinnus heterodon")) %>%
  mutate(`Final ID (BLASTn)` = gsub(" ", "_", `Final ID (BLASTn)`)) %>% # Unindo os nomes com _
  # mutate(`Curated ID` = gsub(" ", "_", `Curated ID`)) %>% # Unindo os nomes com _
  unite(header, `ASV header`, `Final ID (BLASTn)`, sep = "|",`BLASTn pseudo-score`, `OTU`) %>%  # Criando uma nova coluna com o cabeçalho da ASV e o nome de spp atribuido # em 11/23 para obter a arvore de ASVs e definir as IDs
  # unite(header, "ASV header") %>%  # Criando uma nova coluna com o cabeçalho da ASV e o nome de spp atribuido
  reframe(header, `ASV (Sequence)`) %>% 
  unique() 

View(fish_ID_tbl_fasta)

## Criando o FASTA ----
{
  ASVs_all <- c(rbind(fish_ID_tbl_fasta$header, fish_ID_tbl_fasta$`ASV (Sequence)`)) # Unindo as informacoes em formato FASTA
  # write(ASVs_all, "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/alinhamentos/ASVs.fas") # Criando o arquivo
  # write(ASVs_all, "/home/gabriel/projetos/lagoa_ingleses/results/tree/2024/distancia_grupos_arvore/ASVs.fas") # Criando o arquivo  
  write(ASVs_all, "/home/gabriel/projetos/LI_paper/results/tree/distance/ASVs.fas") # Criando o arquivo  
}

# Lendo o FASTA
algn_in_ASVs_nofilt <- readDNAStringSet(filepath = "/home/gabriel/projetos/LI_paper/results/tree/distance/ASVs.fas")
# BrowseSeqs(algn_in_ASVs_nofilt)

{
  # Obtendo as seqs do LGC12sDB 
  
  fish_DB_jun23 <- readDNAStringSet(filepath = "/home/heron/prjcts/eDNA_fish/DB/jun23/DB/LGC12Sdb_261seqs-jun23-pretty_names_noGaps.fasta")

  # Filtrando as refs
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("JQ | Phalloceros sp | 1551", "JQ | Acinocheirodon melanogramma | 1547")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Oreochromis niloticus | 2266", "SF | Coptodon rendalli | 4861")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Eigenmannia virescens | 1118", "SF | Eigenmannia virescens | 1119")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Hemigrammus marginatus | 1382b", "SF | Hemigrammus marginatus | 1382")]
  fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("JQ | Hoplias malabaricus | 7822", "JQ | Hoplias sp | 6630", "SF | Hoplias malabaricus | 1346", "SF | Hoplias malabaricus | 1234", "JQ | Hoplias brasiliensis | 1772", "SF | Hoplias intermedius | 1377")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("Out | Leporellus vittatus | LC104399.1")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Myleus micans | 1162", "SF | Myleus micans | 0763", "SF | Colossoma macropomum | 2872", "SF | Colossoma macropomum | 2871", "JQ | Colossoma macropomum | 2871")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Coptodon rendalli | 4861", "SF | Oreochromis niloticus | 2266")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Pygocentrus piraya | 1229", "SF | Serrasalmus brandtii | 1230", "JQ | Serrasalmus brandtii | 6261", "JQ | Serrasalmus brandtii | 7890", "SF | Serrasalmus brandtii | 0872")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Rhamdia quelen | 1361")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Trichomycterus sp | 1373")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Franciscodoras marmoratus | 0708", "JQ | Wertheimeria maculata | 7817", "JQ | Wertheimeria maculata | 7906")]
  # fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Acestrorhynchus lacustris | 0917", "SF | Acestrorhynchus lacustris | 1136")]
  # fish_DB_filter <- fish_DB_jun23[str_detect(names(fish_DB_jun23), "Leporinus | Megaleporinus")]
  # fish_DB_filter <- fish_DB_jun23[str_detect(names(fish_DB_jun23), "Pimelodus")]
}

fish_DB_filter

# Unindo as refs com as ASVs

fish_ASVs_DB <- c(
                  fish_DB_filter
                  ,
                  algn_in_ASVs_nofilt
                  )

fish_ASVs_DB

# Corrigindo a orientacao 

fish_ASVs_DB_oriented <- fish_ASVs_DB %>% 
  OrientNucleotides() %>% 
  AlignSeqs() 
# %>% 
#   subseq(start = 4,
#          end = 170)
  
# Ver se e necessario trimmar
BrowseSeqs(fish_ASVs_DB_oriented)

# Matriz de distancias

dist_ASVs <-
  DistanceMatrix(fish_ASVs_DB_oriented,
               type = "matrix",
               includeTerminalGaps = FALSE,
               penalizeGapLetterMatches = FALSE,
               # penalizeGapMatches = FALSE,
               # correction = "none",
               correction = "Jukes-Cantor",
               processors = 1,
               verbose = TRUE)

View(dist_ASVs) 

# Salvar o fasta
Biostrings::writeXStringSet(x = fish_ASVs_DB_oriented,
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/distancia_grupos_arvore/Poeciliidae.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_hoplias.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_astyanax.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_moenkhausia.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_hemigrammus.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_rhamdia.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_leporinus.fas"
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_eigenmannia.fas"
)
                            
                            
                            
                            
                            
                            
                            
                            
                            