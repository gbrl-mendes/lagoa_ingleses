---
  title: "Extraindo ASVs da raw_results_tbl para um FASTA"
author: "Gabriel Mendes"
date: "12/2023"
---
  
  ## Atencao! Esse script funciona a partir da extracao das seqs de ASVs contidas no script 012024_post_pipe.r, 
  ## onde estao as analises realizadas com os dados dos sequenciamentos 2, 4, 5 e EM156.
  
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
fish_ID_tbl_fasta <- pre_grouped_by_ID_BLASTid %>%
  filter(`Final ID (BLASTn)`%in% c("Coptodon rendalli", "Oreochromis niloticus", "Oreochromis sp.", "Coptodon zillii")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Hoplias intermedius", "Hoplias malabaricus")) %>% 
  # filter(`Final ID (BLASTn)`%in% c("Psalidodon rivularis", "Astyanax lacustris", "Astyanax paranae", "Astyanax bimaculatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Moenkhausia costae")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Hemigrammus marginatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Rhamdia quelen")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Leporinus steindachneri", "Leporinus crassilabris", "Leporinus octofasciatus", "Leporinus reinhardti", "Leporinus piau", "Megaleporinus garmani", "Leporellus vittatus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Pimelodus albicans","Pimelodus fur", "Pimelodus maculatus", "Pimelodus mysteriosus")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Piabina argentea")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Gymnotus carapo", "Gymnotus sylvius")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Eigenmannia virescens")) %>%
  # filter(`Final ID (BLASTn)`%in% c("Myleus micans")) %>%
  mutate(`Final ID (BLASTn)` = gsub(" ", "_", `Final ID (BLASTn)`)) %>% # Unindo os nomes com _
  # mutate(`Curated ID` = gsub(" ", "_", `Curated ID`)) %>% # Unindo os nomes com _
  unite(header, `ASV header`, `Final ID (BLASTn)`, sep = "|",`BLASTn pseudo-score`, `OTU`) %>%  # Criando uma nova coluna com o cabeçalho da ASV e o nome de spp atribuido # em 11/23 para obter a arvore de ASVs e definir as IDs
  # unite(header, "ASV header") %>%  # Criando uma nova coluna com o cabeçalho da ASV e o nome de spp atribuido
  reframe(header, `ASV (Sequence)`) %>% 
  unique() 

View(fish_ID_tbl_fasta)

## Criando o FASTA ----

ASVs_all <- c(rbind(fish_ID_tbl_fasta$header, fish_ID_tbl_fasta$`ASV (Sequence)`)) # Unindo as informacoes em formato FASTA
# write(ASVs_all, "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/alinhamentos/ASVs.fas") # Criando o arquivo
write(ASVs_all, "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/alinhamentos/ASVs.fas") # Criando o arquivo

# Lendo o FASTA
algn_in_ASVs_nofilt <- readDNAStringSet(filepath = "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/alinhamentos/ASVs.fas")
# BrowseSeqs(algn_in_ASVs_nofilt)

# Obtendo as seqs do LGC12sDB

fish_DB_jun23 <- readDNAStringSet(filepath = "/home/heron/prjcts/eDNA_fish/DB/jun23/DB/LGC12Sdb_261seqs-jun23-pretty_names_noGaps.fasta")

# Filtrando as refs
fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Oreochromis niloticus | 2266", "SF | Coptodon rendalli | 4861")]
# fish_DB_filter <- fish_DB_jun23[names(fish_DB_jun23) %in% c("SF | Myleus micans | 0763", "SF | Myleus micans | 1162")]

# Unindo as refs com as ASVs

fish_ASVs_DB <- c(fish_DB_filter, algn_in_ASVs_nofilt)

# Corrigindo a orientacao e alinhando

fish_ASVs_DB_oriented <- fish_ASVs_DB %>% 
  OrientNucleotides() %>% 
  AlignSeqs()

BrowseSeqs(fish_ASVs_DB_oriented)

# Salvar o alinhamento
Biostrings::writeXStringSet(x = fish_ASVs_DB_oriented,
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/alinh_coptodon.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_hoplias.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_astyanax.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_moenkhausia.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_hemigrammus.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_rhamdia.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_leporinus.fas")
                            # file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/analises_daniel/ASVS_all_ori_eigenmannia.fas")








