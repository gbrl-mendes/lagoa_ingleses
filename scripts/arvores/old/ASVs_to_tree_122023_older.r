---
  title: "Extraindo ASVs da raw_results_tbl para um FASTA"
author: "Gabriel Mendes"
date: "12/2023"
---
  
  ## Atencao! Esse script funciona a partir da extracao das seqs de ASVs contidas no script 102023_tile-plots, 
  ## onde estao as analises realizadas com os dados dos sequenciamoentos 2, 4 e 5.
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
fish_ID_tbl_fasta <- grouped_by_ID_BLASTid %>%
  filter(`Class (BLASTn)` %in% "Actinopteri") %>% 
  mutate(`Final ID (BLASTn)` = gsub(" ", "_", `Final ID (BLASTn)`)) %>% # Unindo os nomes com _
  # mutate(`Curated ID` = gsub(" ", "_", `Curated ID`)) %>% # Unindo os nomes com _
  unite(header, `ASV header`, `Final ID (BLASTn)`, sep = "|",`BLASTn pseudo-score`) %>%  # Criando uma nova coluna com o cabeçalho da ASV e o nome de spp atribuido # em 11/23 para obter a arvore de ASVs e definir as IDs
  # unite(header, "ASV header") %>%  # Criando uma nova coluna com o cabeçalho da ASV e o nome de spp atribuido
  reframe(header, `ASV (Sequence)`) %>% 
  unique() 

View(fish_ID_tbl_fasta)

## Criando o FASTA ----

ASVs_all <- c(rbind(fish_ID_tbl_fasta$header, fish_ID_tbl_fasta$`ASV (Sequence)`)) # Unindo as informacoes em formato FASTA
write(ASVs_all, "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/ASVs_all.fas") # Criando o arquivo


# Criando objeto Biostring para fazer alinhamento 

algn_in_ASVs_nofilt <- readDNAStringSet(filepath = "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/ASVs_all.fas")
# BrowseSeqs(algn_in_ASVs_nofilt)

# Realizando alinhamento

algn_ASVs_out <- algn_in_ASVs_nofilt %>% 
  OrientNucleotides() %>% 
  AlignSeqs()

BrowseSeqs(algn_ASVs_out)

# Construindo um fasta para fazer no MEGAX

oriented_ASVs_out <- algn_in_ASVs_nofilt %>% 
  OrientNucleotides()

# Salvar as seqs orientadas
Biostrings::writeXStringSet(x = oriented_ASVs_out,
                            file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/ASVs_oriented.fas")


## Realizando alinhamento com o LGC12sDB ----

# Obtendo as seqs do LGC12sDB
fish_DB_jun23 <- readDNAStringSet(filepath = "/home/heron/prjcts/eDNA_fish/DB/jun23/DB/LGC12Sdb_261seqs-jun23-pretty_names_noGaps.fasta")

BrowseSeqs(fish_DB_jun23)

# Unindo o objeto das ASVs com o LGC12sDB   
# ASVs_12sDB <- c(algn_in_ASVs_001, fish_DB_jun23)
ASVs_12sDB <- c(algn_in_ASVs_nofilt, fish_DB_jun23)

# Alinhar as ASVs e as seqs do LGC12sDB
algn_ASVs_12sDB <- ASVs_12sDB %>% 
  OrientNucleotides() %>% # Resolvendo a orientacao das seqs
  AlignSeqs() 

BrowseSeqs(algn_ASVs_12sDB)

# Filtrando grupos externos do LGC12sDB
fish_DB_jun23_fil <- fish_DB_jun23
# [!(names(fish_DB_jun23) %>% str_detect(pattern = "Out"))] ## removendo as ASVs com "Out" no header

# Verificar quais sequências foram filtradas
filtrado_db <- names(fish_DB_jun23[!(names(fish_DB_jun23) %in% names(fish_DB_jun23_fil))]) # os nomes foram armazenados no objeto filtrado
list(filtrado_db)

# BrowseSeqs(fish_DB_jun23)
# BrowseSeqs(fish_DB_jun23_fil)

# Agora, devo alinhar as seqs do LGC12sDB sozinhas, sem as ASVs, e em seguida trimar elas
# Após esse alinhamento irei novamente realizar o alinhamento das ASVs.
# Essas etapas foram realizadas para que as ASVs alinhem-se apenas com a regiao homologa das seqs do LGC12sDB.

# Trimmando o LGC12sDB
trim_fish_DB_jun23 <- fish_DB_jun23_fil %>% 
  OrientNucleotides() %>% # Resolvendo a orientacao das seqs
  AlignSeqs() %>% 
  subseq(start = 34,
         end = 246)

# BrowseSeqs(trim_fish_DB_jun23) # Observando resultado do alinhamento do DB trimmado

# Salvar o alinhamento
Biostrings::writeXStringSet(x = trim_fish_DB_jun23,
                            file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/trim_fish_DB_jun23.fas")

# Alinhando as ASVs com o LGC12sDB trimmado
trim_algn_ASVs_12sDB <- c(algn_in_ASVs_nofilt,trim_fish_DB_jun23) %>% 
  DECIPHER::RemoveGaps() %>% 
  OrientNucleotides(reference = 2) %>% 
  OrientNucleotides() %>% 
  AlignSeqs()

BrowseSeqs(trim_algn_ASVs_12sDB)

# Salvar o alinhamento
Biostrings::writeXStringSet(x = trim_algn_ASVs_12sDB,
                            file =  "/home/gabriel/projetos/lagoa_ingleses/results/tree/2023/trim_algn_ASVs_12sDB_nofilt")

# Calcular a matriz de distancia ----
dna_dist <- DECIPHER::DistanceMatrix(myXStringSet = trim_algn_ASVs_12sDB,
                                     includeTerminalGaps = FALSE,      #test both versions on the final output
                                     correction = "Jukes-Cantor",
                                     processors = 20,
                                     verbose = TRUE)


# Construir arvore ----

dna_tree <- ape::njs(dna_dist)

# Salvar arvore ----
# ape::write.tree(phy = dna_tree, file = "~/projetos/lagoa_ingleses/results/tree/2023/asv-LI_db_tree_nofilt.nwk")

# comentario 13_11_2023
# criacao da arvore com os outgroups, sem excluir por tamanho de ASV e tirando as Curated IDs dos nomes

ape::write.tree(phy = dna_tree, file = "~/projetos/lagoa_ingleses/results/tree/2023/20_12_ASV_all.nwk")

# comentario 27_10_2023:
# deu um erro ao executar o comando ape::write.tree()
# e por isso tive que comentar o /correction = "Jukes-Cantor"/
# no codigo que constroi a matriz de distancia DECIPHER::DistanceMatrix()


# Customize (manually)-----

## on iTOL/MEGA software, drag and drop the .nwk file. Use tools to generate a pretty visualization

# Verificando as identificacoes das Seqs ----

# Apos avaliar as relacoes que a arvore formou, e` necessario realizar o alinhamento das sequencias
# cujo agrupamento da arvore formou grupos que nao pertencem ao mesmo clado


# Avaliando Pseudoplastystoma e Pimelodus
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in%
                       c("SF | Pseudoplatystoma corruscans | 1135",
                         "ASV_280_175bp|Pseudoplatystoma_corruscans"
                         # "ASV_385_173bp|Pimelodus_spp.",
                         # "ASV_329_174bp|Pimelodus_pohli",
                         # "ASV_330_174bp|Pimelodus_pohli",
                         # "ASV_386_173bp|Pimelodus_maculatus",
                         # "ASV_387_173bp|Pimelodus_maculatus",
                         # "ASV_388_173bp|Pimelodus_fur",
                         # "ASV_389_173bp|Pimelodus_fur"
                         )] %>% 
  BrowseSeqs()

# Avaliando P. costatus e P. hartii
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(#"ASV_528_170bp|Salmo_salar",
                         "ASV_501_171bp|Prochilodus_costatus",
                         "ASV_500_171bp|Prochilodus_costatus",
                         "ASV_488_171bp|Prochilodus_costatus"
                         ,
                         "JQ | Prochilodus hartii | 1765",
                         "JQ | Prochilodus hartii | 7033"
                         ,
                         "JQ | Prochilodus costatus | 2860",
                         "SF | Prochilodus costatus | 0664",
                         "SF | Prochilodus costatus | 1339",
                         "SF | Prochilodus costatus | 1281"
                         ,
                         "SF | Prochilodus argenteus | 0742",
                         "JQ | Prochilodus argenteus | 5612",
                         "SF | Prochilodus argenteus | 0739",
                         "SF | Prochilodus argenteus | 0779"
                         )] %>% 
  BrowseSeqs()

# Avaliando L. taeniatus
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c("ASV_436_172bp|Leporinus_taeniatus", 
                         "ASV_437_172bp|Leporinus_taeniatus", 
                         "SF | Leporinus taeniatus | 0985", 
                         "JQ | Leporinus taeniatus | 1553", 
                         "SF | Leporinus taeniatus | 1317")] %>% 
  BrowseSeqs()

# Avaliando M. elongatus e M. reinhardti
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(
                         # "SF | Megaleporinus elongatus | 1032",
                         # "SF | Megaleporinus elongatus | 1033",
                         # "SF | Megaleporinus elongatus | 0772",
                         # "SF | Megaleporinus obtusidens | 0758", 
                         # "SF | Megaleporinus obtusidens | 0759", 
                         # "SF | Megaleporinus obtusidens | 0757", 
                         "JQ | Megaleporinus elongatus | 1762",
                         "JQ | Megaleporinus elongatus | 1761",
                         "JQ | Megaleporinus elongatus | 6633",
                         # "SF | Megaleporinus reinhardti | 1344",
                         # "JQ | Megaleporinus garmani | 4290",
                         # "JQ | Megaleporinus garmani | 5636"
                         # ,
                         # "ASV_540_170bp|Leporinus_reinhardti", 
                         # "ASV_539_170bp|Leporinus_reinhardti"
                         # , 
                         "ASV_153_170bp|Megaleporinus_elongatus",
                         "ASV_147_170bp|Megaleporinus_elongatus"
                     )] %>% 
  BrowseSeqs()

# Avaliando H. intermedius e H. spp
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(
                         "ASV_672_167bp|Hoplias_spp.",
                         "ASV_673_167bp|Hoplias_spp.",
                         "SF | Hoplias intermedius | 1377",
                         "JQ | Hoplias brasiliensis | 1772",
                         "ASV_33_202bp|Hoplias_intermedius",
                         "ASV_80_194bp|Hoplias_intermedius",
                         "ASV_686_167bp|Hoplias_intermedius",
                         "ASV_684_167bp|Hoplias_intermedius",
                         "ASV_21_167bp|Hoplias_intermedius",
                         "ASV_29_167bp|Hoplias_intermedius",
                         "ASV_756_167bp|Hoplias_intermedius",
                         "ASV_757_167bp|Hoplias_intermedius",
                         "ASV_14_167bp|Hoplias_intermedius",
                         "ASV_74_169bp|Hoplias_intermedius",
                         "ASV_69_194bp|Hoplias_intermedius",
                         "ASV_705_167bp|Hoplias_intermedius",
                         "ASV_703_167bp|Hoplias_intermedius",
                         "ASV_22_167bp|Hoplias_intermedius",
                         "ASV_31_167bp|Hoplias_intermedius",
                         "ASV_27_167bp|Hoplias_intermedius",
                         "ASV_687_167bp|Hoplias_intermedius"
                         ,
                         "JQ | Hoplias sp | 6630"
                         # ,
                         # "JQ | Hoplias malabaricus | 7822",
                         # "SF | Hoplias malabaricus | 1234",
                         # "SF | Hoplias malabaricus | 1346",
                         # "ASV_15_167bp|Hoplias_malabaricus"
                     )] %>% 
  BrowseSeqs()

# Avaliando H. intermedius e H. spp
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(
                       "SF | Orthospinus franciscensis | 1072", 
                         # "SF | Orthospinus franciscensis | 1071c",
                         "ASV_359_175bp|Moenkhausia_costae",
                         # "ASV_337_175bp|Moenkhausia_costae",
                         # "ASV_346_175bp|Moenkhausia_costae",
                         # "ASV_424_175bp|Moenkhausia_costae",
                         # "ASV_384_175bp|Moenkhausia_costae",
                         # "ASV_352_175bp|Moenkhausia_costae",
                         # "ASV_364_175bp|Moenkhausia_costae",
                         "ASV_287_175bp|Orthospinus_franciscensis"
                         ,
                         # "ASV_326_175bp|Moenkhausia_costae",
                         # "ASV_37_168bp|Moenkhausia_costae",
                         # "ASV_376_175bp|Moenkhausia_costae",
                         # "ASV_473_197bp|Moenkhausia_costae",
                         # "ASV_362_175bp|Moenkhausia_costae",
                         # "ASV_367_175bp|Moenkhausia_costae",
                         # "ASV_343_175bp|Moenkhausia_costae",
                         # "ASV_394_175bp|Moenkhausia_costae",
                         # "ASV_368_175bp|Moenkhausia_costae",
                         # "ASV_356_175bp|Moenkhausia_costae",
                         # "ASV_345_175bp|Moenkhausia_costae",
                         "JQ | Moenkhausia costae | 7812",
                         "SF | Moenkhausia costae | 1074"
                         # ,
                         # "ASV_497_207bp|Moenkhausia_costae"
                         )] %>% 
  BrowseSeqs()

# Avaliando A. spp e A. lacustris
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(
                         "ASV_626_168bp|Astyanax_lacustris", 
                         "ASV_627_168bp|Astyanax_lacustris",
                         "ASV_37_168bp|Astyanax_lacustris",
                         "ASV_34_168bp|Astyanax_spp.",
                         "ASV_643_168bp|Astyanax_lacustris",
                         "ASV_38_168bp|Astyanax_spp.",
                         "ASV_651_168bp|Astyanax_lacustris",
                         "ASV_205_172bp|Astyanax_lacustris",
                         "ASV_652_168bp|Astyanax_lacustris",
                         "JQ | Astyanax lacustris | 7880", 
                         "SF | Astyanax lacustris | 1239", 
                         "JQ | Astyanax cf lacustris | 5672"
                         )] %>% 
  BrowseSeqs()

# Avaliando A. spp e A. fasciatus
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(
                         "ASV_44_168bp|Astyanax_fasciatus", 
                         "ASV_57_168bp|Astyanax_fasciatus",
                         "ASV_657_168bp|Astyanax_spp.",
                         "ASV_658_168bp|Astyanax_spp.",
                         "ASV_35_168bp|Astyanax_fasciatus"
                         )] %>% 
  BrowseSeqs()

# Avaliando Poecilia reticulata e Acinocheirodon melanogramma
trim_algn_ASVs_12sDB[names(trim_algn_ASVs_12sDB) %in% 
                       c(
                         "ASV_530_170bp|Poecilia_reticulata",
                         "ASV_529_170bp|Acinocheirodon_melanogramma",
                         "ASV_156_170bp|Acinocheirodon_melanogramma",
                         "ASV_424_175bp|Acinocheirodon_melanogramma",
                         "JQ | Acinocheirodon melanogramma | 1547"
                         )] %>% 
  BrowseSeqs()

















# Criando uma lista com as spp do DB, para realizar a curadoria no Google Sheets

db_spp <- fish_DB_jun23 %>% names()

db_spp_tbl <-
  db_spp %>%
  str_split("\\|", simplify = TRUE) %>%
  as.data.frame() %>%
  transmute(Bacia = trimws(V1), Species = trimws(V2))

View(db_spp_tbl)
