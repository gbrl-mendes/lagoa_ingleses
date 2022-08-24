---
{library(Biostrings)
  library(DECIPHER)
  library(dplyr)
  library(stringr)
  library(ggtree)
  library(RColorBrewer)
}

## Lendo Fasta DB
db_fasta <- readDNAStringSet(
  filepath = "/home/heron/prjcts/fish_eDNA/DB/mai22/DB/LGC12Sdb_251seqs-mai22-pretty_names_noGaps.fasta") %>% 
  RemoveGaps()

## Lendo Fasta ASVs
asvs_fasta <- readDNAStringSet(
  filepath = "/home/gabriel/projetos/peixes-eDNA/analises/pipeline-modelo/results/pipeline-modelo-all_ASVs_all_primers.fasta") %>%
  RemoveGaps()

BrowseSeqs(db_fasta)
BrowseSeqs(asvs_fasta)

## Removendo os outgrups de db_fasta

db_fasta_filt <- db_fasta[!names(db_fasta) %>% str_detect(pattern = "Out")]

## Conferindo se foram retirados
filt <- names(db_fasta[!names(db_fasta) %in% names(db_fasta_filt)])
filt

## Combinar as seqs

db_asvs_fasta <- c(db_fasta_filt, asvs_fasta) %>% 
  RemoveGaps() %>%
  unique()

## Alinhar as seqs do banco e ASVs

db_asvs_align <- AlignSeqs(myXStringSet = db_asvs_fasta,
                           refinements = 100,
                           iterations = 100,
                           verbose = TRUE)

BrowseSeqs(db_asvs_align)












