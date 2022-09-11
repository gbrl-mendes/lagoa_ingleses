---

## Construcao de Arvores
  
## Pacotes  
{
  library(Biostrings)
  library(DECIPHER)
  library(dplyr)
  library(stringr)
  library(ggtree)
  library(RColorBrewer)
}
----

## Arvore ASVs ----

# Input ASVs
asvs_fasta <- readDNAStringSet(
  filepath = "/home/gabriel/projetos/peixes-eDNA/analises/pipeline-modelo/results/pipeline-modelo-all_ASVs_all_primers.fasta") %>%
  RemoveGaps()

BrowseSeqs(asvs_fasta)

# Alinhando as ASVs
asvs_align <- AlignSeqs(myXStringSet = asvs_fasta,
                        refinements = 100,
                        iterations = 100,
                        verbose = TRUE)
BrowseSeqs(asvs_align)

# Filtrando ASVs muito diferentes

asvs_to_remove <- c("ASV_11_195bp_merged|Cichlidae|NA|NA|Cichlidae",
                    "ASV_204_85bp_merged|NA|NA|NA|NA",
                    "ASV_201_93bp_merged|NA|NA|NA|NA",
                    "ASV_209_81bp_merged|NA|NA|NA|NA",
                    "ASV_208_81bp_merged|NA|NA|NA|NA",
                    "ASV_210_76bp_merged|NA|NA|NA|NA",
                    "ASV_203_86bp_merged|NA|NA|NA|NA",
                    "ASV_200_122bp_merged|NA|NA|NA|NA",
                    "ASV_202_92bp_merged|NA|NA|NA|NA",
                    "ASV_193_153bp_merged|Characidae|NA|NA|Characidae",
                    "ASV_191_159bp_merged|NA|NA|NA|NA",
                    "ASV_188_164bp_merged|NA|NA|NA|NA",
                    "ASV_189_164bp_merged|NA|NA|NA|NA",
                    "ASV_24_191bp_merged|NA|NA|NA|NA",
                    "ASV_31_178bp_merged|NA|NA|NA|NA",
                    "ASV_197_138bp_merged|NA|NA|NA|NA",
                    "ASV_13_194bp_merged|Characidae|Moenkhausia|Moenkhausia costae|Moenkhausia costae",
                    "ASV_14_194bp_merged|Erythrinidae|Hoplias|Hoplias malabaricus|Hoplias malabaricus/sp",
                    "ASV_15_194bp_merged|Erythrinidae|Hoplias|NA|Hoplias brasiliensis/intermedius",
                    "ASV_10_195bp_merged|Pimelodidae|Rhamdia|Rhamdia quelen|Rhamdia quelen",
                    "ASV_23_192bp_merged|Cichlidae|Oreochromis|Oreochromis niloticus|NA niloticus/rendalli",
                    "ASV_4_201bp_merged|NA|NA|NA|NA",
                    "ASV_16_193bp_merged|NA|NA|NA|NA",
                    "ASV_88_169bp_merged|NA|NA|NA|NA",
                    "ASV_211_65bp_merged|NA|NA|NA|NA",
                    "ASV_198_135bp_merged|NA|NA|NA|NA",
                    "ASV_199_130bp_merged|NA|NA|NA|NA",
                    "ASV_195_145bp_merged|NA|NA|NA|NA",
                    "ASV_17_193bp_merged|NA|NA|NA|NA",
                    "ASV_27_182bp_merged|NA|NA|NA|NA",
                    "ASV_207_83bp_merged|Loricariidae|Harttia|NA|Harttia",
                    "ASV_206_83bp_merged|NA|NA|NA|NA",
                    "ASV_194_146bp_merged|NA|NA|NA|NA",
                    "ASV_205_83bp_merged|NA|NA|NA|NA",
                    "ASV_1_211bp_merged|NA|NA|NA|NA",
                    "ASV_2_208bp_merged|NA|NA|NA|NA",
                    "ASV_79_170bp_merged|NA|NA|NA|Hydrochaeris hydrochaeris isol",
                    "ASV_95_169bp_merged|Hominidae|Homo|Homo sapiens|Homo sapiens",
                    "ASV_70_171bp_merged|Hominidae|Homo|Homo sapiens|Homo sapiens",
                    "ASV_87_169bp_merged|Canidae|Canis|Canis familiaris|Canis familiaris",
                    "ASV_89_169bp_merged|Canidae|Canis|Canis familiaris|Canis familiaris",
                    "ASV_71_171bp_merged|NA|NA|NA|Oryctolagus cuniculus mitochon",
                    "ASV_80_170bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_78_170bp_merged|NA|NA|NA|Sus scrofa isolate Europe hapl",
                    "ASV_106_169bp_merged|NA|NA|NA|Sus scrofa breed Andaman Desi ",
                    "ASV_76_170bp_merged|NA|NA|NA|Didelphis albiventris mitochon",
                    "ASV_182_167bp_merged|NA|NA|NA|Cavia magna 12S ribosomal RNA ",
                    "ASV_61_172bp_merged|Hominidae|Homo|Homo sapiens|Homo sapiens",
                    "ASV_73_171bp_merged|Hominidae|Homo|Homo sapiens|Homo sapiens",
                    "ASV_72_171bp_merged|Hominidae|Homo|Homo sapiens|Homo sapiens",
                    "ASV_190_160bp_merged|Hominidae|Homo|Homo sapiens|Homo sapiens",
                    "ASV_99_169bp_merged|NA|NA|NA|NA",
                    "ASV_96_169bp_merged|NA|NA|NA|Hydrochaeris hydrochaeris isol",
                    "ASV_105_169bp_merged|NA|NA|NA|Hydrochaeris hydrochaeris isol",
                    "ASV_108_169bp_merged|NA|NA|NA|Hydrochaeris hydrochaeris isol",
                    "ASV_101_169bp_merged|NA|NA|NA|NA",
                    "ASV_93_169bp_merged|NA|NA|NA|NA",
                    "ASV_91_169bp_merged|NA|NA|NA|NA",
                    "ASV_92_169bp_merged|NA|NA|NA|NA",
                    "ASV_122_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_111_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_112_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_123_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_120_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_132_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_119_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_116_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_118_168bp_merged|Bovidae|Bos|Bos taurus|Bos taurus",
                    "ASV_192_158bp_merged|NA|NA|NA|Human DNA sequence from clone ",
                    "ASV_28_181bp_merged|NA|NA|NA|Progne chalybea mitochondrion,",
                    "ASV_29_180bp_merged|NA|NA|NA|Nannopterum brasilianus mitoch")

asvs_fil <- asvs_fasta[!(names(asvs_fasta) %in% asvs_to_remove)] #removendo as ASVs que estão na lista asvs_to_remove

# Verificar quais sequências foram filtradas
names(asvs_fasta[!(names(asvs_fasta) %in% names(asvs_fil))])

# Verificar as as ASVs que sobraram mas não estão na lista de seqs para filtrar
filtrado[!filtrado %in% asvs_to_remove] #este resultado possui apenas as sequências externas, ja que elas foram retiradas com o str_detect

# Verificar quais nomes estão na lista de seqs para filtrar mas não foram filtradas
asvs_to_remove[!asvs_to_remove %in% filtrado]

BrowseSeqs(dna_fasta_or)
BrowseSeqs(dna_fasta_fil)

# Arvore ASVs + LGC12sDB ----  
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

## Salvar o alinhamento

writeXStringSet(x = db_asvs_align,
                file = "/home/gabriel/projetos/lagoa_ingleses/results/tree/db_asvs_align.fasta")

## Calcular a matriz de distancia
dna_dist <- DistanceMatrix(myXStringSet = db_asvs_align,
                                     # includeTerminalGaps = TRUE,
                                     includeTerminalGaps = FALSE,      #test both versions on the final output
                                     correction = "Jukes-Cantor",
                                     processors = 20,
                                     verbose = TRUE)
