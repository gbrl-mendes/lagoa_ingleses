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

## Arvore ASVs ----

# Sequencias da tabela usada nos plots
asvs_tbl <- read.csv(file = "/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/arvore/ASVs_seqs.csv", sep = ";")

# Criando a coluna cabec unindo o nome das ASVs com o nome das especies
asvs_tbl <- asvs_tbl %>% 
  as_tibble() %>%
  select(inicio, especie, seq) %>% 
  mutate(cabec = paste0(inicio, "|", especie))

# Deixando apenas o nome do cabecalho e as sequencias e reordenando
asvs_tbl <- asvs_tbl[,c(4,3)]

# Criando funcao para transformar a tabela em FASTA
{writeFasta <- function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(data[rowNum,"cabec"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}}

# Criando o FASTA apartir da tabela
writeFasta(asvs_tbl, "/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/arvore/asvs.fasta")

# Input ASVs
asvs_fasta <- readDNAStringSet(
  # filepath = "/home/gabriel/projetos/peixes-eDNA/analises/pipeline-modelo/results/pipeline-modelo-all_ASVs_all_primers.fasta") %>%
  filepath = "/home/gabriel/projetos/lagoa_ingleses/tabelas/raw/arvore/asvs.fasta") %>%
  RemoveGaps()

BrowseSeqs(asvs_fasta)

# Alinhando as ASVs
asvs_align <- AlignSeqs(myXStringSet = asvs_fasta,
                        refinements = 100,
                        iterations = 100,
                        verbose = TRUE)
BrowseSeqs(asvs_align)

# Salvando o alinhamento
writeXStringSet(x = asvs_align,
                file = "~/projetos/lagoa_ingleses/results/tree/setembro/asvs_align.fasta")

# Calculo da matriz de distancia
asvs_dist <- DistanceMatrix(myXStringSet = asvs_align,
                            includeTerminalGaps = FALSE,
                            correction = "Jukes-Cantor",
                            processors = 20,
                            verbose = TRUE)

# Construcao do HeatMap

## Neste caso o uso do heatmap e' mais ludico. Quando foram utilizadas todas ASVs para construcao da arvore, 
## haviam sequencias tao distoantes que entre elas havia uma distancia infinita na matriz.
## O heatmap no momento entao serviu para rastrear quais eram essas sequencias e retira-las manualmente.
## E' bom manter esse trecho de codigo caso seja necessario no futuro.

pheatmap::pheatmap(asvs_dist,
                   fontsize = 2,
                   col = brewer.pal(n = 9,name = "YlOrRd"),
                   treeheight_col = 0,
                   filename = "/home/gabriel/projetos/lagoa_ingleses/results/tree/setembro/heatmap_dist.pdf", ## avaliar quais ASVs estao sendo outliers abrindo o pdf e olhando seus nomes
                   width = 8,
                   height = 8
)

# Construcao da arvore

asvs_tree <- ape::njs(asvs_dist)

# Plotando

print(asvs_tree_plot <- ggtree(tr = asvs_tree) +
        theme_tree2() +
        geom_tiplab(offset = 0,  
                    align = T) +
        xlim(0, 2)
)

# Salvando a figura em pdf

asvs_tree_plot %>% ggsave(filename = "/home/gabriel/projetos/lagoa_ingleses/results/tree/setembro/asvs_tree.pdf", 
                          device = "pdf", 
                          width = 24,
                          height = 24,
                          dpi = 600)

# Salvando a arvore em NWK
ape::write.tree(phy = asvs_tree,
                file = "~/projetos/lagoa_ingleses/tree/asvs_tree.nwk")

# Extraindo as IDs

asvs_ids <- get_taxa_name(asvs_tree_plot) %>% rev()


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
