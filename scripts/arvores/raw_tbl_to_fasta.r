---
  title: "Extraindo ASVs da raw_results_tbl para um FASTA"
author: "Gabriel Mendes"
date: "30/06/2023"
---
  
  ## Atencao! Esse script funciona a partir da extracao das seqs de ASVs
  ## contidas no script 062023_tile-plots, onde estao as analises realizadas
  ## com os dados dos sequenciamoentos 2, 4 e 5.
  ## Para usar com outros dados sera necessario adaptar.

  ## Carregando bibliotecas ----
{
  library(Biostrings)
  library(phyloseq)
  library(ShortRead)
  library(stringr)
  library(tidyverse)
}

  ## Filtrando as ASVs com RRA < 0.01

fish_ID_tbl <- fish_ID_tbl %>% filter(RRA >= 0.01) %>% 
  mutate(`Curated ID` = gsub(" ", "_", `Curated ID`))

# Crie a função para formatar as linhas em formato FASTA
{
  formatar_linha_fasta <- function(`Curated ID`, `ASV header`, `ASV (Sequence)`) {
    cabecalho <- paste(`ASV header`, `Curated ID`, sep = "|")
    sequencia_fasta <- str_wrap(`ASV (Sequence)`, width = 80)
    return(paste(cat(cabecalho, "\n", sequencia_fasta, "\n")))
}

}

# Aplique a função a cada linha da tabela e obtenha as sequências formatadas
asvs_to_fasta <- mapply(formatar_linha_fasta, fish_ID_tbl$`Curated ID`, fish_ID_tbl$`ASV header`, fish_ID_tbl$`ASV (Sequence)`)

sink(file = "/home/gabriel/projetos/lagoa_ingleses/results/tree/ASVs_001.fas")
summary(mapply(formatar_linha_fasta, fish_ID_tbl$`Curated ID`, fish_ID_tbl$`ASV header`, fish_ID_tbl$`ASV (Sequence)`))
sink(file = NULL)
# Salve as sequências em um arquivo FASTA
writeLines(as.character(asvs_to_fasta), "/home/gabriel/projetos/lagoa_ingleses/results/tree/ASVs_001.fas")
