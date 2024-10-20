---
title: "Alpha diversity and NMDS plots Lagoa dos Ingleses "
author: "Gabriel Mendes"
date: "06/2023"
---
  
# Carregando bibliotecas ----
{
  library(Biostrings)
  library(DECIPHER)
  library(factoextra)
  library(future)
  library(ggh4x)
  library(ggord)
  library(ggpubr)
  library(ggvegan)
  library(Matrix)
  library(phyloseq)
  library(ShortRead)
  library(stringr)
  library(vegan)
  library(tidyverse)
}

# Caminhos ----
{
  prjct_path <- "~/projetos/lagoa_ingleses/"
  results_path <- paste0(prjct_path,"/results")
  figs_path <- paste0(results_path,"/figuras")
  tbl_path <- paste0(prjct_path,"/tabelas/raw/run_2_4_5")
  prjct_radical <- "eDNA_Lagoa-dos-Ingleses"
}

# Obtencao de dados ----
{
  raw_results_tbl <- read.csv(paste0(tbl_path,"/","run_2_4_5_lagoa_ingleses_v062023.csv"), sep = ",", check.names = FALSE) %>% tibble()
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Oreochromis niloticus")] <- "Coptodon sp."#"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Coptodon zillii")] <- "Coptodon sp."#"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Tilapia rendalli")] <- "Coptodon sp." #"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
  raw_results_tbl$`Curated ID`[raw_results_tbl$`Curated ID` %in% c("Cichla spp.")] <- "Coptodon sp." #"Tilapia rendalli/Coptodon zillii" # Oreochromis niloticus é Tilapia
}

## Tabelas ----

# Agrupar por ID e Blast ID
{
  grouped_by_ID_BLASTid <- raw_results_tbl %>%
    select(c(
      "Sample",
      "Expedition",
      "Year",
      "Sample.Name",
      "File_name",
      "Curated ID",
      "Relative abundance on sample",
      "Point",
      "New_name",
      "Class",
      "Order",
      "Abundance"
    )) %>% 
    group_by(Sample, `Curated ID`, Expedition, Point, Sample.Name, File_name, ) %>% 
    summarize(`Sample` = unique(Sample),
              `Curated ID` = unique(`Curated ID`),
              `Class`,
              `Order`,
              # `Class` = unique(`Class`),
              # `Order` = unique(`Order`),
              `Expedition` = unique(Expedition),
              `Year` = unique(Year),
              `Point` = unique(Point),
              `New_name` = unique(New_name),
              `Sample.Name` = unique(Sample.Name),
              `File_name` = unique(File_name),
              `Abundance`,
              `RRA` = sum(`Relative abundance on sample`)) %>%
    ungroup()
}


# Organizar ordem
{
  fish_ordens <- c("Characiformes", # ordem das ordens de peixes 
                   "Siluriformes", 
                   "Cichliformes", 
                   "Gymnotiformes", 
                   "Cyprinodontiformes", 
                   "Salmoniformes")
}  

{
  nfish_classes <- c("Aves", # ordem das classes de nao peixes
                     "Actinobacteria",
                     "Mammalia")
}

{
  nfish_ordens <- c("Actinomycetales", # ordem das ordens de nao peixes 
                    "Artiodactyla", 
                    "Carnivora",
                    "Didelphimorphia",
                    "Lagomorpha",
                    "Primates",
                    "Rodentia",
                    "Suliformes")
}

# Definir quais serao as especies e ordenar

{
  spp_fish <- raw_results_tbl %>%  # especies de peixes
    filter(Class == "Actinopteri") %>%
    mutate(words_count = str_count(`Curated ID`, "\\w+")) %>%
    filter(words_count >= 2) %>%
    select(`Curated ID`) %>%
    unique() %>% 
    arrange(`Curated ID`) %>%
    pull()
}

# tabela so peixes
{
  fish_ID_tbl <- grouped_by_ID_BLASTid %>% # ids apenas os peixes a nivel de spp
    mutate(`Curated ID` = factor(`Curated ID`,
                                 levels = rev(spp_fish))) %>% 
    mutate(`Order` = factor(`Order`,
                            levels = fish_ordens)) %>%
    mutate('Nivel' = ifelse(`Year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                            "Cheio", "Vazio")) %>%
    filter(`Curated ID` != "") %>%
    filter(!is.na(`Curated ID`)) %>% 
    arrange(`Curated ID`)
}


## Alpha diversity ----

# Longer table
{
  longer_tbl_IDS <- 
    fish_ID_tbl %>%
    filter(RRA >= 0.01) %>% # definicao de qual threshold de abundancia sera usado
    select(c("Sample", "Curated ID", "Point", "Year", "Nível", "New_name", "RRA")) %>%
    unique() %>% 
    mutate("New_name" = as.factor(New_name)) %>%
    filter(Sample %in% c("LI1-neo-mi", "L2_nov20", "L2_dez20", "L2_nov21", "L3_nov21", "L4_nov21")) %>% #com apenas as amostras que Daniel pediu
    rename("Nivel" = "Nível") %>%
    mutate("Nivel" = as.factor(Nivel))
}
    

## Alpha (within sample) diversity ##

# Common alpha diversity statistics include:
# Shannon: How difficult it is to predict the identity of a randomly chosen individual.
# Simpson: The probability that two randomly chosen individuals are the same species.
# Inverse Simpson: This is a bit confusing to think about. Assuming a theoretically 
# community where all species were equally abundant, this would be the number of 
# species needed to have the same Simpson index value for the community being analyzed.


# Fonte: https://grunwaldlab.github.io/analysis_of_microbiome_community_data_in_r/07--diversity_stats.html

## Calculo usando o tidyverse
{
  # Funcoes para calculo de riqueza, indice de shannon e simpson
  {
    # Funcao para calcular a riqueza de spp
    richness <- function(x){
      
      sum(x > 0)
    }
    
    # Funcao para calcular o indice de shannon
    shannon <- function(x){
      
      rabund <- x[x>0] /sum(x)
      -sum(rabund * log(rabund))
    }
    
    # Funcao para calcular o indice de simpson
    simpson <- function(x){
      
      n <- sum(x)
      
      sum(x*(x-1)/(n*(n-1)))  
    }
    
  }
  
# PLots from Riffomonas Project  
  longer_tbl_IDS %>% 
    group_by(Sample) %>% 
    summarise(sobs = richness(RRA),
              shannon = shannon(RRA)
              ,
              simpson = simpson(RRA),
              invsimpson = 1/simpson,
              n = sum(RRA)) %>% 
    pivot_longer(cols = c(sobs, shannon, invsimpson, simpson), 
                 names_to = "metric") %>% 
    ggplot(aes(x = n,
               y = value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow = 4, scales = "free_y")
  }

## Usando o Vegan
{
  longer_tbl_IDS %>% 
    group_by(Sample) %>% 
    summarise(sobs = specnumber(RRA),
              shannon = diversity(RRA, index = "shannon"),
              simpson = diversity(RRA, index = "simpson"),
              invsimpson = 1/simpson,
              n = sum(RRA)) %>% 
    pivot_longer(cols = c(sobs, shannon, invsimpson, simpson), 
                 names_to = "metric") %>% 
    ggplot(aes(x = n,
               y = value)) +
    geom_point() +
    geom_smooth() +
    facet_wrap(~metric, nrow = 4, scales = "free_y")
}

## Outros plots
{
  index_tbl <- longer_tbl_IDS %>% 
    group_by(Sample) %>% 
    summarise(sobs = richness(RRA),
              shannon = shannon(RRA),
              simpson = simpson(RRA),
              invsimpson = 1/simpson,
              n = sum(RRA))
  
  index_merge_tbl <- longer_tbl_IDS %>% 
    select(Sample, Point, Year, Nivel, New_name) %>% 
    unique() %>% 
    merge(index_tbl, by = "Sample") %>% 
    mutate(Sample = factor(Sample, levels = c("LI1-neo-mi",
                                              "L2_nov20",
                                              "L2_dez20",
                                              "L1_out21",
                                              "L2_out21",
                                              "L3_out21",
                                              "L4_out21",
                                              "L1_nov21",
                                              "L2_nov21",
                                              "L3_nov21",
                                              "L4_nov21")))
  ## Scatterplot
  index_merge_tbl %>% 
    ggplot(aes(x = Sample,
               y = sobs,
               color = Nivel)) +
    geom_point() +
    theme_minimal()
  
  ## Boxplot
  # boxplot_alfa <-
    index_merge_tbl %>% 
    ggplot(aes(x = Nivel,
               y = sobs,
               fill = Nivel)) + 
    geom_boxplot(
      # fill = "blue",
      color = "black",
      size = 0.3) + 
    geom_jitter(height = 0.2, width = 0.1, color = "gray") +
    labs(x = "Nível da lagoa", y = "Riqueza", title = "Riqueza de Espécies: Lagoa cheia vs. nível baixo") +
    theme_minimal() +
    theme(plot.title = element_text(size = 19, face = "bold", hjust = 0),
          axis.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 16, face = "bold"),
          legend.title = element_blank(),
          legend.text = element_text(size = 13),
          legend.position = "right")
  
}  
    # Salvar em pdf
  ggsave(plot =  boxplot_alfa, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/alfa_d/",
                          "alfa-d_boxplot_cheio-vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 30,
         width = 19,
         dpi = 600) 


## Diagrama de Venn
{
  library("ggvenn", lib.loc="/opt/R/4.1.0/lib/R/library")
  
  
  list_nivel_id <- 
    longer_tbl_IDS %>%
    group_by(Nivel) %>% 
    summarise(
      `Curated ID` = unique(`Curated ID`)) %>% 
    with(split(`Curated ID`, 
               factor(Nivel, levels = unique(Nivel))))
  
  # venn_cheio_vazio <-
    list_nivel_id %>% 
    ggvenn(c("Cheio", "Vazio"),
           fill_color = c("#f8766d", "#00bfc4"),
           text_size = 6,
           set_name_size = 18,
           show_elements = TRUE,
           label_sep = "\n"
           # ,
           # auto_scale = TRUE
           )
  
  ggsave(plot = venn_cheio_vazio, 
         filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/alfa_d/",
                          "venn_cheio-vazio", "-", Sys.Date(), ".pdf", sep = ""),
         device = "pdf",
         units = "cm",
         height = 30,
         width = 45,
         dpi = 600)
  }


## NMDS ----

# NMDS (Non-Metric Multidimensional Scaling), ou Escalonamento Multidimensional Não 
# Métrico, é uma técnica de análise multivariada utilizada para visualizar a semelhança 
# ou dissimilaridade entre objetos ou observações em um conjunto de dados. O NMDS é 
# particularmente útil quando os dados não possuem uma estrutura linear clara e quando 
# a distância entre os pontos não pode ser representada de forma métrica.

# Converter planilha de identificações por ponto para o formato Amostras X IDs
{
  # Função necessária para juntar as abds da tabela em IDs
  sum_uniq <- function(vec) {
    if (is.character(vec)) {
      suniq <- BiocGenerics::unique(vec)
    }
    if (is.numeric(vec)) {
      suniq <- sum(vec)
    }
    return(suniq)
  }
  
  # Número de linhas (IDs diferentes) por amostra
  fish_ID_tbl$File_name %>% table()
  
  fish_ID_tbl %>% colnames()
  
  FINAL_tbl_IDs <- fish_ID_tbl %>%
    # filter(Year %in% c("2021")) %>% # definir qual ano entrara na analise
    filter(RRA >= 0.01) %>% # definicao de qual threshold de abundancia sera usado
    # filter(New_name %in% c("Ponte", "Fundacao")) %>% 
    filter(Sample %in% c("LI1-neo-mi", "L2_nov20", "L2_dez20", "L2_nov21", "L3_nov21", "L4_nov21")) %>% #com apenas as amostras que Daniel pediu
    select(c("Sample", "Curated ID", "Point", "Year", "Nível", "New_name", "RRA")) %>% 
    pivot_wider(id_cols = c("Sample", "Point", "Year", "Nível", "New_name"),
                names_from = `Curated ID`,
                values_from = `RRA`,
                values_fn = sum_uniq,
                names_sort = TRUE,
                names_prefix = "ID_") %>% 
    relocate(c("Sample", "Point", "Year", "Nível", "New_name", starts_with("ID_"))) %>%  
    mutate(across(starts_with("ID_"), replace_na, replace = 0)) %>% 
    mutate("New_name" = as.factor(New_name)) %>% 
    mutate("Nível" = as.factor(Nível)) 
  
  # Verificando a soma das linhas
  FINAL_tbl_IDs %>% select(starts_with(match = "ID_")) %>% rowSums(na.rm = TRUE)
  
  # TODO: Verificando a soma das colunas
}

#  Rodando o NMDS
{
  # Preparar dados para entrada no pacote vegan
  colnames(FINAL_tbl_IDs)
  
  FINAL_tbl_IDs$`Sample` %>% unique() %>% sort()
  
  all_IDs_NMDS_tbl <- FINAL_tbl_IDs %>% 
    mutate("Sample number" = 0) %>% 
    relocate("Sample number" )
  
  # Associar números de amostra aos nomes de amostra
  for (sample in 1:nrow(all_IDs_NMDS_tbl)) {
    all_IDs_NMDS_tbl$`Sample number`[sample] <- sample
  }
  
  # Ordenar dataframe usado no NMDS
  all_IDs_NMDS_df <- all_IDs_NMDS_tbl %>% as.data.frame() 
  
  # Nomear linhas como números de amostra e remover coluna
  row.names(all_IDs_NMDS_df) <- all_IDs_NMDS_df$`Sample number`
  
  # Corrigir nomes das espécies para evitar problemas na plotagem
  colnames(all_IDs_NMDS_df)
  colnames(all_IDs_NMDS_df)[7:ncol(all_IDs_NMDS_df)] <- colnames(all_IDs_NMDS_df)[7:ncol(all_IDs_NMDS_df)] %>%
    str_replace_all(pattern = " ", replacement = "_") %>% 
    str_replace_all(pattern = "\\.", replacement = "") %>% 
    str_replace_all(pattern = "\\(", replacement = "") %>% 
    str_replace_all(pattern = "\\)", replacement = "")
  
  # Executando o NMDS
  # Esta e a funcao que faz o NMDS. Para ela fornece-se apenas
  # as colunas relativas as especies nas amostras
  all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)],
                                   comm = all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)],
                                   distance = "bray" # Leva em conta o RRA
                                   # distance = "jaccard" # Considera apenas presença/ausenciaa
                                   )
    
  dim(all_IDs_NMDS_df)

  # Fazer o fit das variáveis ambientais
  meta.envfit <- envfit(all_ps_vegan_ord_meta, all_IDs_NMDS_df[, c("Nível", "Sample")],
                        permutations = 999,
                        na.rm=TRUE) # this fits environmental vectors
  
  # Espécies 
  
  # Fazer o fit das espécies para identificar a significância delas na explicação dos agrupamentos. 
  # Esse é o passo que mais demora quando com muitas amostras e espécies.
  meta.spp.fit <- envfit(all_ps_vegan_ord_meta, 
                         all_IDs_NMDS_df[,7:ncol(all_IDs_NMDS_df)], 
                         permutations = 999) # this fits species vectors
  
  # Obter valores de p para as espécies
  sps_pvals <- tibble("IDs" = names(meta.spp.fit$vectors$pvals),
                      "p-value" = meta.spp.fit$vectors$pvals)
  
  spp.scrs <- as.data.frame(scores(meta.spp.fit, display = "vectors")) %>%
    mutate("IDs" = rownames(.)) %>%
    left_join(y = sps_pvals, by = "IDs")
  
  # Selecionar espécies significativas
  sig.spp.scrs <- spp.scrs
  # %>%
  #   filter(`p-value` <=
  #            0.05) # definir p-value aqui!
  
  # Pontos amostrais
  
  # Definir os valores de NMDS1 e NMDS2, e os metadados de cada amostra/ponto amostral
  site.scrs <- as.data.frame(scores(all_ps_vegan_ord_meta, display = "sites")) %>%
    mutate("Sample number" = as.double(row.names(.))) %>%
    left_join(y = all_IDs_NMDS_df[, c("Sample number",
                                      "Sample",
                                      "Nível",
                                      "New_name"
                                      )], 
              by = "Sample number")
  
  # Determinar centroides
  scrs <- scores(all_ps_vegan_ord_meta, display = "sites")
  
  cent <- aggregate(scrs ~ Nível, data = site.scrs, FUN = "mean")
  
  # Calcular elipses
  NMDS <- data.frame("MDS1" = all_ps_vegan_ord_meta$points[, 1],
                     "MDS2" = all_ps_vegan_ord_meta$points[, 2],
                     "Nível" = as.factor(all_IDs_NMDS_df$Nível), check.names = FALSE)
  
  NMDS.mean <- aggregate(NMDS[, 1:2], list(group = NMDS$Nível), "mean")
}

# Elipses

{
  plot(all_ps_vegan_ord_meta)
  
  # Sobrepor as elipses
  ord <- ordiellipse(ord = all_ps_vegan_ord_meta, 
                     groups = all_IDs_NMDS_df$Nível,
                     display = "sites",
                     kind = "ehull", conf = 0.95, label = T)
  
  
  # Funcao do vegan de calcular elipses
  veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
  
  df_ell <- data.frame()
  
  for(g in levels(NMDS$Nível)){
    print(g)
    df_ell <- 
      rbind(df_ell, 
            cbind(as.data.frame(with(NMDS[NMDS$Nível==g,],
                                                     veganCovEllipse(
                                                       ord[[g]]$cov,
                                                       ord[[g]]$center,
                                                       ord[[g]]$scale))),
                  Nível=g))}
}

# Configurar plot NMDS

# Definir paleta de cores
{
  paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")
  paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[1:3]
  paletteer::paletteer_d("colorBlindness::SteppedSequential5Steps")[c(1:3, 2,7,9,10)]
  
  my_cols <- c("BAM" = "#0F6B99FF",
               "PARAC" = "#6B990FFF",
               "SAM" = "#99540FFF",
               "SFC" = "#A3CC51FF",
               "SFI" = "#7EC3E5FF",
               "SFM" = "#E5B17EFF") 
}

# Plot
{
  NMDS_cheio_vazio <-
    ggplot(data = site.scrs,
           aes(x=NMDS1, 
               y=NMDS2)) +
    
    # Elipses 
    ggforce::geom_mark_ellipse(inherit.aes = FALSE,
                               data = df_ell,
                               aes(x = NMDS1,
                                   y = NMDS2,
                                   group = Nível,
                                   label = Nível,
                                   col = Nível,
                                   fill = Nível
                               ), 
                               alpha=0.10,
                               # n = 200,
                               linetype=2,
                               expand = 0,
                               label.fontsize = 18,
                               con.cap = 0.1
    ) +
    # Niveis amostrais
    geom_point(aes(x=NMDS1,
                   y=NMDS2,
                   fill = Nível,
                   # label = `Sampling unit`,  #descomentar se quiser exibir os nomes dos `Sampling sites`s
                   col=Nível,
                   group = Nível,
                   shape = New_name),
               stroke = 0.5,
               alpha = 0.75,
               size = 3) +
    
    # Nomes dos `Sampling sites`s amostrais
    ggrepel::geom_text_repel(aes(label = Sample),  #descomentar bloco se quiser exibir os nomes dos pontos
                             # hjust=0.5,
                             # vjust=2.75,
                             size=4,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha=0.1,
                             # min.segment.length = 1.5,
                             force = 3,
                             max.overlaps = 100) +
    
    # Vetores das IDs
    geom_segment(data = sig.spp.scrs, aes(x = 0,
                                          xend = NMDS1,
                                          y = 0,
                                          yend = NMDS2),
                 arrow = arrow(length = unit(0.1, "cm")),
                 colour = "grey30",
                 alpha = 0.5,
                 lwd = 0.3) + #add vector arrows of significant species
    
    # Nomes das IDs
    ggrepel::geom_text_repel(data = sig.spp.scrs,
                             aes(x=NMDS1, y=NMDS2, label = IDs),
                             size = 3.5,
                             alpha = 0.75,
                             direction = "both",
                             segment.size = 0.25,
                             segment.alpha = 0.1,
                             max.overlaps = 100) +
    
    # Centroides
    geom_point(data = cent,
               aes(x = NMDS1,
                   y = NMDS2, # colour = Nível,
                   fill = Nível
               ),
               size = 8,
               colour = "#222222",
               alpha = 0.75,
               shape  = 13
    ) + 
    coord_fixed(expand = c(0.5))+
    theme_light() +
    # +
    # scale_colour_manual(values = c("#fcca03","#6b0000","#02cc37"))+
    # scale_colour_manual(values = c("#fcca03","#6b0000","#02cc37"))+
    # scale_fill_manual(values = my_cols)+
    # scale_colour_manual(values = my_cols) +
    # +
    # scale_shape_manual(values = c("MG" = 24,                                    # se quiser formas customizadas para algum metadado
    #                               "RJ" = 21,
    #                               "SP" = 22,
    #                               "Minas Gerais" = 24,
    #                               "Rio de Janeiro" = 21,
  #                               "São Paulo" = 22
  #                               )) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        axis.text = element_text(size = 10)) +
    labs(
      title  = paste0("Composição da ictiofauna: lagoa vazia e lagoa cheia"),
      subtitle = paste0("Stress: ",format(round(all_ps_vegan_ord_meta$stress,4)))
    ) +
    theme(plot.title = element_text(size = 20, face = "bold")) +
    theme(plot.subtitle = element_text(size = 16)) +
    theme(legend.title = element_text(size = 16)) +
    theme(legend.text =  element_text(size = 12)) +
    theme(axis.title = element_text(size = 18),
          axis.text = element_text(size = 16)) +
    theme(legend.position = "bottom") +
    guides(shape= "none")
  
  NMDS_cheio_vazio
}

# Salvar em pdf
ggsave(plot =  NMDS_cheio_vazio, 
       filename = paste("/home/gabriel/projetos/lagoa_ingleses/results/figuras/2023/nmds/",
                        "nmds_cheio_vazio", "-", Sys.Date(), ".pdf", sep = ""),
       device = "pdf",
       units = "cm",
       height = 20,
       width = 30,
       dpi = 600)
                    
                        


