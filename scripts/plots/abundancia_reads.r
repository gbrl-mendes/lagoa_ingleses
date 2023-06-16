---
title: "Análises de Abundancia"
author: "Gabriel Mendes"
date: "22/03/2023"
---
  options(scipen=100)
    grouped_by_ID_BLASTid %>%
      mutate(`Year` = factor(`Year`)) %>% 
      # filter(Year == 2021) %>%
      group_by(Year, Class) %>%
      mutate('Class' = ifelse(`Class` == "" | `Class` == "Actinobacteria",
                                    "Sem identificação", `Class`)) %>%
      summarize(prop = Abundance / sum(Abundance) * 100) %>%
      ungroup() %>% 
    ggplot(aes(x = Year, y = prop, fill = Class)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_y_log10() +
    labs(x = "Year", y = "Proportion (%)", fill = "Class") 
    
    
    # Dividindo por ano
    
    
abundance_curated_ID_year <-
  grouped_by_ID_BLASTid %>%
  group_by(Year) %>%
  mutate("abd total year" =  sum(Abundance)) %>%
  ungroup() %>%
  group_by(Year, `Curated ID`) %>%
  summarize(prop = sum(Abundance) / unique(`abd total year`)*100) %>%
  # mutate('Class' = ifelse(`Class` == "" | `Class` == "Actinobacteria",
  #                         "Sem identificação", `Class`)) %>%
  mutate(`Year` = factor(`Year`)) %>%
  ungroup() %>% 
      ggplot(aes(x = Class, y = log10(prop), fill = Class)) +
      geom_bar(stat = "identity", position = "stack") +
      labs(x = "Year", y = "Proportion (%)", fill = "Class") +
      # scale_y_log10(breaks=c(0.01,0.1,1,10,25,50,100)) +
      facet_grid2(cols = vars(Year))
  
    # Abundancia ASVs

  raw_results_tbl %>% 
    group_by(Class) %>% 
    filter(Year == 2020) %>% 
    summarise("Num ASVs" = length(unique(`ASV (Sequence)`)))









