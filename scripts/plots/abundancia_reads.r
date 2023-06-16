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

## Tabela inf by samples ----
  
  # Agrupamento da ASVs que possuem os mesmos atributos abaixo
  # incluindo as seqs das ASVs
  {
    inf_grouped_by_ID_BLASTid <- raw_results_tbl %>%
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
        "Abundance",
        "ASV (Sequence)",
        "OTU"
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
                `OTU`,
                `ASV (Sequence)`,
                `Abundance`,
                `RRA` = sum(`Relative abundance on sample`)) %>%
      ungroup()
  }
  
  # tabela so peixes
  # incluindo as seqs das ASVs
  {
    inf_fish_ID_tbl <- inf_grouped_by_ID_BLASTid %>% # ids apenas os peixes a nivel de spp
      mutate(`Curated ID` = factor(`Curated ID`,
                                   levels = rev(spp_fish))) %>% 
      mutate(`Order` = factor(`Order`,
                              levels = fish_ordens)) %>%
      mutate('Nível' = ifelse(`Year` == "2020", ## com essa linha a gente inclui o nivel da lagoa
                              "Cheio", "Vazio")) %>%
      filter(`Curated ID` != "") %>%
      filter(!is.na(`Curated ID`)) %>% 
      arrange(`Curated ID`)
  }
  
  inf_tbl <-
    inf_fish_ID_tbl %>%
    mutate(OTU = as.character(OTU)) %>% 
    group_by(Sample, New_name, Year) %>%
    filter(RRA >= 0.01) %>% 
    reframe(Local = New_name,
            Ano = Year,
            Nível,
            Reads = sum(Abundance),
            ASVs = length(unique((`ASV (Sequence)`))),
            OTUs = length(unique(OTU)),
            Spp = length(unique(`Curated ID`)),
    ) %>%
    unique() %>% 
    rename(Amostras = Sample) %>%
    arrange(match(Amostras, c("LI1-neo-mi", "L2_nov20", "L2_dez20", "L1_out21", 
                              "L2_out21", "L3_out21", "L4_out21", "L1_nov21", 
                              "L2_nov21", "L3_nov21",  "L4_nov21"))) 

