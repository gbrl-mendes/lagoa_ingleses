### NMDS

library(vegan)
# data(dune)  
# decorana(dune)

# class(dune)
#1- prepare data for entry in vegan ----

colnames(smp_abd_ID)

raw_results_tbl

all_IDs_NMDS_tbl <- raw_results_tbl %>% 
  
  # mutate(`Curated ID` = factor(`Curated ID`)) %>% 
  select(c(Sample,Primer,Read_origin, `Curated ID`,`Relative abundance on sample`)) %>% 
  filter(Read_origin %in% c("merged")) %>%
  unite(col = "Sample_name_read", c(Sample, Primer, Read_origin), remove = FALSE ) %>% 
  group_by(Sample_name_read, Sample, `Curated ID`,Primer,Read_origin
           # ,Run
  ) %>%
  summarise(`Relative abundance on sample (%)` = sum(`Relative abundance on sample`)) %>% 
  pivot_wider(c(Sample_name_read,Sample,Primer,Read_origin
                # ,Run
  ),names_from =  `Curated ID` ,values_from = `Relative abundance on sample (%)`) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  mutate("Sample number" = 0) %>% 
  ungroup()  %>% 
  select(`Sample number`, 1:(ncol(.)-1)) 

#2- associate sample numbers to sample names ----
for (sample in 1:nrow(all_IDs_NMDS_tbl)) {
  
  all_IDs_NMDS_tbl$`Sample number`[sample] <- sample
  # all_IDs_NMDS_tbl$`Sample number`[sample] <- all_IDs_NMDS_tbl$Sample[sample] %>% stringr::str_remove(pattern = "PP0|PP")
  
}





colnames(all_IDs_NMDS_tbl)
hist(colSums(all_IDs_NMDS_tbl[,-c(1:4)]))
hist(rowSums(all_IDs_NMDS_tbl[,-c(1:4)]))
all_IDs_NMDS_tbl[,-c(1:5)]

all_IDs_NMDS_tbl %>% select(Sample_name_read, `Sample number`) %>% unique()
# all_ps_blst_vegan %>% select(`Sample number`, 1:(ncol(.)-1))

#3- create data.frame of species counts: rownames are Sample numbers ----

all_IDs_NMDS_df <- all_IDs_NMDS_tbl %>% 
  # select(-c("Sample_name_read", "Sample", "Primer","Read_origin"
  # ,"Run"
  # )) %>% 
  select(base::sort(colnames(.))) %>% 
  as.data.frame() 




#4- name rows as Sample numbers and remove column ----
row.names(all_IDs_NMDS_df) <- all_IDs_NMDS_df$`Sample number`

all_IDs_NMDS_df <- all_IDs_NMDS_df %>% 
  select(-c(`Sample number`))
#1- prepare data for entry in vegan ----

colnames(smp_abd_ID)

all_IDs_NMDS_tbl 


colnames(all_IDs_NMDS_tbl)


all_IDs_NMDS_df <- all_IDs_NMDS_tbl %>% 
  select(base::sort(colnames(.))) %>%
  relocate(c("Primer",
             "Sample number",
             "Sample_name_read",
             "Read_origin", 
             "Sample")) %>%
  as.data.frame() 



colnames(all_IDs_NMDS_df)

#4- name rows as Sample numbers and remove column ----
row.names(all_IDs_NMDS_df) <- all_IDs_NMDS_df$`Sample number`

all_IDs_NMDS_df <- all_IDs_NMDS_df %>% 
  select(-c(`Sample number`))


all_IDs_NMDS_df 



# #4- name rows as Sample numbers and remove column ----
# row.names(all_IDs_NMDS_df) <- all_IDs_NMDS_df$`Sample number`
# 
# all_IDs_NMDS_df <- all_IDs_NMDS_df %>% 
#   select(-c(`Sample number`))


colnames(all_IDs_NMDS_df)[6:77] <- colnames(all_IDs_NMDS_df)[6:77] %>% str_replace_all(pattern = " ",replacement = "_")


library(vegan)

all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df[,6:77], comm = all_IDs_NMDS_df[,6:77])
# all_ps_vegan_ord_meta <- prcomp(x = all_IDs_NMDS_df[,6:53])

# all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df[1:9,5:47], comm = all_IDs_NMDS_df[1:9,5:47])
# all_ps_vegan_ord_meta <- metaMDS(veg = all_IDs_NMDS_df[1:9,5:47])


# actually autotransform = FALSE doesn't seem to change the results


all_ps_vegan_ord_meta %>% str()
all_ps_vegan_ord_meta %>% summary()
all_ps_vegan_ord_meta

all_ps_vegan_ord_meta$stress


all_ps_vegan_ord_meta$points



#6b- extract NMDS scores from results

all_vegan_meta <- (vegan::scores(all_ps_vegan_ord_meta) %>% 
                     tidyr::as_tibble(rownames = "Sample number")) %>% 
  mutate(`Sample number` = as.numeric(`Sample number`))
# all_vegan_meta <- as.data.frame(vegan::scores(all_ps_vegan_ord_meta))

#Using the scores function from vegan to extract the site scores and convert to a data.frame

# all_vegan_meta$`Sample number` <- rownames(all_vegan_meta) %>% as.numeric()  

# all_vegan_meta %>% left_join()# create a column of site names, from the rownames of data.scores

# all_vegan_meta <- all_vegan_meta  %>% as_tibble() # create a column of site names, from the rownames of data.scores

#7- bring NMDS scores to complete table

all_vegan_meta_tbl <- left_join(x = unique(all_IDs_NMDS_tbl[,c("Sample number", "Sample_name_read", "Sample","Primer")]),
                                y = all_vegan_meta, by = "Sample number") %>% 
  # mutate(Primer=factor(Primer,levels = c("NeoFish", "MiFish", "COI")),
  # Sample = as.factor(Sample)) %>% 
  select(-c("Sample number"))


library(ggord)
all_IDs_NMDS_tbl$Sample



summary(all_ps_vegan_ord_meta)

all_ps_vegan_ord_meta$species


nmds_PLOT_ord <- ggord(all_ps_vegan_ord_meta, 
                       grp_in = all_IDs_NMDS_tbl$Point, 
                       ellipse = F,
                       size = 10,
                       arrow = 0.5, veccol = "dark grey",
                       txt = 3,
                       repel = T,
                       max.overlaps = 55
                       # ,
                       # # facet=T,
                       # cols = viridis::viridis(option = "turbo",n = nrow(all_IDs_NMDS_df), alpha = 1)
)+
  annotate(geom = "text",
           x=c(0.2),
           y=c(0.20),
           label=c(paste0("Stress: ",format(round(all_ps_vegan_ord_meta$stress,4)))),
           size=5) +
  scale_colour_manual(name = "Samples", values = viridis::viridis(option = "turbo",n = nrow(all_IDs_NMDS_df), alpha = 1))
# + 
#   # labs(title='NEW LEGEND TITLE') +
#   # labs(fill='NEW LEGEND TITLE') +
#   # labs(colour='NEW LEGEND TITLE') + 
#   labs(shape='NEW LEGEND TITLE') 
#   

nmds_PLOT_ord

# nmds_PLOT_ord$guides$shape$title <- "Amostras"
# nmds_PLOT_ord$guides$colour$title <- "Amostras"
nmds_PLOT_ord$guides$colour$title



nmds_PLOT_ord


ggsave(file = paste0(figs_path,"/",prjct_rad,"-",Sys.Date(),"-NMDS_merged.pdf",collapse = ""),
       plot = nmds_PLOT_ord,
       device = "pdf",
       width = 26,
       height = 12,
       dpi = 600)




writexl::write_xlsx(x = all_IDs_NMDS_df,
                    path = paste0(results_path,"/",prjct_rad,"-merged_SPs_by_samples",Sys.Date(),".xlsx"),
                    col_names = TRUE,format_headers = TRUE)


```