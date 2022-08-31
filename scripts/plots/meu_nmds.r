
# Detrended Correspondence Analysis Plot 
{
  # Técnica de estatística Multivariada que permite clusterizar grupos de variaveis
  # em conjuntos de dados muito grandes e com gradiente, permitindo observar tendencias

  NMDS_dec <- decorana(veg = spp_sample_tbl[,4:35])
  NMDS_dec %>% summary()
  NMDS_dec %>% str()
  NMDS_dec$cproj
  
  plot(NMDS_dec)
  plot(NMDS_dec,type = "p")
  plot(NMDS_dec,type = "c") 
  
  points(NMDS_dec, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
  text(NMDS_dec, display = "sites", cex=0.7, col="blue")
  text(NMDS_dec, display = "spec", cex=0.7, col="blue")
}

# NMDS analisys 
{
  # Calculate distances
  
  NMDS_meta <- metaMDS(veg = spp_sample_tbl[,5:36], comm = spp_sample_tbl[,5:36])
  # actually autotransform = FALSE doesn't seem to change the results
  plot(NMDS_meta, type = "t")
  NMDS_meta %>% str()
  NMDS_meta
  plot(NMDS_meta, type = "t")
  plot(NMDS_meta, type = "p")
  
  NMDS_meta$stress
  
  # nmds_plot <- 
  ggord(NMDS_meta, 
        grp_in = spp_sample_tbl$Point, 
        ellipse = F,
        size = 10,
        arrow = 0.5, veccol = "dark grey",
        txt = 3,
        repel = T,
        max.overlaps = 55
        ) +
    annotate(geom = "text",
             x=c(0.2),
             y=c(0.20),
             label=c(paste0("Stress: ",format(round(NMDS_meta$stress,4)))),
             size=5
             ) +
  scale_shape_manual(name = "Expedition",
                     values = spp_sample_tbl$form)
  # +
  # scale_fill_manual(name = "Expedition",
  #                   values = Exp_point_colors)
  # +
  # scale_color_manual(name = "Expedition2",
  #                    values = Exp_point_colors)
  
  
  ### Tirar as especies com ids muito genericas e fazer esse grafico apenas
  ### para 2021 separando os clusters por pontos coletados (L1, L2, L3 e L4)
  
  
  nmds_PLOT <- all_vegan_meta_tbl %>% 
    # filter(Run %in% c("LGC_MiniSeq_1", "LGC_MiniSeq_2")) %>% 
    ggplot(aes(x = NMDS1,y = NMDS2, col = Sample,
               shape = Year,
               label = Sample,
               Group = Year))+
    geom_point(size = 11)+
    theme_linedraw(base_size = 18) +
    theme(legend.position="bottom") +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = viridis::viridis(option = "turbo",n = 30, alpha = 1, begin = 0, end = 1, direction = 1)
    ) +
    annotate(geom = "text",
             x=c(0.30),y=c(-0.3),label=c(paste0("Stress: ",fact_NMDS_ps_vegan_ord_meta$stress)),size=5)  
  
  nmds_PLOT
  
  
  ggsave(file = paste0(figs_path,"/",prjct_radical,"_NMDS.pdf"),
         plot = nmds_PLOT,
         device = "png",
         width = 31,
         height =20,
         units = "cm",
         dpi = 300)