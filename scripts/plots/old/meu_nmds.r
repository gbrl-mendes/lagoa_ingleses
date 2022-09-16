
# DCA plot
{
  # DCA significa Detrended Correspondence Analysis.
  # Técnica de estatística Multivariada que permite clusterizar grupos de variaveis
  # em conjuntos de dados muito grandes e com gradiente, permitindo observar tendencias

  # Apenas peixes
  {
    NMDS_dec_f <- decorana(veg = spp_sample_tbl_f[,5:36])
    NMDS_dec_f %>% summary()
    NMDS_dec_f %>% str()
    NMDS_dec_f$cproj
    
    plot(NMDS_dec_f)
    plot(NMDS_dec_f,type = "p")
    plot(NMDS_dec_f,type = "c") 
    
    points(NMDS_dec_f, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
    text(NMDS_dec_f, display = "sites", cex=0.7, col="blue")
    text(NMDS_dec_f, display = "spec", cex=0.7, col="blue")
  }
  
  # Todas ids
  {
    NMDS_dec_all <- decorana(veg = spp_sample_tbl_all[,5:44])
    NMDS_dec_all %>% summary()
    NMDS_dec_all %>% str()
    NMDS_dec_all$cproj
    
    plot(NMDS_dec_all)
    plot(NMDS_dec_all,type = "p")
    plot(NMDS_dec_all,type = "c") 
    
    points(NMDS_dec_all, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
    text(NMDS_dec_all, display = "sites", cex=0.7, col="blue")
    text(NMDS_dec_all, display = "spec", cex=0.7, col="blue")
  }
  
  # Todas ids e samples
  {
    NMDS_dec_tudo <- decorana(veg = spp_sample_tbl_tudo[,5:51])
    NMDS_dec_tudo %>% summary()
    NMDS_dec_tudo %>% str()
    NMDS_dec_tudo$cproj
    
    plot(NMDS_dec_tudo)
    plot(NMDS_dec_tudo,type = "p")
    plot(NMDS_dec_tudo,type = "c") 
    
    points(NMDS_dec_tudo, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
    text(NMDS_dec_tudo, display = "sites", cex=0.7, col="blue")
    text(NMDS_dec_tudo, display = "spec", cex=0.7, col="blue")
  }
  
}

# NMDS Plot
{
  # Apenas peixes
  {
    # Calcular distancias
    NMDS_meta_f <- metaMDS(veg = NMDS_dec_f, comm = spp_sample_tbl_f[,5:36])
    plot(NMDS_meta_f, type = "t")
    NMDS_meta_f %>% str()
    NMDS_meta_f
    plot(NMDS_meta_f, type = "t")
    plot(NMDS_meta_f, type = "p")
    
    NMDS_meta_tudo$stress
    
    print(nmds_plot_f <- 
            ggord(NMDS_meta_f,
                  grp_in = spp_sample_tbl_f$Point,
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
                     label=c(paste0("Stress: ",format(round(NMDS_meta_f$stress,4)))),
                     size=5
            ) +
            scale_shape_manual(name = "Expedition",
                               values = spp_sample_tbl_f$Shape
            ))
    
    ### Tirar as especies com ids muito genericas e fazer esse grafico apenas
    ### para 2021 separando os clusters por pontos coletados (L1, L2, L3 e L4)
    
    
    nmds_PLOT <- all_vegan_meta_tbl %>% 
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
}
  
  # Todas ids
  {
    # Calcular distancias
    NMDS_meta_all <- metaMDS(veg = NMDS_dec_tudo, comm = spp_sample_tbl_all[,5:51])
    plot(NMDS_meta_all, type = "t")
    NMDS_meta_all %>% str()
    NMDS_meta_all
    plot(NMDS_meta_all, type = "t")
    plot(NMDS_meta_all, type = "p")
    
    NMDS_meta_all$stress
    
    print(nmds_plot_all <- 
            ggord(NMDS_meta_all,
                  grp_in = spp_sample_tbl_all$Point,
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
                     label=c(paste0("Stress: ",format(round(NMDS_meta_all$stress,4)))),
                     size=5
            ) +
            scale_shape_manual(name = "Expedition",
                               values = spp_sample_tbl_all$Shape
            ))
    
    ### Tirar as especies com ids muito genericas e fazer esse grafico apenas
    ### para 2021 separando os clusters por pontos coletados (L1, L2, L3 e L4)
    
    
    nmds_PLOT <- all_vegan_meta_tbl %>% 
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
}
  
  # Todas ids e samples
  {
    # Calcular distancias
    NMDS_meta_tudo <- metaMDS(veg = NMDS_dec_tudo, comm = spp_sample_tbl_tudo[,5:51])
    plot(NMDS_meta_tudo, type = "t")
    NMDS_meta_tudo %>% str()
    NMDS_meta_tudo
    plot(NMDS_meta_tudo, type = "t")
    plot(NMDS_meta_tudo, type = "p")
    
    NMDS_meta_tudo$stress
    
    print(nmds_plot_tudo <- 
            ggord(NMDS_meta_tudo,
                  grp_in = spp_sample_tbl_tudo$Point,
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
                     label=c(paste0("Stress: ",format(round(NMDS_meta_tudo$stress,4)))),
                     size=5
            ) +
            scale_shape_manual(name = "Expedition",
                               values = spp_sample_tbl_tudo$Shape
            ))
    
    ### Tirar as especies com ids muito genericas e fazer esse grafico apenas
    ### para 2021 separando os clusters por pontos coletados (L1, L2, L3 e L4)
    
    
    nmds_PLOT <- all_vegan_meta_tbl %>% 
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
  }
}
  


  