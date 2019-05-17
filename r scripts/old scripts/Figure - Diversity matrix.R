
library(gridExtra)

abund_thresh = 5 ## Lowest abundance to exclude


#### MATRIX PLOTS with text in square
## Species diversity 
species_mat <- dat_sp %>%
    filter(abundance >= abund_thresh) %>% # only species with positive abundance 
    group_by(mean_cndd, mean_gndd) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = mean_gndd, y = mean_cndd, 
               fill = sp_shannon)) +
   scale_x_continuous(breaks = unique(dat_sp$mean_gndd)) +
  scale_y_continuous(breaks = unique(dat_sp$mean_cndd)) +
    geom_tile(colour = "black") + theme_bw()  + scale_fill_distiller(palette = "Spectral", 
                                                     limits = c(0, 4)) + 
    geom_text(aes(x = round(mean_gndd, 2), y = round(mean_cndd, 2), 
                  label = sprintf("%.2f", round(sp_shannon, 2)))) + 
    labs(x = expression('b'['gen']), y = expression('b'['con']),
         fill = "Shannon's index", title = "Species diversity") +
  theme(panel.border = element_blank())

species_mat


## Allelic
allelic_mat <- dat_sp %>%
        filter(abundance >= abund_thresh) %>% # only species with positive abundance
        group_by(mean_cndd, mean_gndd) %>%
        summarise_all(mean) %>%
        ggplot(aes(x = round(mean_gndd, 2), y = round(mean_cndd, 2),
                   fill = allelic_shannon)) +
      scale_x_continuous(breaks = unique(dat_sp$mean_gndd)) +
      scale_y_continuous(breaks = unique(dat_sp$mean_cndd)) +
        geom_tile(colour = "black") + theme_bw()  + scale_fill_distiller(palette = "Spectral", 
                                                         limits = c(0, 4)) +
        geom_text(aes(x = round(mean_gndd, 2), y = round(mean_cndd, 2),
                      label = sprintf("%.2f", round(allelic_shannon, 2)))) + 
  labs(x = expression('b'['gen']), y = expression('b'['con']),
       fill = "Shannon's index", title = "Genetic diversity") + 
  theme(panel.border = element_blank())

allelic_mat

pdf("./figures/Figure 2 - Diversity matrix.pdf", height = 5, width = 15)
grid.arrange(species_mat, allelic_mat, ncol = 2)
dev.off()
