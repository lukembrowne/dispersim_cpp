
library(gridExtra)

abund_thresh = 5 ## Lowest abundance to exclude

## SGDC scatter plot - averaging across GNDD
sgdcs <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  group_by(mean_gndd) %>%
  mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate SGDC
  # group_by(mean_gndd) %>%
  # summarise_all(mean) %>%
  ggplot() +
    geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), size = I(5), pch = 19) + 
    geom_hline(aes(yintercept = 0)) + theme_bw() +
    ylab("R (SGDC)") + xlab(expression('b'['gen']))

sgdcs


## Min GNDD
min_gndd <- dat_sp %>%
    filter(abundance >= abund_thresh) %>% # only species with positive abundance
     filter(mean_cndd != 1) %>% # Filter out neutral models
     filter(mean_gndd == min(mean_gndd)) %>%
    group_by(sim_id) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = sp_shannon, y = allelic_shannon, fill = mean_cndd)) + 
    geom_point(shape = 21, colour = "black", size = 2, stroke = 1) + 
    theme_bw() + 
    scale_fill_distiller(palette = "Spectral", 
                       limits = c(0, 1.0)) + 
    xlab("Species Shannon's diversity") + ylab("Genotype Shannon's diversity") +
  labs(fill = expression('b'['con']))

min_gndd

## Max GNDD
max_gndd <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models
   filter(mean_gndd == max(mean_gndd)) %>%
  group_by(sim_id) %>%
  summarise_all(mean) %>%
  ggplot(aes(x = sp_shannon, y = allelic_shannon, fill = mean_cndd)) +
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) + 
  theme_bw() + 
  scale_fill_distiller(palette = "Spectral", 
                       limits = c(0, 1.0)) + 
  xlab("Species Shannon's diversity") + ylab("Genotype Shannon's diversity") + 
  labs(fill = expression('b'['con']))

max_gndd


### Save to file and edit in inkscape

# pdf("./figures/Figure 2a - SGDCs by NDD.pdf", width = 6, height = 4)
# sgdcs
# dev.off()
# 
# pdf("./figures/Figure 2b - SGDCs by NDD.pdf", width = 5, height = 4)
# min_gndd
# dev.off()
# 
# pdf("./figures/Figure 2c - SGDCs by NDD.pdf", width = 5, height = 4)
# max_gndd
# dev.off()


pdf("./figures/Figure 1 - SGDCs due to conspecific NDD.pdf", width = 12, height = 4)
grid.arrange(min_gndd, sgdcs, max_gndd, ncol = 3)
dev.off()



# Comparing across variation in gndd --------------------------------------

## SGDC scatter plot - averaging across CNDD
sgdcs <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_gndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  group_by(mean_cndd) %>%
  mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate SGDC
  ggplot() +
  geom_point(aes(x = factor(round(mean_cndd, 2)), y = sgdc), size = I(5), pch = 19) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() +
  ylab("R (SGDC)") + xlab(expression('b'['con']))

sgdcs

## Min CNDD
min_cndd <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_gndd != 1) %>% # Filter out neutral models
  filter(mean_cndd == min(mean_cndd)) %>%
  group_by(sim_id) %>%
  summarise_all(mean) %>%
  ggplot(aes(x = sp_shannon, y = allelic_shannon, fill = mean_gndd)) + 
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) + 
  theme_bw() +
  xlab("Species Shannon's diversity") + ylab("Genotype Shannon's diversity") +
  labs(fill = expression('b'['gen']))

min_cndd

## Max CNDD
max_cndd <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_gndd != 1) %>% # Filter out neutral models
  filter(mean_cndd == max(mean_cndd)) %>%
  filter(mean_cndd != 0.9) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>%
  ggplot(aes(x = sp_shannon, y = allelic_shannon, fill = mean_gndd)) +
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) + 
  theme_bw() + 
  xlab("Species Shannon's diversity") + ylab("Genotype Shannon's diversity") + 
  labs(fill = expression('b'['gen']))

max_cndd

pdf("./figures/Figure S - SGDCs due to genotypic NDD.pdf", width = 12, height = 4)
grid.arrange(min_cndd, sgdcs, max_cndd, ncol = 3)
dev.off()




# Pop size ----------------------------------------------------------------

single_species <- dat_sp %>%
  #filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1.0) %>% # Filter out neutral models
  filter(sp == 0 | sp == 49) %>% # Just most abundant species
  filter(mean_gndd == 1) %>% # Filter out neutral models
  ggplot(aes(x = factor(mean_cndd), y = abundance, fill = factor(sp))) + 
  geom_boxplot() + 
  geom_jitter(shape = 21, colour = "black", size = 2, stroke = 1, height = 0) + 
  theme_bw() + 
 # scale_fill_distiller(palette = "Spectral", 
 #                      limits = c(0, 1.0)) + 
  xlab("Strength of NFD") + ylab("Abundance") +
  labs(fill = expression('b'['con']))

single_species

pdf("./figures/Population size single species.pdf", height = 5, width =8)
single_species
dev.off()

### SGDC single species
single_species_sgdc <- dat_sp %>%
  #filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(sp == 0 | sp == 49) %>%
  filter(mean_cndd != 1) %>% # Filter out neutral models
  filter(mean_gndd == 1.0) %>%
  ggplot(aes(x = sp_shannon, y = allelic_shannon, fill = factor(sp))) + 
  geom_point(shape = 21, colour = "black", size = 2, stroke = 1) + 
  theme_bw() + 
  xlab("Species Shannon's diversity") + ylab("Genotype Shannon's diversity") +
  labs(fill = expression('b'['con']))

single_species_sgdc

pdf("./figures/SGDC single species.pdf", height = 3, width =6)
single_species_sgdc
dev.off()

