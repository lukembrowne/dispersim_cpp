
### Histogram of SGDC values by species
sgdc_hist <- dat_sp %>%
  filter(abundance > abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models
  # group_by(abundance_rank, mean_gndd) %>% ## Based on abundance rank
  group_by(sp, mean_gndd) %>% ## Based on species ID
  mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>%
  ggplot(aes(x = sgdc)) + 
  geom_histogram(col = "black", fill = "steelblue2") +
  xlim(-1, 1) + 
  theme_bw() + 
  labs(x = "R(SGDC)", y = "Frequency") + geom_vline(xintercept = 0, lty = 2)

sgdc_hist

pdf("./figures/Figure 4 - SGDC histogram.pdf", height = 3, width = 5)
sgdc_hist
dev.off()




### Plotting species abundance by species ID

sp_by_rank <- dat_sp %>%
  filter(abundance > abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1 & mean_gndd != 1) %>%
  ggplot(aes(x = sp + 1, y = abundance_rank)) + geom_point() + geom_smooth() + 
  theme_bw() + 
  labs(x = "Species ID", y = "Rank abundance")
 
sp_by_rank

pdf("./figures/Figure S - species id by rank abundance.pdf", width = 5, height = 5)
sp_by_rank
dev.off()


## Calculating correlation
temp <- dat_sp %>%
  filter(abundance > abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1 & mean_gndd != 1)

cor.test(temp$sp, temp$abundance_rank)






