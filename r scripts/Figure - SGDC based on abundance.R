
abund_thresh = 5 ## Lowest abundance to exclude

## Based on SPECIES ID - CNDD
sgdc_by_sp <- dat_sp %>%
    filter(abundance >= abund_thresh) %>% # only species with positive abundance
    filter(mean_cndd != 1) %>% # Filter out neutral models
  #  filter(mean_gndd == 1.0) %>%
    mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
    group_by(sp, mean_gndd) %>%
    mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
    summarise_all(mean) %>% ## To get average abundance_rank 
    mutate(gndd_label = paste("Bgen =", mean_gndd)) %>%
    ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
    geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
    geom_hline(aes(yintercept = 0)) + facet_grid(. ~ gndd_label) + 
    labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance")

sgdc_by_sp

pdf("./figures/Figure 5 - SGDC by sp.pdf", height = 4, width = 6)
sgdc_by_sp
dev.off()


#########
### Showing changes in sgdc by changes in GNDD by species

# sgdc_by_sp <- dat_sp %>%
#   filter(abundance >= abund_thresh) %>% # only species with positive abundance
#   filter(mean_cndd != 1) %>% # Filter out neutral models
#   mutate(sp = sp + 1) %>%
#   group_by(sp, mean_gndd) %>%
#   mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
#   summarise_all(mean) %>%
#   #filter(mean_gndd == 0.05 | mean_gndd == 1.0) %>%
#   group_by(sp) %>%
#   mutate(sgdc_pos = ifelse(sgdc > 0, "Yes", "No")) %>%
#   mutate(sgdc_diff = sgdc[mean_gndd == 0.15] - sgdc[mean_gndd == 1.0]) %>%
#   mutate(sgdc_strength = ifelse(sgdc_pos == "Yes" & sgdc_diff > 0,
#                                 "Strengthen", ifelse(sgdc_pos == "No" & sgdc_diff > 0, "Strengthen", "Weaken"))) %>%
#   ggplot(aes(x = sp, y = abs(sgdc_diff), fill = factor(sgdc_strength))) + 
#   geom_point(size = 5, pch = 21) + 
#   scale_fill_manual(
#     values = c("Strengthen" = "steelblue2","Weaken" = "black")) + 
#   theme_bw() + 
#   geom_hline(aes(yintercept = 0))   + 
#   labs(x = "Species ID", y = "Difference in R(SGDC)", fill = "Positive SGDC \nunder no genotypic NDD")
# 



sgdc_by_sp <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models
  mutate(sp = sp + 1) %>%
  group_by(sp, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  filter(mean_gndd == 0.05 | mean_gndd == 1.0) %>%
  summarise_all(mean) %>%
  # group_by(sp) %>%
  mutate(sgdc_pos = ifelse(sgdc > 0, "Yes", "No")) %>%
  mutate(sgdc_diff = sgdc[mean_gndd == 1.0] - sgdc[mean_gndd == 0.05]) %>%
  mutate(sgdc_strength = ifelse(sgdc_pos == "Yes" & sgdc_diff < 0,
                                "Positive", ifelse(sgdc_pos == "No" & sgdc_diff < 0, "Positive", "Negative"))) %>%
  ggplot(aes(x = sp, y = abs(sgdc_diff), fill = factor(sgdc_strength))) + 
  geom_point(size = 5, pch = 21) + 
  scale_fill_manual(
    values = c("Positive" = "steelblue2","Negative" = "black")) + 
  theme_bw() + 
  geom_hline(aes(yintercept = 0))   + 
  labs(x = "Species ID", y = "Difference in R(SGDC)", fill = "Direction of change \nfrom no genotypic NDD \nto strong genotypic NDD")


sgdc_by_sp

pdf("./figures/Figure S - changes in SGDC by sp.pdf", height = 5, width = 8)
sgdc_by_sp
dev.off()



# sgdc_by_sp <- dat_sp %>%
#   filter(abundance >= abund_thresh) %>% # only species with positive abundance
#   filter(mean_cndd != 1) %>% # Filter out neutral models
#   #mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
#   mutate(sp = sp + 1) %>%
#   group_by(sp, mean_gndd) %>%
#   mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
#   summarise_all(mean) %>% ## To get average abundance_rank 
#   mutate(gndd_label = paste("Bgen =", mean_gndd)) %>%
#   ggplot(aes(x = mean_gndd, y = sgdc, fill = abundance)) + 
#   geom_point(aes(size = abundance), pch = 21) + 
#   #geom_smooth() + 
#   theme_bw() + scale_fill_distiller(palette = "Spectral") +
#   geom_hline(aes(yintercept = 0))  + facet_wrap(~ sp) + 
#   labs(x = "Species ID", y = "R(SGDC)", fill = "Abundance")
# 
# sgdc_by_sp
# 
# pdf("./figures/Figure S - changes in SGDC by sp.pdf", height = 12, width = 16)
# sgdc_by_sp
# dev.off()
# 


## Based on SPECIES ID - GNDD
sgdc_by_sp <- dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_gndd != 1) %>% # Filter out neutral models
  mutate(mean_cndd = factor(round(mean_cndd, 2))) %>%
  group_by(sp, mean_cndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  mutate(gndd_label = paste("Bcon =", mean_cndd)) %>%
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + facet_grid(. ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance")

sgdc_by_sp

pdf("./figures/Figure S - SGDC by sp - cndd.pdf", height = 4, width = 12)
sgdc_by_sp
dev.off()




## For individual species

dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
   filter(mean_cndd != 1) %>% # Filter out neutral models
  filter(sp %in% c(0, 5, 10, 15, 49)) %>%
  group_by(mean_gndd) %>%
  ggplot(aes(x = sp_shannon, y =allelic_shannon, fill = mean_cndd, size = abundance)) + 
  geom_point(pch = 21) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_smooth(method='lm',formula=y~x) + 
  facet_grid(sp ~ mean_gndd)  

## GNDD Vs abundance
dat_sp %>%
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models
  filter(sp %in% c(0, 5, 10, 15, 49)) %>%
  group_by(mean_gndd) %>%
  ggplot(aes(x = mean_cndd, y =abundance, fill = mean_cndd, size = abundance)) + 
  geom_point(pch = 21) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_smooth(method='lm',formula=y~x) + 
  facet_grid(sp ~ mean_gndd)  




