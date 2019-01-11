## Setwd, load libraries, and load functions
setwd("~/Dropbox/dispersim_cpp")

library(tidyverse)
library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors
library(viridis)
library(gridExtra)

source("./r scripts/functions.R")

## Set run name
run_name <- "run_031"

## Load in data
dat_sp_sens = read_summary_out(path_to_folder = 
                            paste("./data/", run_name,"/summary_out/", sep = ""),
                          type = "by_sp")
dat_sp_sens

## Merge with parameter values
params <- read_csv(paste("./data/", run_name, "/", run_name, ".csv", sep = ""))

dat_sp_sens <- left_join(dat_sp_sens, params, by = "sim_id")

glimpse(dat_sp_sens)

## Rank species by abundance
## If abundance is a tie, will continute sequential counts of ranks
## Eg. 0, 0, 0, would be 97..98..99 ranks rather than all be a tie

## Loop through sim_ids, then through steps, then assign rank abundance

dat_sp_sens$abundance_rank <- NA

for(id in unique(dat_sp_sens$sim_id)){
  
  cat("On sim id:", id, "... \n")
  
  for(step in unique(dat_sp_sens$step)){
    
    indices <- which(dat_sp_sens$sim_id == id & dat_sp_sens$step == step)
    
    # ## Length should equal number of species
    # if(length(indices) != dat_sp_sens$n_sp_init[1]){
    #   cat("STOPPING \n")
    #   stop()
    # }
    
    if(length(indices) == 0){
      next
    }
    
    dat_sp_sens$abundance_rank[indices[order(dat_sp_sens$abundance[indices], decreasing = T)]] <- 1:dat_sp_sens$n_sp_init[indices][1] ## Make sure it goes up to max number  of species
    
  } # End step loop
  
} # End sim_id loop


summary(dat_sp_sens$sim_id)



# Begin plots -------------------------------------------------------------



##################
## Difference in steps

table(dat_sp_sens$steps)

abund_thresh = 5 ## Lowest abundance to exclude

################
## STEPS

## Community_averaged
com_steps <- dat_sp_sens %>%
  
  ### Change values here
  filter(steps == 5000 | steps == 20000) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(steps, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("Steps = ", steps)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Steps")

com_steps

## BY species
sp_steps <- dat_sp_sens %>%
  ### Change values here
  filter(steps == 5000 | steps == 20000) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(steps, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, steps, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Steps = ", steps)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Steps")

sp_steps




#####################
## SP AND ALLELE INT

table(dat_sp_sens$n_sp_init)

## Community_averaged
  com_sp_init <- dat_sp_sens %>%
    
    ### Change values here
    filter(n_sp_init == 25 | n_sp_init == 75) %>% 
    
    filter(abundance >= abund_thresh) %>% # only species with positive abundance 
    filter(mean_cndd != 1) %>% # Filter out neutral models
    group_by(sim_id) %>%
    summarise_all(mean) %>% ## Average across simulation
    
    ## Update grouping here
    group_by(n_sp_init, mean_gndd) %>%
    mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
    
    ## Update label here
    mutate(sens_label = paste("Sp = ", n_sp_init, "; Alleles = ", n_alleles_init)) %>%
    
    ggplot() +
    geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
    geom_hline(aes(yintercept = 0)) + theme_bw() + 
    facet_grid(sens_label ~ .) + 
    ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("No. initial species and alleles")
  
  com_sp_init


## By species
  sp_sp_init <- dat_sp_sens %>%
    ### Change values here
    filter(n_sp_init == 25 | n_sp_init == 75) %>% 
    
    filter(abundance >= abund_thresh) %>% # only species with positive abundance
    filter(mean_cndd != 1) %>% # Filter out neutral models
  
    mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
    
    ## Update grouping here
    group_by(n_sp_init, mean_gndd) %>%
    mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
    
    ## Update grouping here
    group_by(sp, n_sp_init, mean_gndd) %>%
    
    mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
    summarise_all(mean) %>% ## To get average abundance_rank 
    
    ## Update label here
    mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Sp = ", n_sp_init, "; Alleles = ", n_alleles_init)) %>%
    
    ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
    geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
    geom_hline(aes(yintercept = 0)) + 
    #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
    facet_grid(sens_label ~ gndd_label) + 
    labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("No. initial species and alleles")
   
  sp_sp_init


  
  
  
################
## Seed dispersal distance

table(dat_sp_sens$seed_disp_dist)
  
## Community_averaged
com_disp_dist <- dat_sp_sens %>%
  
  ### Change values here
  filter(seed_disp_dist == 2 | seed_disp_dist == 20) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(seed_disp_dist, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("Seed dispersal dist = ", seed_disp_dist)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Seed dispersal distance")

com_disp_dist
 
## By species 
sp_disp_dist <- dat_sp_sens %>%
  
  ### Change values here
  filter(seed_disp_dist == 2 | seed_disp_dist == 20) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(seed_disp_dist, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, seed_disp_dist, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Seed dispersal dist = ", seed_disp_dist)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Seed dispersal distance")

sp_disp_dist




################
## Neighborhood radius

table(dat_sp_sens$neighbor_radius)

## Community_averaged
com_neighbor <- dat_sp_sens %>%
  
  ### Change values here
  filter(neighbor_radius == 10 | neighbor_radius == 30) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(neighbor_radius, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("Neighborhood radius = ", neighbor_radius)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Neighborhood radius")

com_neighbor




## By species
sp_neighbor <- dat_sp_sens %>%
  ### Change values here
  filter(neighbor_radius == 10 | neighbor_radius == 30) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(neighbor_radius, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, neighbor_radius, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Neighborhood radius = ", neighbor_radius)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Neighborhood radius")

sp_neighbor


################
## NDD Range

table(dat_sp_sens$range_cndd)

com_range_cndd <- dat_sp_sens %>%
  
  ### Change values here
  filter(range_cndd == 0.1 | range_cndd == 0.4) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(range_cndd, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("NDD range = ", range_cndd)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Range of NDD values")

com_range_cndd


## By species
sp_range_cndd <- dat_sp_sens %>%
  ### Change values here
  filter(range_cndd == 0.1 | range_cndd == 0.4) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(range_cndd, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, range_cndd, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("NDD range = ", range_cndd)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Range of NDD values")

sp_range_cndd


################
## Seeds per adult

table(dat_sp_sens$seeds_per_adult)

## Community averaged
com_fecundity <- dat_sp_sens %>%
  
  ### Change values here
  filter(seeds_per_adult == 1000 | seeds_per_adult == 10000) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(seeds_per_adult, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("Fecundity = ", seeds_per_adult)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Fecundity")
 
com_fecundity


sp_fecundity <- dat_sp_sens %>%
  ### Change values here
  filter(seeds_per_adult == 1000 | seeds_per_adult == 10000) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(seeds_per_adult, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, seeds_per_adult, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Fecundity = ", seeds_per_adult)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Fecundity")

sp_fecundity




################
## Landscape size

table(dat_sp_sens$height)


## Community averaged
com_landscape <- dat_sp_sens %>%
  
  ### Change values here
  filter(height == 50 | height == 200) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(height, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("Landscape height = ", height)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Landscape size")

com_landscape

## By species
sp_landscape <- dat_sp_sens %>%
  ### Change values here
  filter(height == 50 | height == 200) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(height, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, height, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Landscape height = ", height)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Landscape size")

sp_landscape



################
## Migration rate

table(dat_sp_sens$migration_rate)

## Community averaged
com_migration <- dat_sp_sens %>%
  
  ### Change values here
  filter(migration_rate == 0.01 | migration_rate == 0.00001) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance 
  filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(sim_id) %>%
  summarise_all(mean) %>% ## Average across simulation
  
  ## Update grouping here
  group_by(migration_rate, mean_gndd) %>%
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update label here
  mutate(sens_label = paste("Migration rate = ", migration_rate)) %>%
  
  ggplot() +
  geom_point(aes(x = factor(round(mean_gndd, 2)), y = sgdc), fill = "steelblue2", size = I(5), pch = 21) + 
  geom_hline(aes(yintercept = 0)) + theme_bw() + 
  facet_grid(sens_label ~ .) + 
  ylab("R (SGDC)") + xlab(expression('b'['gen'])) + ggtitle("Migration rate")

com_migration



sp_migration <- dat_sp_sens %>%
  ### Change values here
  filter(migration_rate == 0.01 | migration_rate == 0.00001) %>% 
  
  filter(abundance >= abund_thresh) %>% # only species with positive abundance
  filter(mean_cndd != 1) %>% # Filter out neutral models

  mutate(mean_gndd = factor(round(mean_gndd, 2))) %>%
  
  ## Update grouping here
  group_by(migration_rate, mean_gndd) %>%
  mutate(sgdc_overall = cor(sp_shannon, allelic_shannon)) %>% ## Calculate overall SGDC for horizontal line
  
  ## Update grouping here
  group_by(sp, migration_rate, mean_gndd) %>%
  
  mutate(sgdc = cor(sp_shannon, allelic_shannon)) %>%
  summarise_all(mean) %>% ## To get average abundance_rank 
  
  ## Update label here
  mutate(gndd_label = paste("Bgen =", mean_gndd), sens_label = paste("Migration rate = ", migration_rate)) %>%
  
  ggplot(aes(x = sp, y = sgdc, fill = abundance_rank)) + 
  geom_point(pch = 21, size = 2) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
  geom_hline(aes(yintercept = 0)) + 
  #geom_hline(aes(yintercept = sgdc_overall), lty = 2) + 
  facet_grid(sens_label ~ gndd_label) + 
  labs(x = "Species ID", y = "R(SGDC)", fill = "Rank abundance") + ggtitle("Migration rate")

sp_migration


##### Putting it all together into a master plot(s)


pdf("./figures/Figure S - Model robustness pt1.pdf", height = 10, width = 15)

grid.arrange(com_range_cndd, sp_range_cndd,
             com_fecundity, sp_fecundity,
             ncol = 2, widths = 1:2)

dev.off()

pdf("./figures/Figure S - Model robustness pt2.pdf", height = 10, width = 15)

grid.arrange(com_landscape, sp_landscape,
             com_migration, sp_migration,
             ncol = 2, widths = 1:2)

dev.off()


pdf("./figures/Figure S - Model robustness pt3.pdf", height = 10, width = 15)

grid.arrange(com_neighbor, sp_neighbor,
             com_sp_init, sp_sp_init,
             ncol = 2, widths = 1:2)

dev.off()

pdf("./figures/Figure S - Model robustness pt4.pdf", height = 10, width = 15)

grid.arrange(com_disp_dist, sp_disp_dist,
             com_steps, sp_steps, 
             ncol = 2, widths = 1:2)

dev.off()





