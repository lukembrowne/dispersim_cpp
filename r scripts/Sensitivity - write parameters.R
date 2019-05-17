
### Writing simulation parameters to tab delimited text file so that a BASH 
## script can read this file and set up jobs on cluster

library(tidyverse)
library(lhs)



set.seed(90210)

n_samples <- 100 # agner (1995) demonstrated that a sample of 2000 parameter combinations is adequate for GSA of models with 10–20 parameters; and that additional samples do not significantly improve accuracy. - will ultimately be 10x this if there is a sequence of 10 values along CNDD

lhs_raw <- randomLHS(n = n_samples, # Number of rows / samples
               k = 9) # Number of columns or parameter variables

head(lhs_raw)


# Initialize DF
lhs_df <- setNames(data.frame(matrix(ncol = 17, nrow = n_samples)),
                   c("sim_id" , "save_every_n_steps", "verbose", "steps",    
                      "mortality_rate", "n_sp_init" , "n_alleles_init", "seed_disp_dist", 
                      "neighbor_radius", "seeds_per_adult", "range_cndd",
                      "mean_gndd", "range_gndd", "width", "height",            
                      "migration_rate",     "dispersal_mode" ))

head(lhs_df)


# Set variables that will not change across sensitivity analysis
  lhs_df$sim_id <- 1:nrow(lhs_df)
  lhs_df$verbose <- 0
  lhs_df$mortality_rate <- 0.1
  lhs_df$mean_gndd <- 0
  lhs_df$range_gndd <- 0
  lhs_df$dispersal_mode <- 0
  
  
# Set variables that will change
  
# Steps  
  lhs_df$steps <- round(qunif(lhs_raw[, 1], min = 5000, max = 20000))
  lhs_df$save_every_n_steps <- lhs_df$steps # So that it saves at last step
  
# N sp init
  lhs_df$n_sp_init <- round(qunif(lhs_raw[, 2], min = 25, max = 75))
  
# N alleles init
  lhs_df$n_alleles_init <- round(qunif(lhs_raw[, 3], min = 5, max = 20))
  
# Seed disp dist
  lhs_df$seed_disp_dist<- round(qunif(lhs_raw[, 4], min = 1, max = 15))
  
# Neighbor radius?
  lhs_df$neighbor_radius <- round(qunif(lhs_raw[, 5], min = 10, max = 30))
  
# Seeds per adult
  lhs_df$seeds_per_adult <- round(qunif(lhs_raw[, 6], min = 1000, max = 20000))
  
# Mean cndd - 
#  lhs_df$mean_cndd <- qunif(lhs_raw[, 7], min = -.2, max = .9)
  
# Range cndd
  lhs_df$range_cndd <- qunif(lhs_raw[, 7], min = .1, max = .4)
  
# Width / height
  lhs_df$width <- round(qunif(lhs_raw[, 8], min = 50, max = 200))
  lhs_df$height <- lhs_df$width
  
# Migration rate  
  lhs_df$migration_rate <- qunif(lhs_raw[, 9], min =  0.00001, max = 0.01)
  
  
  
  head(lhs_df)
  
  
# Set range of CNDD values for each sim  
  cndd_df = expand.grid(sim_id = lhs_df$sim_id,
                        mean_cndd = c(seq(-.2, .89, length.out = 10)))
  
# Join with lhs dataframe  
  lhs_df <- left_join(lhs_df, cndd_df)
  
  # Reset sim_id
  lhs_df$sim_id <- 1:nrow(lhs_df)
  
  
  # Make sure there are 10 counts of each
  lhs_df %>%
    group_by(steps, n_sp_init, n_alleles_init, seed_disp_dist, neighbor_radius,
             seeds_per_adult, range_cndd, width, height, migration_rate) %>%
    count() %>%
    ungroup() %>%
    select(n) # Should be 10 each
  
  
  # Make sure columns are in the right order
  lhs_df <- lhs_df %>%
    select(sim_id, save_every_n_steps, verbose, steps, mortality_rate,
           n_sp_init, n_alleles_init, seed_disp_dist, neighbor_radius, seeds_per_adult,
           mean_cndd, range_cndd, mean_gndd, range_gndd, width, height, migration_rate,
           dispersal_mode)
  

 # pairs(lhs_df)


## Write to file
write.csv(lhs_df, file = "./parameter files/run_lhs_002.csv", row.names = FALSE)








