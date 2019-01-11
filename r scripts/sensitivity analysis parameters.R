
### Writing simulation parameters to tab delimited text file so that a BASH 
## script can read this file and set up jobs on cluster

library(dplyr)


setwd("~/Dropbox/dispersim_cpp/r scripts")

reps = 5 ## How many reps per parameter set? 

### Grid of NDD values
# This one is important! Sets mean CNDD across simulations
ndd_df = expand.grid( mean_cndd = c(seq(.05, .85, length.out = 5), 1), # Go up to 1 for neutral dynamics
                      mean_gndd = c(seq(.05, .85, length.out = 5), 1)) 

ndd_df

### Add in parameters to vary

ndd_df_base <- ndd_df

## Baseline/default values

# Set 0
  ndd_df_base$steps <- 10000
## Set 1
  ndd_df_base$n_sp_init <- 50
  ndd_df_base$n_alleles_init <- 10
# Set 2  
  ndd_df_base$seed_disp_dist <- 5
# Set 3
  ndd_df_base$neighbor_radius <- 20
# Set 4
  ndd_df_base$range_cndd <- .2
  ndd_df_base$range_gndd <- 0.0
# Set 5
  ndd_df_base$seeds_per_adult <- 5000
# Set 6
  ndd_df_base$width <- 100
  ndd_df_base$height <- 100
# Set 7
  ndd_df_base$migration_rate <- 0.001

  
## Making master ndd_df  
  
## Set 0 - Steps
  ndd_df1 <- mutate(ndd_df_base, steps = 5000)
  ndd_df2 <- mutate(ndd_df_base, steps = 20000)
  
  ndd_df <- bind_rows(ndd_df1, ndd_df2)

# Set 1 - Sp and allele init
  ndd_df1 <- mutate(ndd_df_base, n_sp_init = 25, n_alleles_init = 5)
  ndd_df2 <- mutate(ndd_df_base, n_sp_init = 75, n_alleles_init = 20)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)
  
  
# Set 2 - Seed dispersal distance
  ndd_df1 <- mutate(ndd_df_base, seed_disp_dist = 2)
  ndd_df2 <- mutate(ndd_df_base, seed_disp_dist = 20)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)

# Set 3 - Neighborhood radius
  ndd_df1 <- mutate(ndd_df_base, neighbor_radius = 10)
  ndd_df2 <- mutate(ndd_df_base, neighbor_radius = 30)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)
  
# Set 4 - NDD range
  ndd_df1 <- mutate(ndd_df_base, range_cndd = .1, range_gndd = .1)
  ndd_df2 <- mutate(ndd_df_base, range_cndd = .4, range_gndd = .4)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)
  
# Set 5 - Seeds per adult
  ndd_df1 <- mutate(ndd_df_base, seeds_per_adult = 1000)
  ndd_df2 <- mutate(ndd_df_base, seeds_per_adult = 10000)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)
    
# Set 6 - Landscape size
  ndd_df1 <- mutate(ndd_df_base, height = 50, width = 50)
  ndd_df2 <- mutate(ndd_df_base,  height = 200, width = 200)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)

# Set 7 - Migration rate
  ndd_df1 <- mutate(ndd_df_base, migration_rate = 0.01)
  ndd_df2 <- mutate(ndd_df_base,  migration_rate = 0.00001)
  
  ndd_df <- bind_rows(ndd_df, ndd_df1, ndd_df2)
  
  
#### Non-varying columns
  ndd_df$save_every_n_steps = 2000
  ndd_df$verbose = 0
  ndd_df$mortality_rate = 0.1
  ndd_df$dispersal_mode = 0
    

out <- ndd_df

## Make sure that neutral communities have range 0

if(any(out$mean_cndd == 1)){
  out$range_cndd[out$mean_cndd == 1] <- 0
  out$range_gndd[out$mean_gndd == 1] <- 0
}


### Repeat each parameter set based on number of reps
out <- out[sort(rep(row.names(out), reps)), ]



## Add sim_id as first column
out <- cbind(sim_id = 1:nrow(out), out)

## Re-arrange column order
out <- out[, c("sim_id", "save_every_n_steps", "verbose",          
                "steps", "mortality_rate", "n_sp_init",         
               "n_alleles_init",  "seed_disp_dist", "neighbor_radius",   
                "seeds_per_adult", "mean_cndd",  "range_cndd",        
                "mean_gndd",  "range_gndd", "width" ,            
                "height",  "migration_rate","dispersal_mode")]


### Check to make sure cndd or gndd values don't go over 1 - positive density dependence

if(any((out$mean_cndd + (out$range_cndd/2)) > 1.0)){
  cat("WARNING!!!! Positive density dependence detected... \n
      Need to reduce either mean or range of CNDD values")
}

if(any((out$mean_gndd + (out$range_gndd/2)) > 1.0)){
  cat("WARNING!!!! Positive density dependence detected... \n
      Need to reduce either mean or range of GNDD values")
}


## Printing and plotting
cat("Total simulations:", nrow(out), "\n")

cat("Expected number of hours:", nrow(out)/96*1.25, "\n")



plot(jitter(out$mean_cndd), jitter(out$mean_gndd))


### Arrow plot of GNDD and CNDD values

plot(sort(unique(out$mean_cndd)), 1:length(unique(out$mean_cndd)))
arrows(x0 = sort(unique(out$mean_cndd)) - out$range_cndd[1]/2,
       x1 = sort(unique(out$mean_cndd)) + out$range_cndd[1]/2,
       y0 = 1:length(unique(out$mean_cndd)),
       code = 3, angle = 90, length = 0.1)






## Example of what matrix plot would look like
out %>%
  # filter(mean_cndd != 1) %>% # Filter out neutral models
  group_by(mean_cndd, mean_gndd) %>%
  summarise_all(mean) %>%
  ggplot(aes(x = mean_gndd, y = mean_cndd, fill = mean_gndd)) +
  geom_tile() + theme_bw()  + scale_fill_distiller(palette = "Spectral") + 
  geom_text(aes(x = mean_gndd, y = mean_cndd, label = round(mean_gndd,1)))



## Write to file
write.csv(out, file = "run_031.csv", row.names = FALSE)








