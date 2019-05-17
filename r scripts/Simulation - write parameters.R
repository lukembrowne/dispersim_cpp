
### Writing simulation parameters to tab delimited text file so that a BASH 
## script can read this file and set up jobs on cluster

reps = 2 ## How many reps per parameter set? 

### Grid of NDD values
# This one is important! Sets mean CNDD across simulations
ndd_df = expand.grid(#mean_cndd = c(seq(-.2, .9, length.out = 20), 1), # Go up to 1 for neutral dynamics
                    mean_cndd = c(seq(-.2, .89, length.out = 40), 1),
                     mean_gndd = c(1),
                     seed_disp_dist = c(5, 10, 20),
                     dispersal_mode = c(0, 1)) # Adding in global dispersal

# Remove simulations with varying dispersal distance from global dispersal
# Would need to adjust threshold for seed dispersal distance if that changes

if(any(ndd_df$dispersal_mode > 0)){
  ndd_df <- ndd_df[-which(ndd_df$dispersal_mode == 1 & ndd_df$seed_disp_dist > 
                            min(ndd_df$seed_disp_dist)), ]
  
  ## For global dispersal, change dispersal distance to -1 for plotting ease
  ### THIS NEEDS TO BE UPDATED
  ndd_df$seed_disp_dist[ndd_df$seed_disp_dist == min(ndd_df$seed_disp_dist) & 
                          ndd_df$dispersal_mode == 1] <- -1
}

ndd_df

out <- data.frame(
  
  ###1 - How often to save out data save_every_n_steps
  save_every_n_steps = 2000,
  
## 2 - Whether printing is verbose or not - do not run verbose on cluster
  verbose = 0,

## 3 -  Number of steps
## Each step = 5-10 years with 0.10 % mortality rate
  steps = 10000,

## 4 - Mortality rate
## Proportion of landscape dying per step
  mortality_rate = 0.1,

## 5 - Number of initial species
  n_sp_init = 50,

## 6 - Number of initial alleles per species
  n_alleles_init = 10,

## 7 - Seed dispersal distance
## In units of cells
  seed_disp_dist = ndd_df$seed_disp_dist,

## 8 - Neighbor radius for dispersal
## In units of cells
## Miranda = 5; Banitz = 20
  neighbor_radius = 20,

## 9 - Seeds per adult
## Equal to fecundity..
  seeds_per_adult = 10000,

## 10, 11 - CNDD parameters
## Values closer to 1 = weaker NDD
## Range is the max cndd value - min cndd value. 
 mean_cndd = ndd_df$mean_cndd,
 range_cndd = 0.20,

## 12, 13 - GNDD parameters
## Values closer to 1 = weaker NDD
  mean_gndd = ndd_df$mean_gndd,
  range_gndd = 0.00,

## 14, 15  - Landscape parameters - height, width
# In units of cells, also width x height equals total number of adults in simulation
  width = 100,
  height = 100,

## 16 - Migration rate
## Immigrant per recruit (~1 in 10,000) is from BCI paper; 1 in 9000 used in Muller Landau 2007
  migration_rate = 0.0001,

## 17 - Dispersal mode, 1 == global; 0 == local
  dispersal_mode = ndd_df$dispersal_mode

) # End data frame

out

## Make sure that neutral communities have range 0

if(any(out$mean_cndd == 1)){
  out$range_cndd[out$mean_cndd == 1] <- 0
  out$range_gndd[out$mean_gndd == 1] <- 0
}


### Repeat each parameter set based on number of reps
out <- out[sort(rep(row.names(out), reps)), ]



## Add sim_id as first column
out <- cbind(sim_id = 1:nrow(out), out)


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


### Arrow plot of GNDD and CNDD values
unique_cndd <- unique(out$mean_cndd)


plot(sort(unique_cndd), 1:length(unique_cndd), 
     xlim = c(-.25, 1.1),
     pch = 19, xlab = "Strength of NDD", bty = "l",
      ylab = "Community ID", yaxt = "n")
axis(side = 2, at = 1:10, labels = as.character(1:10), las = 1)
arrows(x0 = sort(unique_cndd) - out$range_cndd[1]/2,
       x1 = sort(unique_cndd) + out$range_cndd[1]/2,
       y0 = 1:length(unique_cndd),
       code = 3, angle = 90, length = 0.1)
abline(v = 1)


## Write to file
#write.csv(out, file = "./parameter files/run_036.csv", row.names = FALSE)








