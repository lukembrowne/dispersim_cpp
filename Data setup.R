## Setwd, load libraries, and load functions
setwd("~/Dropbox/dispersim_cpp")

library(tidyverse)
library(RColorBrewer) #to use brewer.pal
library(fields) #to use designer.colors
library(viridis)

source("./r scripts/functions.R")

## Set run name
run_name <- "run_030"

## Load in data
dat_sp = read_summary_out(path_to_folder = 
                            paste("./data/", run_name,"/summary_out/", sep = ""),
                          type = "by_sp")
dat_sp

## Merge with parameter values
params <- read_csv(paste("./data/", run_name, "/", run_name, ".csv", sep = ""))

dat_sp <- left_join(dat_sp, params, by = "sim_id")

glimpse(dat_sp)


## Rank species by abundance
## If abundance is a tie, will continute sequential counts of ranks
## Eg. 0, 0, 0, would be 97..98..99 ranks rather than all be a tie

## Loop through sim_ids, then through steps, then assign rank abundance

dat_sp$abundance_rank <- NA

for(id in unique(dat_sp$sim_id)){
  
  cat("On sim id:", id, "... \n")
  
  for(step in unique(dat_sp$step)){
    
    indices <- which(dat_sp$sim_id == id & dat_sp$step == step)
    
    ## Length should equal number of species
    if(length(indices) != dat_sp$n_sp_init[1]){
      cat("STOPPING \n")
      stop()
    }
    
    dat_sp$abundance_rank[indices[order(dat_sp$abundance[indices], decreasing = T)]] <- 1:dat_sp$n_sp_init[1]
    
  } # End step loop
  
} # End sim_id loop


## Subset down to last step, or not
dat_sp_allsteps <- dat_sp

dat_sp <- filter(dat_sp, step == max(step))
nrow(dat_sp)/dat_sp$n_sp_init[1] ## Number of simulations

# dat_sp %>%
#   filter(mean_cndd != 1) %>% # Filter out neutral models
# ggplot(aes(x = mean_cndd, y = abundance)) + geom_point() + facet_grid(.~sp + mean_gndd)


