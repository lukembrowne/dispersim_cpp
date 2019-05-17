## Setwd, load libraries, and load functions
  library(tidyverse)
  library(RColorBrewer) #to use brewer.pal
  library(fields) #to use designer.colors
  library(viridis)
  
  

# Read in simulation output -----------------------------------------------

  
## Set run name
  run_name <- "run_036"
  run_name <- "run_lhs_001"
  
  
  path_to_summaries_by_sp <- grep(pattern = "by_sp",
                            x = list.files(paste0("./data/", 
                                                  run_name, "/summary_out"),
                                           full.names = TRUE),
                            value = TRUE)
                            
  dat_sp <- plyr::ldply(path_to_summaries_by_sp, read_tsv)
  
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
  
  table(dat_sp$sim_id)
  
  # Which sim ids are missing?
  which(!1:max(dat_sp$sim_id) %in% dat_sp$sim_id)
  
  
  
  
#######################################  
### Plot patterns averaged across species

## Sp shannon
dat_sp %>%
  group_by(sim_id) %>% # Average across species in a simulation
  summarise_all(mean) %>%
  ggplot(aes(x = mean_cndd, y = sp_shannon, 
             fill = factor(seed_disp_dist),
             group = factor(seed_disp_dist))) +
  geom_point(pch = 21) + theme_bw() +
#    geom_smooth(aes(col = factor(seed_disp_dist)), se = FALSE) +
 # facet_wrap( ~ seed_disp_dist)
 # scale_fill_distiller(palette = "Spectral")
    NULL
  
  
## Allelic shannon
  dat_sp %>%
    group_by(sim_id) %>% # Average across species in a simulation
    summarise_all(mean) %>%
    ggplot(aes(x = mean_cndd, y = allelic_shannon, fill = factor(seed_disp_dist))) +
    geom_point(pch = 21) + theme_bw() +
    #   facet_wrap( ~ seed_disp_dist)
   # facet_wrap( ~ mean_gndd) +
  #  scale_fill_distiller(palette = "Spectral") +
    NULL
  
  
  
  
### Plotting species vs genetic diversity averaged across community
  
  ### Averaging across CNDD - This plot doesn't make sense to me anymore...
  # dat_sp %>%
  #   filter(abundance > 0) %>% # only species with positive abundance 
  #   #  filter(mean_gndd != 1) %>%
  #   #  mutate(neutral = ifelse(mean_gndd == 1, 1, 0)) %>%
  #   group_by(sim_id) %>%
  #   summarise_all(mean) %>%
  #   group_by(mean_cndd) %>%
  #   mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>%
  #   group_by(mean_cndd) %>%
  #   summarise_all(mean) %>%
  #   ggplot() +
  #   geom_point(aes(x = mean_cndd, y = sgdc), size = I(5), pch = 19) + 
  #   geom_hline(aes(yintercept = 0)) + theme_bw() 
  # 
  # 
  
  
## Scatterplot of species vs. genetic diversity at community level
  dat_sp %>%
    filter(abundance > 0) %>% # only species with positive abundane
     filter(mean_cndd != 1) %>% # Filter out neutral models
    group_by(sim_id) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = sp_shannon, y = allelic_shannon,
               fill = factor(seed_disp_dist))) + 
    geom_point(aes(col = factor(seed_disp_dist))) + theme_bw() +
    NULL
  
  
  
## Just neutral models
  dat_sp %>%
    filter(abundance > 0) %>% # only species with positive abundane
    filter(mean_cndd == 1) %>% # Only neutral models
    group_by(sim_id) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = sp_shannon, y = allelic_shannon)) + geom_point() + theme_bw() +
    facet_wrap( ~ mean_gndd ) + geom_smooth()  
  
  
  
  ####################
  ####################
  
  
  ### Calculate SGDC and plot based on rank abundance
  
  ## Grouped by gndd
  dat_sp %>%
    # filter(abundance > 0) %>% # only species with positive abundance
    # filter(mean_cndd != 1) %>% # Filter out neutral models
    #filter(mean_gndd != 1) %>% 
    group_by(abundance_rank, mean_gndd, seed_disp_dist) %>%
    mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>%
    ggplot(aes(x = abundance_rank, y = sgdc, fill = sgdc)) + 
    geom_point(pch = 21) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
    geom_hline(aes(yintercept = 0)) + facet_grid(. ~ seed_disp_dist)
  
  
  ## Based on SPECIES ID
  dat_sp %>%
    filter(abundance > 0) %>% # only species with positive abundance
    filter(mean_cndd != 1) %>% # Filter out neutral model
    group_by(sp, mean_gndd) %>%
    mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>%
    ggplot(aes(x = sp, y = sgdc, fill = sgdc)) + 
    geom_point(pch = 21) + theme_bw() + scale_fill_distiller(palette = "Spectral") +
    geom_hline(aes(yintercept = 0)) + facet_grid(. ~ seed_disp_dist)
  
  
  ### Histogram of SGDC values by species
  dat_sp %>%
    filter(abundance > 0) %>% # only species with positive abundance
    # filter(mean_cndd != 1) %>% # Filter out neutral models
    # filter(mean_gndd != 1) %>% 
    group_by(abundance_rank, mean_gndd) %>%
    mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>%
    ggplot(aes(x = sgdc)) + geom_histogram() + theme_bw()
  
  
  
  ## Based on SPECIES ID
  dat_sp %>%
    filter(abundance > 0) %>% # only species with positive abundance
    filter(mean_cndd == 1) %>% # Just neutral models
    group_by(sp, seed_disp_dist) %>%
    mutate( sgdc = cor(sp_shannon, allelic_shannon)) %>%
    ggplot(aes(x = sp, y = sgdc, col = abundance, size = abundance)) + 
    geom_point() + theme_bw() + 
    geom_hline(aes(yintercept = 0)) + 
    facet_grid(. ~ seed_disp_dist)
  
  dat_sp %>%
    filter(abundance > 0) %>%
    filter(mean_cndd != 1) %>%
    group_by(sp, seed_disp_dist) %>%
    filter(sp == 0 | sp == 49 | sp == 25) %>%
    ggplot(aes(x = sp_shannon, y = allelic_shannon)) +
    geom_point() + geom_smooth() +
    facet_wrap( ~ sp + seed_disp_dist) + theme_bw()
  
  
  dat_sp %>%
    filter(abundance > 0) %>%
    filter(mean_cndd != 1) %>%
    group_by(sp, seed_disp_dist) %>%
    filter(sp <= 2) %>%
    ggplot(aes(x = sp_shannon, y = allelic_shannon)) +
    geom_point() + geom_smooth() +
    facet_wrap( ~ sp + seed_disp_dist) + theme_bw()
  
  
  ## Of just specific
  dat_sp %>%
    filter(abundance > 0) %>%
    filter(mean_cndd != 1) %>%
    # filter(seed_disp_dist == 100) %>%
    filter(sp == 1) %>%
    group_by(sp) %>%
    ggplot(aes(x = sp_shannon, y = allelic_shannon, col = abundance, size = abundance)) +
    geom_point() + geom_smooth() +
    facet_wrap( ~ sp + seed_disp_dist) + theme_bw()
  
  
  

# Plot changes in values over time ----------------------------------------

  ####################
  ## Plotting changes in values over time
  
  # Species shannon
  # By migration rate
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(migration_rate, step) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = sp_shannon, color = factor(migration_rate))) + geom_line()
  
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(migration_rate, step) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = sp_richness, color = factor(migration_rate))) + geom_line()
  
  # By dispersal distance
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(seed_disp_dist, step) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = sp_shannon, color = factor(seed_disp_dist))) + geom_line()
  
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(seed_disp_dist, step, mean_cndd) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = sp_shannon, color = factor(seed_disp_dist))) + geom_line() +
    facet_grid(. ~ mean_cndd)
  
  
  ## Allelic shannon
  
  # By migration rate
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(migration_rate, step) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = allelic_shannon, color = factor(migration_rate))) + 
    geom_line()
  
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(migration_rate, step) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = allelic_richness, color = factor(migration_rate))) + 
    geom_line()
  
  # By dispersal distance
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(seed_disp_dist, step) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = allelic_shannon, color = factor(seed_disp_dist))) + geom_line()
  
  dat_sp_allsteps %>%
    filter(abundance > 0) %>%
    group_by(seed_disp_dist, step, mean_cndd) %>%
    summarise_all(mean) %>%
    ggplot(aes(x = step, y = allelic_shannon, color = factor(seed_disp_dist))) + geom_line() +
    facet_grid(. ~ mean_cndd)

  
  
  
  
  

# Plot landscape ----------------------------------------------------------

  sim_id <- "099"
  
  sp_id <- "000"
  
  landscape_species <- read.table(paste("./data/", run_name,"/landscape_out/landscape_species_sim_", sim_id, "_FINAL.txt", sep = ""), header = FALSE)
  
  ## For species
  # landscape_species <- read.table(paste("./data/", 
  #                                       run_name,"/landscape_out/landscape_genotype_sim_", sim_id,
  #                                       "_sp_", sp_id, "_FINAL.txt", sep = ""), header = FALSE)
  
  
  landscape_species <- as.matrix(landscape_species)
  
  
  # # Plot relative abundance
  # 
  barplot(sort(table(landscape_species), decreasing = TRUE), las = T, col = "skyblue2")
  
  sort(table(landscape_species), decreasing = TRUE)
  
  ## Plot landscape
  ColRamp <- rev(designer.colors(n=50, col=brewer.pal(11, "Spectral")))
  
  image(landscape_species, asp = 1, col = ColRamp, bty = "n", xaxt = "n", yaxt = "n")
  
  
  
  
  

# Plot dispersal kernel ---------------------------------------------------

  ## Imaginging a plot of squares with center point as focal point and then each square is colored as the relative probability of dispersal. Faceted plot showing the difference in dispersal probability based on neighborhood radius and average seed dispersal distance
  
  ## Or a 3d graph
  
  
  # alpha = 1 / mean dispersal distance
  # dij = distance between source and endpoint
  # R = Radius of cells from source cell
  
  
  
  neg_expo_discrete <- function(alpha, dij, R){
    
    numerator = exp(-alpha * dij)
    
    denominator <- 0
    
    for(i in -R:R){
      for(j in -R:R){
        #   cat(i,j, "\n")
        denominator = denominator +  exp(-alpha * sqrt(i^2 + j^2));
      }
    }
    
    probability <- numerator / denominator
    
    return(probability)
  }
  
  
  ## Fill up matrix of probabilities
  
  R = 20
  
  seed_disp_dist = c(5)
  
  
  probs <- matrix(NA, nrow = (2*R + 1), ncol = (2*R + 1))
  probs_df <- data.frame(x = NA, y = NA, prob = NA, seed_disp_dist)
  
  probs
  
  center <- ceiling(nrow(probs)/2)
  
  x = 1
  
  for(disp_dist in seed_disp_dist){
    
    for(i in 1:nrow(probs)){
      for(j in 1:ncol(probs)){
        #cat(i,j, "\n")
        
        ## Distance to center of matrix 
        dij <- sqrt((i - center)^2 + (j - center)^2)
        
        # print(dij)
        
        prob <- neg_expo_discrete(alpha = 1/disp_dist, dij = dij, R = R)
        
        # probs[i, j] <- prob
        
        probs_df[x, "x"] <- i
        probs_df[x, "y"] <- j
        probs_df[x, "prob"] <- prob
        probs_df[x, "seed_disp_dist"] <- disp_dist
        x = x + 1
      }
    }
    
    # Set center as 0
    probs_df[which(probs_df$x == center & probs_df$y == center), "prob" ] <- NA
    
  }
  
  
  ## 3d plot
  
  par(mfrow = c(length(seed_disp_dist), 2),
      mar = c(1.5, 1.5, 1.5, 1.5))
  
  ymax = 0.01
  yseedsmax = 40
  fecundity = 5000
  
  for(disp_dist in seed_disp_dist){
    
    probs_df_sub <- subset(probs_df, probs_df$seed_disp_dist == disp_dist)
    
    probs_df_sub <- probs_df_sub[order(probs_df_sub$x, probs_df_sub$y),]
    
    # Transform to matrix
    
    probs_mat <- matrix(NA, ncol = max(probs_df_sub$x), nrow = max(probs_df_sub$x))
    seeds_mat <- matrix(NA, ncol = max(probs_df_sub$x), nrow = max(probs_df_sub$x))
    
    for(row in 1:nrow(probs_df_sub)){
      probs_mat[probs_df_sub$x[row], probs_df_sub$y[row]] <- probs_df_sub$prob[row]
      seeds_mat[probs_df_sub$x[row], probs_df_sub$y[row]] <- rbinom(1, fecundity, 
                                                                    probs_df_sub$prob[row])
    }
    
    cat("Probability sum:", sum(probs_mat, na.rm = TRUE), "... \n")
    
    summary(c(probs_mat))
    summary(c(seeds_mat))
    
    ## Probability
    ## From example
    nrz <- nrow(probs_mat)
    ncz <- ncol(probs_mat)
    
    nbcol <- 100
    color <- viridis(nbcol, begin = 0.5)
    zfacet <- probs_mat[-1, -1] + probs_mat[-1, -ncz] + probs_mat[-nrz, -1] + probs_mat[-nrz, -ncz]
    facetcol <- cut(zfacet, nbcol)
    
    persp(z = probs_mat, zlim = c(0, ymax), theta = 90, phi = 20, expand = .5, 
          col = color[facetcol], 
          ltheta = 120, shade = 0.05, ticktype = "simple",
          xlab = "", ylab = "", zlab = "Dispersal \n probability",
          main = disp_dist, box = FALSE)
    
    ## Number of seeds
    ## From example
    nrz <- nrow(seeds_mat)
    ncz <- ncol(seeds_mat)
    
    nbcol <- 100
    color <- viridis(nbcol, begin = 0.5)
    zfacet <- seeds_mat[-1, -1] + seeds_mat[-1, -ncz] + seeds_mat[-nrz, -1] + seeds_mat[-nrz, -ncz]
    facetcol <- cut(zfacet, nbcol)
    
    persp(z = seeds_mat, theta = 90, phi = 20, expand = .25, 
          zlim = c(0, yseedsmax),
          col = color[facetcol], 
          ltheta = 120, shade = 0.05, ticktype = "detailed",
          xlab = "", ylab = "", zlab = "Seeds",
          main = paste(disp_dist, " | Sum:", sum(seeds_mat, na.rm = TRUE)))
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

