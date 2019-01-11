

##### ANIMATION #######

library(animation)
library(RColorBrewer)
library(fields)
### Testing out animation

path_to_landscapes <- "./data/animation_ndd/"

file_list <- list.files(path = path_to_landscapes,
                        recursive = TRUE)

sim = "028"


landscape_files <- file_list[grep(paste("landscape_species_sim_", sim, "_step", sep = ""), file_list)]

landscape_files


# ## Order landscape files
landscape_files

steps <- unlist(lapply(strsplit(unlist(strsplit(landscape_files, ".txt")), "_"), function(x) x[6]))

steps




## Colors for plot
ColRamp <- rev(designer.colors(n=50, col=brewer.pal(11, "Spectral")))


## Loop through files and animate

saveGIF({

#oopt = ani.options(interval = 0.1, nmax = 50)
x = 1
for(landscape in landscape_files){

  
  print(landscape)
  
  landscape_species <- read.table(paste(path_to_landscapes, 
                                        landscape, sep = ""), header = FALSE)

  landscape_species <- as.matrix(landscape_species)

  image(landscape_species, zlim = c(0, 49), asp = 1, col = ColRamp,
        bty = "n", xaxt = "n", yaxt = "n", main = paste("Step:", steps[x]))

  ani.pause()
x = x + 1
}
#ani.options(oopt)

},
movie.name = "NDD.gif", interval = 0.1, nmax = 50, ani.width = 750, ani.height = 750
)
##############


