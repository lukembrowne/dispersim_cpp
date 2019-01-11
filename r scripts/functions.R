#### Function to read in overall summary files


## Reads in all the text files and consolidates into one file


read_summary_out <- function(path_to_folder, type){
  
  # Type should == either "overall" or "by_sp"
  
  file_list <- list.files(path_to_folder, full.names = T)
  
  ## Only overall summaries
  file_list_overall <- file_list[grep(type, file_list)]
  
  x = 1
  
  for(file in file_list_overall){
    
    temp_df <- read_tsv(file)
    
    ## If on first file, set up output DF, or add it to output df
    if(x == 1){
      out <- temp_df
    } else {
      out <- bind_rows(out, temp_df)
    }
    
    x = x + 1
    
  } # End file loop
  
  
  return(out)
  
}











# 
# 
# 
# 
# 
# library(RColorBrewer) #to use brewer.pal
# library(fields) #to use designer.colors
# 
# setwd("/Users/lukebrowne/Library/Developer/Xcode/DerivedData/dispersim_cpp-dmesgaqjucaepucwfndyaqnypvka/Build/Products/Debug")
# 
# 
# 
# ########
# ## Read in data - no header
# 
# landscape_species <- read.table("./landscape_out/landscape_species_FINAL.txt", header = FALSE)
# 
# landscape_species <- as.matrix(landscape_species)
# 
# 
# # Plot relative abundance
# 
# # barplot(sort(table(landscape_species), decreasing = TRUE), las = T, col = "skyblue2")
# # 
# # sort(table(landscape_species), decreasing = TRUE)
# # 
# 
# 
# ## Plot landscape
# ColRamp <- rev(designer.colors(n=50, col=brewer.pal(11, "Spectral")))
# 
# 
# image(landscape_species, asp = 1, col = ColRamp, bty = "n", xaxt = "n", yaxt = "n")
# 
# 
# 
# 













