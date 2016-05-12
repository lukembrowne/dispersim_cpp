//
//  utils.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/29/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "utils.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include "summary.hpp"


// Returns probability of dispersal for negative exponential, discretized
// Following Banitz et al. 2008
// alpha = 1 / mean dispersal distance
// dij = distance between source and endpoint
// R = Radius of cells from source cell
float neg_expo_discrete(float alpha, float dij, int R){
    
    float numerator = exp(-alpha * dij);
    
    float denominator{0.0};
    
    for(int i = -R; i < R; i++){
         for(int j = -R; j < R; j++){
             denominator += exp(-alpha * sqrt(pow(i, 2) + pow(j, 2)));
         }
    }
    
    float prob;
    prob = numerator / denominator;
    
    return(prob);
    
}

// Write summary stats to file
// Accepts a vector of summary objects and writes them to file

void write_summary(std::vector<Summary_step>& summary_over_time){
    
    // Start file stream
    std::ofstream out_summary("summary_over_time_overall.txt");
    
    // Set column names
    out_summary << "step \t sp_richness \t allelic_richness_avg \n";
    
    // Loop through summaries by step
    for(int i = 0; i < summary_over_time.size(); i++){
        out_summary << summary_over_time[i].step << "\t";
        out_summary << summary_over_time[i].sp_richness << "\t";
        out_summary << summary_over_time[i].allelic_richness_avg << "\t";
        out_summary << "\n"; // End line
    }
    
    out_summary.close();
    std::cout << "Successfully wrote summary_over_time_overall.txt to file... \n";

}



// Write landscape for both species and genotypes to file
// Export as tab delimited file

void write_landscape(std::vector<int>& sp,
                     std::vector<int>& gen,
                     int height, int width,
                     std::string species_filename){
    
    std::ofstream out_species(species_filename.c_str());
   // std::ofstream out_genotypes("landscape_genotypes.txt");
    
    for(int row = 0; row < height; row++){
        for(int col = 0; col < width; col++){
            out_species << sp[row*width + col] << "\t";
     //       out_genotypes<< gen[i*width + j] << "\t";

        } // End column loop
        
        out_species << "\n";
     //   out_genotypes << "\n";
        
    } // End row loop
    
    out_species.close();
   // out_genotypes.close();
    
    std::cout << "Successfully wrote " << species_filename << " to file... \n";
    
}

