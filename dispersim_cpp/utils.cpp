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
    std::ofstream out_summary_overall("./summary_out/summary_overall.txt");
    std::ofstream out_summary_by_sp("./summary_out/summary_by_sp.txt");
    
    // Set column names
    out_summary_overall << "step \t sp_richness \t sp_shannon \t allelic_richness_avg \t allelic_shannon_avg\n";
    
    out_summary_by_sp << "step \t sp_richness \t sp_shannon \t sp \t abundance \t allelic_richness \t allelic_shannon\n";

    
    // Loop through summaries by step
    for(int i = 0; i < summary_over_time.size(); i++){
        
        
        // Overall
        out_summary_overall << summary_over_time[i].step << "\t";
        out_summary_overall << summary_over_time[i].sp_richness << "\t";
        out_summary_overall << summary_over_time[i].sp_shannon << "\t";
        out_summary_overall << summary_over_time[i].allelic_richness_avg << "\t";
        out_summary_overall << summary_over_time[i].allelic_shannon_avg << "\t";
        out_summary_overall << "\n"; // End line
        
        // By species
        
        // Loop over species
        for(int j = 0; j < summary_over_time[i].allelic_richness_by_sp.size(); j++){
            
            out_summary_by_sp << summary_over_time[i].step << "\t";
            out_summary_by_sp << summary_over_time[i].sp_richness << "\t";
            out_summary_by_sp << summary_over_time[i].sp_shannon << "\t";
            out_summary_by_sp << j << "\t"; // Species ID
            out_summary_by_sp << summary_over_time[i].abundance_by_sp[j] << "\t";
            out_summary_by_sp << summary_over_time[i].allelic_richness_by_sp[j] << "\t";
            out_summary_by_sp << summary_over_time[i].allelic_shannon_by_sp[j] << "\t";

            out_summary_by_sp << "\n"; // End line


        } // End species loop
    } // End summary step loop
    
    // Close files
    out_summary_overall.close();
    out_summary_by_sp.close();
    
    std::cout << "Successfully wrote summaries to file... \n";

}



// Write landscape for both species and genotypes to file
// Export as tab delimited file

void write_landscape(std::string species_filename, Sim& sim, Params& params){
    
    std::ofstream out_species(species_filename.c_str());
   // std::ofstream out_genotypes("landscape_genotypes.txt");
    
    for(int row = 0; row < params.height; row++){
        for(int col = 0; col < params.width; col++){
            out_species << sim.sp[row * params.width + col] << "\t";
     //       out_genotypes<< gen[i*width + j] << "\t";

        } // End column loop
        
        out_species << "\n";
     //   out_genotypes << "\n";
        
    } // End row loop
    
    out_species.close();
   // out_genotypes.close();
    
    std::cout << "Successfully wrote " << species_filename << " to file... \n";
    
}

