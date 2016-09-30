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
#include <sstream>
#include <iomanip>


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

void write_summary(std::vector<Summary_step>& summary_over_time, Params& params){
    
    // Start file stream
    
    // Set file names
    // Write landscape of species to tab delimited .txt file
    std::string summary_overall_filename = "./summary_out/summary_overall_";
    std::string summary_overall_sp_filename = "./summary_out/summary_by_sp_";
    std::string suffix = ".txt";
    
    // Write to buffer to add leading 0s
    std::stringstream overall_buffer;
    std::stringstream overall_by_sp_buffer;
    
    // Adding sim_id and step number to file name, setw and setfill add leading for sorting
    overall_buffer << summary_overall_filename << "sim_" << std::setw(3) << std::setfill('0') << params.sim_id << suffix;
    overall_by_sp_buffer << summary_overall_sp_filename << "sim_" << std::setw(3) << std::setfill('0') << params.sim_id << suffix;

    
    std::ofstream out_summary_overall(overall_buffer.str().c_str());
    std::ofstream out_summary_by_sp(overall_by_sp_buffer.str().c_str());
    
    // Set column names
    out_summary_overall << "sim_id \t step \t sp_richness \t sp_shannon \t allelic_richness_avg \t allelic_shannon_avg\n";
    
    out_summary_by_sp << "sim_id \t step \t sp_richness \t sp_shannon \t sp \t abundance \t allelic_richness \t allelic_shannon\n";

    
    // Loop through summaries by step
    for(int i = 0; i < summary_over_time.size(); i++){
        
        
        // Overall
        out_summary_overall << params.sim_id << "\t";
        out_summary_overall << summary_over_time[i].step << "\t";
        out_summary_overall << summary_over_time[i].sp_richness << "\t";
        out_summary_overall << summary_over_time[i].sp_shannon << "\t";
        out_summary_overall << summary_over_time[i].allelic_richness_avg << "\t";
        out_summary_overall << summary_over_time[i].allelic_shannon_avg << "\t";
        out_summary_overall << "\n"; // End line
        
        // By species
        
        // Loop over species
        for(int j = 0; j < summary_over_time[i].allelic_richness_by_sp.size(); j++){
            
            out_summary_by_sp << params.sim_id << "\t";
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


// Write parameter values to file
// Export as tab delimited file
// 2 columns - one with parameter name, one with value

void write_params(std::string params_filename, Params& params){
    
    std::ofstream out_params(params_filename.c_str());

    // Set up column wrow
    out_params << "param \t value \n";
    
    out_params << "sim_id \t" << params.sim_id << "\n";
    out_params << "steps \t" << params.steps << "\n";
    out_params << "mortality_rate \t" << params.mortality_rate << "\n";
    out_params << "n_sp_init \t" << params.n_sp_init << "\n";
    out_params << "n_alleles_init \t" << params.n_alleles_init << "\n";
    out_params << "seed_disp_dist \t" << params.seed_disp_dist << "\n";
    out_params << "neighbor_radius \t" << params.neighbor_radius << "\n";
    out_params << "seeds_per_adult \t" << params.seeds_per_adult<< "\n";
    
    out_params << "mean_cndd \t" << params.mean_cndd << "\n";
    out_params << "range_cndd \t" << params.range_cndd << "\n";
    out_params << "mean_gndd \t" << params.mean_gndd << "\n";
    out_params << "range_gndd \t" << params.range_gndd << "\n";
    
    out_params << "width \t" << params.width << "\n";
    out_params << "height \t" << params.height << "\n";
    out_params << "area \t" << params.area << "\n";
    
    out_params << "migration_rate \t" << params.migration_rate << "\n";
    
    out_params << "dispersal_mode \t" << params.dispersal_mode << "\n";
    
    // Close out file
    out_params.close();

    
    std::cout << "Successfully wrote " << params_filename << " to file... \n";
    
}

