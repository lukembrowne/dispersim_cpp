//
//  summary.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/6/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "summary.hpp"
#include <vector>
#include <iostream>





// CONSTRUCTOR FUNCTION

// Calculate summary statistics when object is created

Summary_step::Summary_step(Sim& sim, Params& params, int step){
    
    // Save step
    Summary_step::step = step;
    
    // Species richness
    Summary_step::calc_sp_richness(sim.sp);
    
    // Allelic richness
    Summary_step::calc_allelic_richness(params.n_sp_init, sim.sp, sim.gen);
    
    // Average allelic richness
    allelic_richness_avg = Summary_step::calc_avg(allelic_richness_by_sp);
    
    // Print results
    Summary_step::print();
    
}


// SPECIES RICHNESS
void Summary_step::calc_sp_richness(std::vector<int>& sp){
    
    // Copy and sort species vector before eliminating duplicates
    std::vector<int> sp_copy(sp);
    
    std::sort(sp_copy.begin(), sp_copy.end());
    
    // Eliminate duplicates - size is now # of unique species
    sp_copy.erase(std::unique(sp_copy.begin(), sp_copy.end()), sp_copy.end());
    
    // Assign to Summary field
    sp_richness = sp_copy.size();

    
}


// ALLELIC RICHNESS
void Summary_step::calc_allelic_richness(int n_sp_init, std::vector<int>& sp, std::vector<int>& gen){
    
    std::vector<float> allelic_richness(n_sp_init, 0.0);

    // Loop through species
    for(int i = 0; i < n_sp_init; i ++){
        
        // Temp vector to hold alleles for each species
        std::vector<int> allele_temp;
        
        for(int j = 0; j < sp.size(); j++){ // Loop through species vetcor
            
            if(i == sp[j]){ // If species match
                allele_temp.push_back(gen[j]);
            } // End if
        }
        
        // Sort vector
        std::sort(allele_temp.begin(), allele_temp.end());
        
        // Eliminate duplicates
        allele_temp.erase(std::unique(allele_temp.begin(), allele_temp.end()), allele_temp.end());
        
        // Size == number of alleles
        // Probably a more efficient way of doing this??e
        allelic_richness_by_sp.push_back(allele_temp.size());
        
        } // End species loop

}

// CALCULATE AVERAGE OF A METRIC
// Used for averaging allelic richness across species, etc
float Summary_step::calc_avg(std::vector<int> metric){
    
    float sum = 0.0;
    
    for(int i = 0; i < metric.size(); i++){
        sum += metric[i];
    }
    
    float avg = sum / metric.size();
    
    return(avg);
}


// PRINTING FUNCTION
void Summary_step::print(){
    
    std::cout << "Step: " << step << " | Sp. richness: " << sp_richness <<
    " | Allelic richness (avg.): " << allelic_richness_avg << "\n";
    
}




