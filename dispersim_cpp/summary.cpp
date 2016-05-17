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
    
    // Species shannon diversity
    Summary_step::calc_sp_shannon(sim, params);
    
    // Allelic richness
    Summary_step::calc_allelic_richness_shannon(sim, params);
    
    // Average allelic richness
    allelic_richness_avg = Summary_step::calc_avg(allelic_richness_by_sp);
    
    allelic_shannon_avg = Summary_step::calc_avg(allelic_shannon_by_sp);
    
    // Print results
    Summary_step::print();
    
}


// SPECIES RICHNESS
void Summary_step::calc_sp_richness(std::vector<int>& sp){
    
    // Copy and sort species vector before eliminating duplicates
    std::vector<int> sp_copy(sp);
    
    std::sort(sp_copy.begin(), sp_copy.end());
    
    // Eliminate duplicates - size is now # of unique species
    sp_copy.erase(std::unique(sp_copy.begin(), sp_copy.end()),
                  sp_copy.end());
    
    // Assign to Summary field
    sp_richness = sp_copy.size();
    
}


// SHANNON DIVERSITY - SPECIES

void Summary_step::calc_sp_shannon(Sim& sim, Params& params){
    
    // Need to calculate relative abundance of each species
    
    std::vector<float> shannon_holder(params.n_sp_init, 0.0); // Temp vector
    
    float shannon{0.0};
    
    // Loop over species vector
    for(int i = 0; i < sim.sp.size(); i++){
        shannon_holder[sim.sp[i]] += 1.0;
    }
    
    // Divide by total to get relative abundance
    for(int i = 0; i < shannon_holder.size(); i++){
        shannon_holder[i] = shannon_holder[i] / params.area;
        
        // Can't take log of 0.. so skip if 0
        if(shannon_holder[i] == 0) continue;
        
        // Then multiple by natural log of that proportion
        shannon_holder[i] = shannon_holder[i] *
                                std::log(shannon_holder[i]);
        // Sum
        shannon += shannon_holder[i];
    }
    
    // Multiply by -1 and assign to class object
    sp_shannon = -shannon;
    
}



// ALLELIC RICHNESS
void Summary_step::calc_allelic_richness_shannon(Sim& sim,
                                                 Params& params){
    
    std::vector<float> allelic_richness(params.n_sp_init, 0.0);
    
    std::vector<float> shannon_holder(params.n_alleles_init, 0.0); // Temp vector
    
    float shannon{0.0};

    // Loop through species
    for(int i = 0; i < params.n_sp_init; i ++){
        
        // Temp vector to hold alleles for each species
        std::vector<int> allele_temp;
        
        // Loop through species vector
        for(int j = 0; j < sim.sp.size(); j++){
            
            // Probably a more efficient way of doing this??
            if(i == sim.sp[j]){ // If species match
                allele_temp.push_back(sim.gen[j]);
            } // End if
        }
        
        // Calculating both richness and shannon in same function
        // To save on computation of getting alleles for each species
        
        /////
        // SHANNON DIVERSITY
        
        // Loop over allele temp vector
        for(int i = 0; i < allele_temp.size(); i++){
            shannon_holder[allele_temp[i]] += 1.0;
        }
        
        // Divide by total to get relative abundance
        for(int i = 0; i < shannon_holder.size(); i++){
            shannon_holder[i] = shannon_holder[i] / allele_temp.size();
            
            // Can't take log of 0.. so skip if 0
            if(shannon_holder[i] == 0) continue;
            
            // Then multiple by natural log of that proportion
            shannon_holder[i] = shannon_holder[i] *
                                std::log(shannon_holder[i]);
            // Sum
            shannon += shannon_holder[i];
        }
        
        // Multiply by -1 and assign to class object
        allelic_shannon_by_sp.push_back(-shannon);

        // Reset vars
        shannon = 0;
        std::fill(shannon_holder.begin(),
                  shannon_holder.end(), 0.0);
        
        
        /////
        // RICHNESS
        
        // Sort vector
        std::sort(allele_temp.begin(), allele_temp.end());
        
        // Eliminate duplicates
        allele_temp.erase(std::unique(allele_temp.begin(), allele_temp.end()), allele_temp.end());
        
        // Size == number of alleles
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


// CALCULATE AVERAGE OF A METRIC
// Used for averaging allelic richness across species, etc
float Summary_step::calc_avg(std::vector<float> metric){
    
    float sum = 0.0;
    
    for(int i = 0; i < metric.size(); i++){
        sum += metric[i];
    }
    
    float avg = sum / metric.size();
    
    return(avg);
}


// PRINTING FUNCTION
void Summary_step::print(){
    std::cout << "Step: " << step <<
    " | Sp. richness: " << sp_richness <<
    " | Sp. shannon: " << sp_shannon <<
    " | Allelic richness (avg.): " << allelic_richness_avg <<
    " | Allelic shannon (avg.): " << allelic_shannon_avg <<

    "\n";
    
}




