//
//  summary.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/6/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include <stdio.h>
#include <vector>
#include "sim.hpp"
#include "params.hpp"



// Summary class object that calculates summary statistics of the simulation at a given step
class Summary_step {
    
public:
    
    // Fields
    int step;
    
    // Species diversity
    int sp_richness;
    float sp_shannon;
    
    // Genetic diversity
    std::vector<int> allelic_richness_by_sp; // Richness
    float allelic_richness_avg;
    std::vector<float> allelic_shannon_by_sp; // Shannon diversity
    float allelic_shannon_avg;
    
    //Constructor function
    Summary_step(Sim& sim, Params& params, int step);
    
    // Member functions
    void calc_sp_richness(std::vector<int>& sp);
    void calc_sp_shannon(Sim& sim, Params& params);
    void calc_allelic_richness_shannon(Sim& sim, Params& params);
    float calc_avg(std::vector<int> metric);
    float calc_avg(std::vector<float> metric); // Overloaded for float
    void print(); // Print to console
};




#endif /* summary_hpp */
