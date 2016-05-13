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
    
    // Genetic diversity
    std::vector<int> allelic_richness_by_sp;
    float allelic_richness_avg;
    
    //Constructor function
    Summary_step(Sim& sim, Params& params, int step);
    
    // Member functions
    void calc_sp_richness(std::vector<int>& sp);
    void calc_allelic_richness(int n_sp_init, std::vector<int>& sp,
                               std::vector<int>& gen);
    float calc_avg(std::vector<int> metric); // Maybe need to overload another for float?
    void print(); // Print to console
};




#endif /* summary_hpp */
