//
//  sim.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/12/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef sim_hpp
#define sim_hpp

#include <stdio.h>
#include "params.hpp"
#include <vector>
#include <random>
#include <boost/random/mersenne_twister.hpp>


// Sim class that holds information and parameters about simulation

class Sim {
    
public:
    
    // Constructor function
    Sim(Params& params, int n_dead_per_step);
    
    
    // Member variables
    std::vector<int> sp; // Holds species ID for each cell
    std::vector<int> gen; // Holds genotype for each cell
    
    std::vector<int> empty_cell_indices; // Holds indices of cells about to die
    
    std::vector<float> cndd_sp; // Holds CNDD values for each species
    std::vector<float> gndd_sp; // Holds GNDD values for each species
    
    // RNGS
    std::random_device device;
    std::mt19937 generator;
    boost::mt19937 boost_generator;
    std::uniform_int_distribution<int> cell_rng;
    std::uniform_real_distribution<> migration_rng;
    std::uniform_int_distribution<int> species_rng;
    std::uniform_int_distribution<int> gen_rng;

    
};




#endif /* sim_hpp */
