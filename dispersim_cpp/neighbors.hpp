//
//  neighbors.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/6/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef neighbors_hpp
#define neighbors_hpp

#include <stdio.h>
#include <random>
#include <vector>
#include "params.hpp"
#include "sim.hpp"
#include <boost/random/binomial_distribution.hpp>



// Create a neighbors class that is re-used during demographic stage that holds
// Species ID of neighbors
// Genotypes of neighbors
// Gen 1d index
// Seeds by genotype
// RNG for dispersal for each neighbor that's called each iteration
// size will be dependent on R = or how many cells in radius are searched for potential dispersal

// Can have a printing function that prints out ## of seeds per species per genotype, etc


class Neighbors {
    
public:
    
    // Fields
    int n_neighbors;
    std::vector<int> nn_sp;
    std::vector<int> nn_index;
    std::vector<int> nn_gen;
    std::vector<int> nn_gen_1d_index;
    std::vector<bool> nn_gen_1d_index_dupe; // TRUE if is a duplicate
    std::vector<float> probabilities; // Recruitment probabilities
    
   // std::vector<std::binomial_distribution<int> > seed_rng;
    std::vector<boost::binomial_distribution<> > seed_rng;
   
    std::vector<float> seeds_by_sp; // Initialize to 0
    std::vector<float> dead_seeds_sp;
    std::vector<float> seeds_by_gen;
    float seeds_total;    

    //Constructor function
    Neighbors(Params& params);
    
    //Initialize seed RNGS based on distance
    void initSeedRNG(Params& params);
    
    // Find 1d indices of neighhbors
    void getNeighborIndex(int focal_cell, Params& params);
    
    // Updates nn_sp, nn_gen - calls getNeighborIndex
    void updateNeighbors(int focal_cell, Params& params, Sim& sim);
    
    // Generate numeber of seeds contributed by each neighbor
    // Based on distance and neg exponential dispersal kernel
    // Updates seeds_by_sp and seeds_by_gen
    void disperseSeeds(std::mt19937& generator);
    
    
    // Genotype dependent density dependence
    void GNDD(std::vector<float>& gndd_sp);
    
    // Conspecific NDD
    void CNDD(std::vector<float>& cndd_sp);
    
    // NDD simulaneous for GNDD and CNDD
    void NDD(std::vector<float>& gndd_sp, std::vector<float>& cndd_sp);
    
    // Count total number of seeds
    void totalSeeds();
    
    // Choose recruit based on relative frequency in cell after GNDD and CNDD
    void chooseWinner(int focal_cell, Params& params, Sim& sim);
    
    // Print status
    void printStatus(Params& params);
    
    // Reset - set values back to 0
    void reset();
    
};



#endif /* neighbors_hpp */
