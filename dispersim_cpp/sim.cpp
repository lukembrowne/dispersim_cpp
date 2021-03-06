//
//  sim.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/12/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "sim.hpp"
#include "params.hpp"
#include <iostream>
#include <random>
#include <boost/random/mersenne_twister.hpp>




Sim::Sim(Params& params, int n_dead_per_step) :

// INIT RNGS

// Initialize Random number generators
// Good info on generating random numbers in C++
// http://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c

generator(device()),

// RNG that randomly chooses cell in landscape
cell_rng(0, (params.area-1)),

// Migration RNG - chooses between 0 and 1
migration_rng(0, 1),

// RNG for species - ranging from 0 to # of species - 1 (species ID used for indexing)
species_rng(0, (params.n_sp_init - 1)),

// Initialize genotype RNG
gen_rng(0, (params.n_alleles_init-1))

{
    
    // Stores index for cells that die each step
    empty_cell_indices.resize(n_dead_per_step);
    
    
    // Initialize species vector and place species randomly throughout landscape
    sp.resize(params.area); // Init vector that will hold species
    
        for (auto& iter : sp){
            iter = species_rng(generator);
            assert(iter <= params.n_sp_init);
        }
    
    
    // Initialize genotype vector and place genotypes randomly throughout landscape
    gen.resize(params.area); // Init vector that will hold species
    
        for (auto& iter : gen){
            iter = gen_rng(generator);
        }
    
    /////////////////////////////////
    // Assign CNDD values to species - conspecific density dependence
    // Lower species id == Stronger NDD
    cndd_sp.resize(params.n_sp_init);
    
        // Will divide by 0 here if only one species
        double cndd_increment = params.range_cndd/((double)params.n_sp_init - 1.0);
        int i = 0;
    
        // Figure out max and min cndd values
    double max_cndd = params.mean_cndd + (params.range_cndd/2);
    double min_cndd = params.mean_cndd - (params.range_cndd/2);
    
    
    for(double cndd_val = max_cndd; i <= params.n_sp_init; cndd_val -= cndd_increment){
        
        if(cndd_val >= 1.0){
            cndd_sp[i] = 1.0;
        } else {
            cndd_sp[i] = cndd_val;
        }
        
        i++;
        if(i == params.n_sp_init) break; // Make sure it doensn't go out of bounds
    }
    
    
    /////////////////////////////////
    // GENOTYPE DEPENDENT NDD (GNDD)
    // Lower species id == Stronger NDD
    gndd_sp.resize(params.n_sp_init);
    
        // Potential bug Will divide by 0 here if only one species
        double gndd_increment = params.range_gndd/((double)params.n_sp_init - 1.0);
        i = 0;
        
    
    // Figure out max and min cndd values
    double max_gndd = params.mean_gndd + (params.range_gndd/2);
    double min_gndd = params.mean_gndd - (params.range_gndd/2);

    
    for(double gndd_val = max_gndd; i <= params.n_sp_init; gndd_val -= gndd_increment){
        
        if(gndd_val >= 1.0){
            gndd_sp[i] = 1.0;
        } else {
            gndd_sp[i] = gndd_val;
        }
        i++;
        if(i == params.n_sp_init) break; // Make sure it doensn't go out of bounds
    }


    
    
    /////////////////////////////////
    // Print out information about simulation at initialization
    std::cout << "| ---------------------------- | \n";
    std::cout << "Steps: " << params.steps << " or ~ " << params.steps * 5 <<
    " years .. " << params.steps / 10 << " generations \n";
    std::cout << "Area: " << (params.width * 10 * params.height * 10)/10000 << " ha \n";
    std::cout << "Individuals: " << params.area << " .. " << params.area/params.n_sp_init << " per species \n";
    
}



