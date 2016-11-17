//
//  neighbors.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/6/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "neighbors.hpp"
#include "utils.hpp"
#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <boost/random/discrete_distribution.hpp>


// Constructor function - initialize member variables
// The functions following the colon is a initialization list
Neighbors::Neighbors(Params& params){
    
    
    // Calculate number of neighbors based on radius -Minus one to exclude self
    n_neighbors = (((2 * params.neighbor_radius) + 1) * ((2 * params.neighbor_radius) + 1)) - 1;
    
    // Resize / initialize vectors
    seed_rng.resize(n_neighbors);
    nn_sp.resize(n_neighbors);
    nn_index.resize(n_neighbors);
    nn_gen.resize(n_neighbors);
    nn_gen_1d_index.resize(n_neighbors);
    nn_gen_1d_index_dupe.resize(n_neighbors);
    
    
    seeds_by_sp.resize(params.n_sp_init);
    dead_seeds_sp.resize(params.n_sp_init);
    seeds_by_gen.resize(params.n_sp_init * params.n_alleles_init);
    probabilities.resize(params.n_sp_init * params.n_alleles_init);
    
    seeds_total = {0.0};
    
    // Initialize seed RNGS
    initSeedRNG(params);

}


//Initialize seed RNGS based on distance
void Neighbors::initSeedRNG(Params& params){
    
    
    /////////////////////////////////
    // Create RNG for seed dipsersal
    
    // alpha = 1 / mean dispersal distance in units of cells
    // dij = distance between source and endpoint
    // R = Radius of cells from source cell

    // Starting loop in top left corner and calculate distance
    
    int i = 0; // For counter
    
    float distance;
    
    for(int row = 0; row < (2 * params.neighbor_radius + 1); row++){
        
        for(int col = 0; col < (2 * params.neighbor_radius + 1); col++){
            
            
            if(row == params.neighbor_radius & col == params.neighbor_radius){
            //    std::cout << "XXXXXX" << "\t";
                continue; // Skip self
            }
           
            
            // If Global disperal, all individuals in 'neighborhood' (which are actually a sample of adults from the entire landscape), have an equal distance and probability of dispersal
            if(params.dispersal_mode == 1){
                
                distance = 1; //Distance is equal for all NN
                // Turn average seed dispersal distance to neighborhood radius so kernel is calculating with -1, which would make all probabilities low
                params.seed_disp_dist = params.neighbor_radius;
                
            } else { // For local dispersal
                
                distance = sqrt((params.neighbor_radius - row)*(params.neighbor_radius - row) + (params.neighbor_radius - col)*(params.neighbor_radius - col));
            }
            
            float disp_prob = neg_expo_discrete(1.0/params.seed_disp_dist,
                                                distance,
                                                params.neighbor_radius);
            
            std::binomial_distribution<int> seed_rng2(params.seeds_per_adult, disp_prob);
            
            seed_rng[i] = seed_rng2;

            std::cout << disp_prob << "\t";
            i++; // Increment count
            
        } // End col loop
        
        std::cout << "\n";
    } // End row loop
    
    
    // If global dispersal, set seed disp distance back to -1 to differentiate in output
    if(params.dispersal_mode == 1){
        params.seed_disp_dist = -1;
    }
    
    
}




// Updates neighbor indices in
void Neighbors::getNeighborIndex(int focal_cell, Params& params){
    // col = focal cell % width
    // row = focal cell / width
    // To go from x, y to array index: array_index = row * width + col
    
    // Initialize variables
    
    int i = 0; // Used for looping through nn_index
    int NN_index_temp;
    int row;
    int col;
    
    int width = params.width;
    int height = params.height;
    int neighbor_radius = params.neighbor_radius;
    
    // Row and col start is the row and col number of the TL neighbor
    int row_start = (focal_cell / width - neighbor_radius);
    int col_start = (focal_cell % width - neighbor_radius);
    
    // Start at Top left corner
    // Loop over rows and cols of neighborhood radius
    for(int row_ctr = 0; row_ctr < (2 * neighbor_radius + 1); row_ctr++){
        
        // Row = row number
        row = (row_start + row_ctr) % height;
        
        // Boundary check
        if(row < 0) row = row + height;
        
        // Loop through columns
        for(int col_ctr = 0; col_ctr < (2 * neighbor_radius + 1); col_ctr++){
            
            // Col = col number -- % width is boundary check
            col = (col_start + col_ctr) % width;
            
            // Boundary check
            if(col < 0) col = col + width;
            
            // 1d index of neighbor
            NN_index_temp = row * width + col;
           
            // Skip focal cell
            if(NN_index_temp == focal_cell){
                continue;
            }
         
            // Assignment to class
            nn_index[i] = NN_index_temp;
            
           // std::cout << NN_index_temp << " \t";
            
            i++; // Increment counter
       
            // Make sure NN is not out of bounds
                assert(NN_index_temp < params.area);
                assert(NN_index_temp >= 0);
        }
        
     //   std::cout << " \n";
    }
}


// Update neighbors information

void Neighbors::updateNeighbors(int focal_cell, Params& params, Sim& sim){
    
    // Update indices
    // Adds 1d array indices of neighbors to member var nn_index
    
    /////////////////
    // LOCAL DISPERSAL
    if(params.dispersal_mode == 0){
        Neighbors::getNeighborIndex(focal_cell, params);
    }
    
    /////////////////
    // GLOBAL DISPERAL
    // Will assign random 'neighbors' from anywhere on landscape, regardless of distance
    if(params.dispersal_mode == 1){
        for(auto& iter : nn_index){
            // Assign 'neighbor indices' from random cells on landscape
            iter = sim.cell_rng(sim.generator);
        }
    }
    
    // Use indices to update nn_sp and nn_gen
    int i = 0;
    for(auto& iter : nn_index){
        
        // Make sure sp number isn't crazy
        nn_sp[i] = sim.sp[iter];
        assert(sim.sp[iter] <= params.n_sp_init);
        
        // Make sure genotype is within range
        nn_gen[i] = sim.gen[iter];
        assert(sim.gen[iter] <= params.n_alleles_init);
        
        // Index for seeds_by_gen array
        nn_gen_1d_index[i] = sim.sp[iter] * params.n_alleles_init + sim.gen[iter];
        i++;
    }
    
    
}


// Disperse seeds from neighbors, based on negative exponential dispersal kernel
void Neighbors::disperseSeeds(std::mt19937& generator){
    
    int i = 0; // Counter
    int seeds_to_add;
    
    // Loop through neighbors
    for(auto& sp_index_iter : nn_sp ){
        
        seeds_to_add = seed_rng[i](generator);

        //Check to see if this has already been added
        // If it has, mark it as a duplicate to make CNDD and GNDD calculations
        // easier
        if(seeds_by_gen[nn_gen_1d_index[i]] > 0 ){
            nn_gen_1d_index_dupe[i] = true;
        }
        
        // Add seeds from that species and genotype combination
        seeds_by_gen[nn_gen_1d_index[i]] += seeds_to_add; // Add seeds for that specific genotype
        seeds_by_sp[nn_sp[i]] += seeds_to_add;
        
        i++; // Increment counter
        
    }
    
}



/////////
// GNDD
// Genotype dependent negative density dependence
// Reduces densities of seeds in seed_by_gen based on genotype density

void Neighbors::GNDD(std::vector<float>& gndd_sp){
    
    // Loop over species by gen counts and reduce densities
    // total number of iterations will equal number of neighbors
    
    for(int i = 0; i < nn_gen_1d_index.size(); i++){
        
        // If this is a duplicate - already reduced densities
        // Skip to avoid double reducing densities
        if(nn_gen_1d_index_dupe[i]){
            continue;
        }
        
        // Following equation of Harms et al. 2000
        // Where ## of recruits is a function of beta exponent and log + 1 density of seeds
        // Here, converting predicted log +1 density of recruits back to actual number of recruits with std::exp
        // Adding + 1 to recruit density makes it so theres always at least 2.7 possible recruits
        // Assume intercept (alpha) = 0

        
        seeds_by_gen[nn_gen_1d_index[i]] = std::exp(gndd_sp[nn_sp[i]] *
                                                         std::log(seeds_by_gen[nn_gen_1d_index[i]]));
        
       
        
    }
   
}





/////////
/// CNDD
// Based on densities of conspecifics

void Neighbors::CNDD(std::vector<float>& cndd_sp){
    
    
    // Loop over species and count how many seeds should die
    
    for(int i = 0; i < seeds_by_sp.size(); i++){
        
        // Check to make sure seeds don't equal 0
        // or else will take log of 0 and cause issues
        
        if(seeds_by_sp[i] == 0){
            dead_seeds_sp[i] = 0;
       
        } else {
   
        dead_seeds_sp[i] = seeds_by_sp[i] -
                            std::exp(cndd_sp[i] *
                            std::log(seeds_by_sp[i]));
        }
        
    }
    

    // Loop over species by gen counts and reduce densities
    // total number of iterations will equal number of neighbors
    
    for(int i = 0; i < nn_gen_1d_index.size(); i++){
        
        // If this is a duplicate - already reduced densities
        // Skip to avoid double reducing densities
        if(nn_gen_1d_index_dupe[i]){
            continue;
        }
        
        // Check for 0's or else might divide by 0
        
        if(seeds_by_sp[nn_sp[i]] == 0){
            continue;
        } else {
        
        seeds_by_gen[nn_gen_1d_index[i]] -= seeds_by_gen[nn_gen_1d_index[i]]/seeds_by_sp[nn_sp[i]] * dead_seeds_sp[nn_sp[i]];
        }

    }
} // End function


/////////
// NDD
// Simultaneous NDD that reduces density based on GNDD and CNDD in one formula

void Neighbors::NDD(std::vector<float>& gndd_sp, std::vector<float>& cndd_sp){
    
    // Loop over species by gen counts and reduce densities
    // total number of iterations will equal number of neighbors
    
    for(int i = 0; i < nn_gen_1d_index.size(); i++){
        
        // If this is a duplicate - already reduced densities
        // Skip to avoid double reducing densities
        if(nn_gen_1d_index_dupe[i]){
            continue;
        }
        
        // Check for 0's or else might divide by 0 and return NAN
        if(seeds_by_sp[nn_sp[i]] == 0 || seeds_by_gen[nn_gen_1d_index[i]] == 0){
            continue;
        } else {
        
        // Following equation of Harms et al. 2000
        // Where ## of recruits is a function of beta exponent and log density of seeds
        // Here, converting predicted log density of recruits back to actual number of recruits with std::exp
        // Assume intercept (alpha) = 0
        // based off this formula in R -- R = (bgi * log(sgi) + bci * log(sci) + log(sgi/sci)) / 2
        //seeds_by_gen[nn_gen_1d_index[i]] = std::exp((gndd_sp[nn_sp[i]] *
        //                                            std::log(seeds_by_gen[nn_gen_1d_index[i]]) +
        //                                            cndd_sp[nn_sp[i]] * std::log(seeds_by_sp[nn_sp[i]]) +
        //                                             std::log(seeds_by_gen[nn_gen_1d_index[i]]/seeds_by_sp[nn_sp[i]])) / 2);
            
        // Standard harms formula, with rounding down to the nearest seed
            seeds_by_gen[nn_gen_1d_index[i]] = std::floor(std::exp(cndd_sp[nn_sp[i]] * std::log(seeds_by_sp[nn_sp[i]]) +
                                                     std::log(seeds_by_gen[nn_gen_1d_index[i]]/seeds_by_sp[nn_sp[i]])));
            
        } // End else statement
    } // End loop over neighbors
} // End function definition


////
// Count the total number of seeds in the cell

void Neighbors::totalSeeds(){
    
    for(int i = 0; i < seeds_by_gen.size(); i++){
    
        seeds_total += seeds_by_gen[i];
    }
    
    // Make sure there at least some seeds
    assert(seeds_total > 0);
}


// Choose winner
// Lottery based on relative frequency in cell
// Then assigns species and genotype to focal cell

void Neighbors::chooseWinner(int focal_cell, Params& params, Sim& sim){

    
    // Assign probs based on relative frequency
    int i = 0;
    
    for(auto& iter : nn_gen_1d_index){
        
        // Skip if duplicate to avoid double counting
        if(nn_gen_1d_index_dupe[i]){
            i++;
            continue;
        }
        
        probabilities[nn_gen_1d_index[i]] = seeds_by_gen[nn_gen_1d_index[i]] / seeds_total;
        
            i++;
    }
    
    
    // Choose species to establish in empty cell
    // Choose a winner / seed to establish based on relative frequency and weighted probability
    boost::random::discrete_distribution<> seed_winner_rng(probabilities.begin(), probabilities.end());
    
    
    // Winner index is index of probabilities - which is species x gen
    int winner_index; // Move up and out of loop?
    winner_index = seed_winner_rng(sim.generator);
    
    // Reassign Species ID and genotype for 'winner' then continue looping through empty cells
    
    // Error check
    assert((winner_index / params.n_alleles_init) <= params.n_sp_init);
    assert((winner_index % params.n_alleles_init) <= params.n_alleles_init);
    
    // Have to do modulus and division operation in indexing to find species and genotype
    
    sim.sp[focal_cell]  = winner_index / params.n_alleles_init;
    sim.gen[focal_cell] = winner_index % params.n_alleles_init;
    
//    // Print winner
//    std::cout<< "Winner - species: " << nn_sp[winner_index / n_alleles_init] << " \n";
    
}





// Print status of number of seeds per neighbor

void Neighbors::printStatus(Params& params){
    
    
    ////////////////////////
    // Grid of neighbors
    std::cout << "\n| ----- Neighbors (species.gen) -------- | \n";
    
    int i = 0;
    for(int row_ctr = 0; row_ctr < (2 * params.neighbor_radius + 1); row_ctr++){
       
        for(int col_ctr = 0; col_ctr < (2 * params.neighbor_radius + 1); col_ctr++){
            
            // For self
            if(col_ctr == params.neighbor_radius & row_ctr == params.neighbor_radius){
                std::cout << "-- \t";
                continue;
            }
            
            std::cout << nn_sp[i] <<  "\t";
            //"." << nn_gen[i] << ";

            i++;
        }
        
        std::cout << "\n";
    } // End neighbor print
    
    
    ////////////////////////
    // Seeds by species
    std::cout << "\n| ----- Seeds by species -------- | \n";
    i = 0;
    for(auto& iter : seeds_by_sp){
        std::cout << "Species " << i << ": " << iter << "\n";
        i++;
    }
    
    
    ////////////////////////
    // Seeds by gen
    std::cout << "\n| ----- Seeds by species by gen -------- | \n";
    i = 0;
    
    // Column labels
    std::cout<< "Gen: \t\t";
    for(int col = 0; col < params.n_alleles_init; col++){
        std::cout << col << " \t";
        }
    std::cout << "\n";

    for(int row = 0; row < params.n_sp_init; row++){
        std::cout << "Species " << row << ": \t";
        
        for(int col = 0; col < params.n_alleles_init; col++){
            std::cout << seeds_by_gen[i] << "\t";
            
            i++;
        }
            std::cout << "\n";
    }
  
}


// Reset - reset values back to 0

void Neighbors::reset(){
    
    std::fill(nn_gen_1d_index_dupe.begin(),
              nn_gen_1d_index_dupe.end(),
              false);
    std::fill(seeds_by_sp.begin(), seeds_by_sp.end(), 0);
    std::fill(dead_seeds_sp.begin(), dead_seeds_sp.end(), 0);
    std::fill(seeds_by_gen.begin(), seeds_by_gen.end(), 0);
    std::fill(probabilities.begin(), probabilities.end(), 0);
    
    seeds_total = 0;
    
}







