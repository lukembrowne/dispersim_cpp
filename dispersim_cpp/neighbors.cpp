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
Neighbors::Neighbors(int neighbor_radius, int n_sp_init, int n_alleles_init){
    
    n_neighbors = (((2 * neighbor_radius) + 1) * ((2 * neighbor_radius) + 1)) - 1; // Minus one to exclude self
    
    // Resize / initialize vectors
    seed_rng.resize(n_neighbors);
    nn_sp.resize(n_neighbors);
    nn_index.resize(n_neighbors);
    nn_gen.resize(n_neighbors);
    nn_gen_1d_index.resize(n_neighbors);
    nn_gen_1d_index_dupe.resize(n_neighbors);
    
    
    seeds_by_sp.resize(n_sp_init);
    dead_seeds_sp.resize(n_sp_init);
    seeds_by_gen.resize(n_sp_init * n_alleles_init);
    probabilities.resize(n_sp_init * n_alleles_init);
    
    seeds_total = {0.0};

}


//Initialize seed RNGS based on distance
void Neighbors::initSeedRNG(int neighbor_radius,
                            int seed_disp_dist,
                            int seeds_per_adult){
    
    
    /////////////////////////////////
    // Create RNG for seed dipsersal
    
    // alpha = 1 / mean dispersal distance in units of cells
    // dij = distance between source and endpoint
    // R = Radius of cells from source cell

    // Starting loop in top left corner and calculate distance
    
    int i = 0; // For counter
    
    for(int row = 0; row < (2 * neighbor_radius + 1); row++){
        
        for(int col = 0; col < (2 * neighbor_radius + 1); col++){
            
            
            if(row == neighbor_radius & col == neighbor_radius){
                continue; // Skip self
            }
           
            float distance = sqrt((neighbor_radius - row)*(neighbor_radius - row) + (neighbor_radius - col)*(neighbor_radius - col)); // ...
            
            float disp_prob = neg_expo_discrete(1.0/seed_disp_dist,
                                                distance,
                                                neighbor_radius);
            
            std::binomial_distribution<int> seed_rng2(seeds_per_adult, disp_prob);
            
            seed_rng[i] = seed_rng2;

          //  std::cout << disp_prob << "\t";
            i++; // Increment count
            
        } // End col loop
        
       // std::cout << "\n";
    } // End row loop
    
    int test = 0;
   

}




// Updates neighbor indices in
void Neighbors::getNeighborIndex(int focal_cell, int height,
                                 int width, int area,
                                 int neighbor_radius){
    // col = focal cell % width
    // row = focal cell / width
    // To go from x, y to array index: array_index = row * width + col
    
    // Initialize variables
    
    int i = 0; // Used for looping through nn_index
    int NN_index_temp;
    int row;
    int col;
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
                assert(NN_index_temp < area);
                assert(NN_index_temp >= 0);
        }
        
     //   std::cout << " \n";
    }
}


// Update neighbors information
// Could clean up if we had a Sim class that held sp and gen and information about simulation like height width, area, neighborhood radius etc

void Neighbors::updateNeighbors(int focal_cell, int height,
                                int width, int area,
                                int neighbor_radius,
                                std::vector<int>& sp,
                                std::vector<int>& gen,
                                int n_alleles_init,
                                int n_sp_init){
    
    // Update indices
    // Adds 1d array indices of neighbors to member var nn_index
    Neighbors::getNeighborIndex(focal_cell, height, width, area, neighbor_radius);
    
    // Use indices to update nn_sp and nn_gen
    int i = 0;
    for(auto& iter : nn_index){
        
        // Make sure sp number isn't crazy
        nn_sp[i] = sp[iter];
        assert(sp[iter] <= n_sp_init);
        
        // Make sure genotype is within range
        nn_gen[i] = gen[iter];
        assert(gen[iter] <= n_alleles_init);
        
        // Index for seeds_by_gen array
        nn_gen_1d_index[i] = sp[iter] * n_alleles_init + gen[iter];
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
        
        // Add to total seeds per species
        
        seeds_by_sp[nn_sp[i]] += seeds_by_gen[nn_gen_1d_index[i]];
        
        
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

////
// Count the total number of seeds in the cell

void Neighbors::totalSeeds(){
    
    for(int i = 0; i < seeds_by_gen.size(); i++){
    
        seeds_total += seeds_by_gen[i];
    }
    
    assert(seeds_total >= 0);
}


// Choose winner
// Lottery based on relative frequency in cell
// Then assigns species and genotype to focal cell

void Neighbors::chooseWinner(std::mt19937& generator,
                             int focal_cell,
                             std::vector<int>& sp,
                             std::vector<int>& gen,
                             int n_sp_init,
                             int n_alleles_init){

    
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
//    
//    std::discrete_distribution<float> test{probabilities.begin(), probabilities.end()};
//    
//    for(int i = 0; i < 100; i++){
//        std::cout << seed_winner_rng(generator) << "\n";
//    }
//    
    
    
    // Winner index is index of probabilities - which is species x gen
    int winner_index; // Move up and out of loop?
    winner_index = seed_winner_rng(generator);
    
    // Reassign Species ID and genotype for 'winner' then continue looping through empty cells
    
    // Error check
    assert((winner_index / n_alleles_init) <= n_sp_init);
    assert((winner_index % n_alleles_init) <= n_alleles_init);
    
    // Have to do modulus and division operation in indexing to find species and genotype
    
    sp[focal_cell]  = winner_index / n_alleles_init;
    gen[focal_cell] = winner_index % n_alleles_init;
    
//    // Print winner
//    std::cout<< "Winner - species: " << nn_sp[winner_index / n_alleles_init] << " \n";
    
}





// Print status of number of seeds per neighbor

void Neighbors::printStatus(int neighbor_radius, int n_sp_init, int n_alleles_init){
    
    
    // Graph of neighbors
    std::cout << "\n| ----- Neighbors (species.gen) -------- | \n";
    
    
    int i = 0;
    for(int row_ctr = 0; row_ctr < (2 * neighbor_radius + 1); row_ctr++){
       
        for(int col_ctr = 0; col_ctr < (2 * neighbor_radius + 1); col_ctr++){
            
            // For self
            if(col_ctr == neighbor_radius & row_ctr == neighbor_radius){
                std::cout << "-- \t";
                continue;
            }
            
            std::cout << nn_sp[i] <<  "\t";
            //"." << nn_gen[i] << ";

            i++;
        }
        
        std::cout << "\n";
    } // End neighbor print
    
    // Seeds by species
    std::cout << "\n| ----- Seeds by species -------- | \n";
    i = 0;
    for(auto& iter : seeds_by_sp){
        std::cout << "Species " << i << ": " << iter << "\n";
        i++;
    }
    
    // Seeds by gen
    std::cout << "\n| ----- Seeds by species by gen -------- | \n";
    i = 0;
    
    // Column labels
    std::cout<< "Gen: \t\t";
    for(int col = 0; col < n_alleles_init; col++){
        std::cout << col << " \t";
        }
    std::cout << "\n";

    for(int row = 0; row < n_sp_init; row++){
        std::cout << "Species " << row << ": \t";
        
        for(int col = 0; col < n_alleles_init; col++){
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







