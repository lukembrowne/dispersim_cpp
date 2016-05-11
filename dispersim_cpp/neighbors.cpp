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
            
            this->seed_rng[i] = seed_rng2;

            
            i++; // Increment count
        }
        
    }

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
            this->nn_index[i] = NN_index_temp;
            
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
                                int n_alleles_init){
    
    // Update indices
    Neighbors::getNeighborIndex(focal_cell, height, width, area, neighbor_radius);
    
    // Use indices to update nn_sp and nn_gen
    int i = 0;
    for(auto& iter : this->nn_index){
        this->nn_sp[i] = sp[iter];
        this->nn_gen[i] = gen[iter];
        // Index for seeds_by_gen array
        this->nn_gen_1d_index[i] = sp[iter] * n_alleles_init + gen[iter];
        i++;
    }
    
    
}


// Disperse seeds from neighbors, based on negative exponential dispersal kernel
void Neighbors::disperseSeeds(std::mt19937& generator){
    
    int i = 0; // Counter
    int seeds_to_add;
    
    // Loop through neighbors
    for(auto& sp_index_iter : this->nn_sp ){
        
        seeds_to_add = this->seed_rng[i](generator);

        //Check to see if this has already been added
        // If it has, mark it as a duplicate to make CNDD and GNDD calculations
        // easier
        if(this->seeds_by_gen[this->nn_gen_1d_index[i]] > 0 ){
            this->nn_gen_1d_index_dupe[i] = true;
        }
        
        // Add seeds from that species and genotype combination
        this->seeds_by_gen[this->nn_gen_1d_index[i]] += seeds_to_add; // Add seeds for that specific genotype
        
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
    
    for(int i = 0; i < this->nn_gen_1d_index.size(); i++){
        
        // If this is a duplicate - already reduced densities
        // Skip to avoid double reducing densities
        if(this->nn_gen_1d_index_dupe[i]){
            continue;
        }
        
        // Following equation of Harms et al. 2000
        // Where ## of recruits is a function of beta exponent and log + 1 density of seeds
        // Here, converting predicted log +1 density of recruits back to actual number of recruits with std::exp
        // Adding + 1 to recruit density makes it so theres always at least 2.7 possible recruits
        // Assume intercept (alpha) = 0

        
        this->seeds_by_gen[this->nn_gen_1d_index[i]] = std::exp(gndd_sp[this->nn_sp[i]] *
                                                         std::log(this->seeds_by_gen[this->nn_gen_1d_index[i]]));
        
        // Add to total seeds per species
        
        this->seeds_by_sp[this->nn_sp[i]] += this->seeds_by_gen[this->nn_gen_1d_index[i]];
        
        
    }
   
}





/////////
/// CNDD
// Based on densities of conspecifics

void Neighbors::CNDD(std::vector<float>& cndd_sp){
    
    
    // Loop over species and count how many seeds should die
    
    for(int i = 0; i < seeds_by_sp.size(); i++){
        
        dead_seeds_sp[i] = this->seeds_by_sp[i] -
                            std::exp(cndd_sp[i] *
                            std::log(seeds_by_sp[i]));
        
    }
    

    // Loop over species by gen counts and reduce densities
    // total number of iterations will equal number of neighbors
    
    for(int i = 0; i < this->nn_gen_1d_index.size(); i++){
        
        // If this is a duplicate - already reduced densities
        // Skip to avoid double reducing densities
        if(this->nn_gen_1d_index_dupe[i]){
            continue;
        }
        

        this->seeds_by_gen[this->nn_gen_1d_index[i]] -= this->seeds_by_gen[this->nn_gen_1d_index[i]]/seeds_by_sp[this->nn_sp[i]] * dead_seeds_sp[this->nn_sp[i]];

    }
} // End function

////
// Count the total number of seeds in the cell

void Neighbors::totalSeeds(){
    
    for(int i = 0; i < this->seeds_by_gen.size(); i++){
    
        this->seeds_total += this->seeds_by_gen[i];
    }
}


// Choose winner

int Neighbors::chooseWinner(std::mt19937& generator){

    
    // Assign probs based on relative frequency

    int i = 0;
    
    for(auto& iter : this->nn_gen_1d_index){
        
        // Skip if duplicate to avoid double counting
        if(this->nn_gen_1d_index_dupe[i]){
            i++;
            continue;
        }
        
        
        this->probabilities[this->nn_gen_1d_index[i]] = this->seeds_by_gen[this->nn_gen_1d_index[i]] / this->seeds_total;
        
            i++;
    }
    
    
    // Choose species to establish in empty cell
    // Choose a winner / seed to establish based on relative frequency and weighted probability
    boost::random::discrete_distribution<> seed_winner_rng(this->probabilities);
    
    int winner_index; // Move up and out of loop?
    winner_index = seed_winner_rng(generator);

    
    // Print winner to console

        std::cout << "\n\n Winner is.. Sp: " << this->nn_sp[winner_index] <<
        " | Genotype: " << this->nn_gen[winner_index] << "\n\n\n" <<
        "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n\n";

    
    return(winner_index);
    
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
            
            std::cout << this->nn_sp[i] << "." << this->nn_gen[i] << " \t";

            i++;
        }
        
        std::cout << "\n";
    } // End neighbor print
    
    // Seeds by species
    std::cout << "\n| ----- Seeds by species -------- | \n";
    i = 0;
    for(auto& iter : this->seeds_by_sp){
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
            std::cout << this->seeds_by_gen[i] << "\t";
            
            i++;
        }
            std::cout << "\n";
    }
  
}





