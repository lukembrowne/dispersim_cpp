//
//  main.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/28/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include "spatial.hpp"
#include "utils.hpp"


int main(int argc, const char * argv[]) {

////////////////////////
////////////////////////
// PARAMETERS
    
    // Number of generations
    int steps = 1000; // Each step = 5-10 years with 0.10 % mortality rate
    float mortality_rate = 0.1; // Proportion of landscape dying per step
    
    // Species parameters
    int n_sp_init = 100;
    int n_alleles_init = 10;
    float seed_disp_dist = 5; // In units of cells
    int seeds_per_adult = 500; // Equal to fecundity..
    
    // NDD parameters
    float min_ndd = 1; // Min must be greater than max
    float max_ndd = 0;
    
    
    // Landscape parameters
    int width  = 250;
    int height = 250;
    int area = width * height;
    
    int dispersal_mode = 1; // 1 == global; 0 == local
    
    int n_dead_per_step = mortality_rate * area;
    int empty_cell_indices[n_dead_per_step]; // Initialize
    
    
    // Lagniappe parameters
    
    int print_every_n_steps = 50;
    
   /// NEED TO ADD IN MIGRATION
    
    // boost::random::uniform_01 - generate uniform
    // Can set migration based on BCI paper in Plos one condit 2012

    
////////////////////////
////////////////////////
// Initialize simulation
    
    
    // Print out information about simulation
    std::cout << "| ---------------------------- | \n";
    std::cout << "Steps: " << steps << " or ~ " << steps * 5 << " years .. " << steps / 10 <<" generations \n";
    std::cout << "Area: " << (width * 5 * height * 5)/10000 << " ha \n";

    
    
    
    // Initialize Random number generators
    // Good info on generating random numbers in C++
    // http://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c
    
    std::random_device device;
    std::mt19937 generator(device());
    
    // RNG that randomly chooses cells in landscape - used for mortality algorithm
        std::uniform_int_distribution<int> cell_rng(0, (area-1));
    
    // RNG for species - ranging from 0 to # of species - 1 (species ID used for indexing)
        std::uniform_int_distribution<int> species_rng(0, n_sp_init - 1);
        
        std::cout << "Individuals: " << area << " .. " << area/n_sp_init << " per species \n";
        
    // Initialize species vector and place species randomly throughout landscape
        std::vector<int> sp(area); // Init vector that will hold species
        
        for (auto& iter : sp){
            iter = species_rng(generator);
        }

    
    // Add genotypes randomly to species
    
        // Initialize genotype RNG
        std::uniform_int_distribution<int> gen_rng(0, n_alleles_init);

        std::vector<int> gen(area); // Init vector that will hold species
    
        for (auto& iter : gen){
            iter = gen_rng(generator);
        }
    
    // Assign NDD values to species
    // Lower species id == Stronger NDD
    
    std::vector<float> ndd_sp(n_sp_init);
    
    float ndd_increment = std::abs((max_ndd-min_ndd))/(float)n_sp_init;
    int i = 0;
    
    // If min and max are the same...
    if(min_ndd == max_ndd){
        for (auto& iter : ndd_sp){
            iter = max_ndd;
        }
    } else {
    
    // If variation in NDD (min_ndd != max_ndd)...
    for(float ndd_val = max_ndd; ndd_val< min_ndd; ndd_val += ndd_increment){
        ndd_sp[i] = ndd_val;
        i++;
        if(i == (n_sp_init-1)) break;
        }
    }
    
    // Create RNG for seed dipsersal
    
        // Assumes all neighboring cells have the same probability..
        // which will probably need to be changed later
        // alpha = 1 / mean dispersal distance
        // dij = distance between source and endpoint
        // R = Radius of cells from source cell

    float disp_prob = neg_expo_discrete(1/seed_disp_dist, 1, 1);
    
    std::binomial_distribution<int> seed_rng(seeds_per_adult, disp_prob);
    
//    for(int i = 0; i < 50; i ++){
//        std::cout << seed_rng(generator) << "\t";
//    }
    
    //Initialize vectors for survival calculations

    std::vector<int> neighbors(8); // Change from 8 if doing more than 8 nearest neighbors
    std::vector<int> seeds(n_sp_init, 0.0); // Initialize to 0
    int seeds_total = {0};
    
////////////////////////
////////////////////////
    
// Begin simulation
    
    std::cout << "\n\n";
    
  
    for (int step = 1; step <= steps; step++){
        
        
    // Death of adults
    // A set proportion dies every step, dependent on mortality rate
        
        
        // Choose which cells will die
        for(auto& empty_cell_iter : empty_cell_indices) {
            empty_cell_iter = cell_rng(generator);
        }
        
        // Loop through empty cells and assign new species
        
        for(auto& empty_cell_iter : empty_cell_indices){
            
            // GLOBAL DISPERAL
            
                // Assigns species to empty cell based on relative frequency in population
                // May need to exclude dead trees from being chosen..
                if(dispersal_mode == 1){
                    sp[empty_cell_iter] = sp[cell_rng(generator)]; // Random species
                    continue;
                }
        
            // LOCAL DISPERSAL
            if(dispersal_mode == 0){

              
                // Find surrounding 8 neighbors
                // Returns position in array
                // Use position to look up species later
                 neighbors = findNN(empty_cell_iter, height, width, area);
                
                
                // Disperse seeds from 8 neighbors into empty cell
                
                // Loop over neighbors
                std::vector<int> nn_sp_key(8);
                int i = 0;
                for(auto iter : neighbors){
                    
                    // Generate species key
                    // Hopefully more efficient than a map
                    nn_sp_key[i] = sp[iter];
                    
                    // Vector holding number of seeds contributed by each NN
                    
                    seeds[sp[iter]] += seed_rng(generator);
                    //std::cout << nn_seeds[i] << " ";
                    i++;
                   
                }
                   //  std::cout << "\n";
                
                // Remove duplicate species in species key list
                std::sort(nn_sp_key.begin(), nn_sp_key.end());
                nn_sp_key.erase(std::unique(nn_sp_key.begin(), nn_sp_key.end()), nn_sp_key.end());
                
                
                
                
                // NDD algorithm
                // Reduce number of seeds based on density
                seeds_total = 0; // Reset total to 0
                
                for(auto& nn_sp : nn_sp_key){
                    seeds[nn_sp] = seeds[nn_sp] * ndd_sp[nn_sp]; // Reduce by half
                    
                    seeds_total += seeds[nn_sp]; // Add to total seeds
                    
                }

                
                    // Assign probs based on relative frequency
                std::vector<float> probabilities(nn_sp_key.size());
                i = 0;
                
                for(auto& iter : nn_sp_key){
                    probabilities[i] = (float)seeds[iter]/(float)seeds_total; // Need to cast as float or probability is 0
                    
                    seeds[iter] = 0; // Reset seeds count
                    
                    i++;
                }
                
                
                // Choose species to establish in empty cell
                // Choose a winner / seed to establish based on relative frequency and weighted probability
                   boost::random::discrete_distribution<> seed_winner_rng(probabilities);
                   sp[empty_cell_iter] = nn_sp_key[seed_winner_rng(generator)]; // Reassign species

                
            } // End looping over empty cells
            
        } // End local dispersal if

        if(step % print_every_n_steps == 0 | step == 1 ){
            std::cout << "Step:" << step << " | " << calcSpeciesRichness(sp) ;
        }
                       
    } // End step loop
    
    writeLandscape(sp, height, width);


    
    return 0;
}
