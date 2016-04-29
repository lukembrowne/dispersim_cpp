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
    
    
    // Landscape parameters
    int width  = 1000;
    int height = 1000;
    int area = width * height;
    
    int dispersal_mode = 0; // 1 == global; 0 == local
    
    int n_dead_per_step = mortality_rate * area;
    int empty_cell_indices[n_dead_per_step];
    
    
    // Lagniappe parameters
    
    int print_every_n_steps = 50;
    
    

  
////////////////////////
////////////////////////
// Initialize simulation
    
    
    std::cout << "| ---------------------------- | \n";
    std::cout << "Steps: " << steps << " or ~ " << steps * 5 << " years \n";
    std::cout << "Area: " << (width * 5 * height * 5)/10000 << " ha \n";

    
    
    
    // Random number generators
    // Good info on generating random numbers in C++
    // http://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c
    
    std::random_device device;
    std::mt19937 generator(device());
    
    // RNG for cell
    std::uniform_int_distribution<int> cell_rng(0, (area-1));
    
    // RNG
    
    // Add species to landcape
    std::uniform_int_distribution<int> species_rng(0, n_sp_init - 1); // Init RNG
    
    std::cout << "Individuals: " << area << " \n";
    std::vector<int> sp(area); // Init vector that will hold species
    
    for (auto& iter : sp){
        iter = species_rng(generator);
    }

    
    // Add genotypes to species
    
    //...
    
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
    
    
    
////////////////////////
////////////////////////
    
// Begin simulation
    
    std::cout << "\n\n";
    
  
    for (int step = 1; step <= steps; step++){
        
        
    // Death of adults
    // A set proportion dies every step, dependent on mortality rate

        for(auto& empty_cell_iter : empty_cell_indices) {
            empty_cell_iter = cell_rng(generator);
        }
        
        // Loop through empty cells and assign new species
        
        for(auto& empty_cell_iter : empty_cell_indices){
            
            // GLOBAL DISPERAL
            if(dispersal_mode == 1){
                sp[empty_cell_iter] = sp[cell_rng(generator)]; // Random species
                continue;
            }
        
            // LOCAL DISPERSAL
            if(dispersal_mode == 0){

              
                // Find surrounding 8 neighbors
                // Returns position in array
                // Use position to look up species later
                std::vector<int> neighbors = findNN(empty_cell_iter,
                                                    height, width, area);
                
                
                // Disperse seeds from 8 neighbors into empty cell
                
                // Make a map that holds species - ## of seeds pair
                std::map <int, int> sp_seed;
            
                double seeds_total{0}; // Sum up total seeds to calculate relative frequency later
                
               for(auto iter : neighbors){
                   
                   sp_seed[sp[iter]] = seed_rng(generator) + sp_seed[sp[iter]];
                   
                   seeds_total = seeds_total + sp_seed[sp[iter]];
                   
                } // End looping over neighbors
                
                
                // HERE IS WHERE NDD CODE WOULD GO
                // REDUCING NUMBER OF SEEDS BASED ON DENSITY AND GENOTYPE
                // REDUCE DENSITY BASED ON NEARBY ADULTS OR BASED ON SEED DENSITY
                
                
                // Choose species to establish in empty cell
            
                
                    // Assign probs based on relative frequency
                std::vector<double> probabilities(sp_seed.size());
                int sp_names[sp_seed.size()];
                int i = 0;
                for(auto& iter : sp_seed){
                    probabilities[i] = iter.second / seeds_total;
                    sp_names[i] = iter.first;
                    i++;
               //     std::cout << probabilities[i-1] << " ";
                }
                
            // Choose a winner / seed to establish based on relative frequency and weighted probability
               boost::random::discrete_distribution<> seed_winner_rng(probabilities);
               sp[empty_cell_iter] = sp_names[seed_winner_rng(generator)]; // Reassign species
                
            } // End looping over empty cells
            
        } // End local dispersal if

        if(step % print_every_n_steps == 0 | step == 1 ){
            std::cout << "Step:" << step << " | " << calcSpeciesRichness(sp) ;
        }
                       
    } // End step loop


    
    return 0;
}
