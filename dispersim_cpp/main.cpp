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
#include "spatial.hpp"

int main(int argc, const char * argv[]) {

////////////////////////
////////////////////////
// PARAMETERS
    
    // Number of generations
    int steps = 1000;
    float mortality_rate = 0.1; // Proportion of landscape dieing per step
    
    // Species parameters
    int n_sp_init = 100;
    int n_alleles_init = 10;
    
    
    // Landscape parameters
    int width  = 3;
    int height = 3;
    int area = width * height;
    
    int n_dead_per_step = mortality_rate * area;
    int empty_cell_indices[n_dead_per_step];

  
////////////////////////
////////////////////////
// Initialize simulation
    
    // Random number generators
    // Good info on generating random numbers in C++
    // http://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c
    
    std::random_device device;
    std::mt19937 generator(device());
    
    // RNG for cell
    std::uniform_int_distribution<int> cell_rng(0, area);
    
    // Add species to landcape
    std::uniform_int_distribution<int> species_rng(0, n_sp_init - 1); // Init RNG
    
    std::cout << "Initialize species vector with " << area << " inds \n";
    std::vector<int> sp(area); // Init vector that will hold species
    
    for (auto& iter : sp){
        iter = species_rng(generator);
    }

    
    // Add genotypes to species
    
    //...
    
    
    
    
    
    
////////////////////////
////////////////////////
    
// Begin simulation
    
    
    for (int step = 0; step < steps; step++){
        
        
    // Death of adults
    // A set proportion dies every step, dependent on mortality rate

        for(auto& empty_cell_iter : empty_cell_indices) {
            empty_cell_iter = cell_rng(generator);
        }
        
        // Loop through empty cells and assign new species
        
        for(auto& empty_cell_iter : empty_cell_indices){
            
            
           // sp[empty_cell_iter] = species_rng(generator);
            
            
            // Disperse seeds into each empty cell
            
            //
            
            
            
        }
        
        
      //  std::cout << "Species in cell 0 – " << sp[0] << "\n";
        
    } // End step loop

    findNN(0, height, width, area);

    
    return 0;
}
