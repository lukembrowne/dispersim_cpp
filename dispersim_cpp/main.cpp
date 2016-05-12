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
#include "utils.hpp"
#include "summary.hpp"
#include "neighbors.hpp"
#include "sim.hpp"
#include <iomanip>
#include <sstream>



/// MAIN

int main(int argc, const char * argv[]) {

////////////////////////
////////////////////////
// PARAMETERS
    
    // Something strange happening where landscapes get vertical striping in species, ex when radius = 5 and disp = 3
    
    // Looks find when radius = 1, but banding happens when radius = 2
    
    // Lagniappe parameters
    int save_every_n_steps = 25;
    bool verbose = false;
    
    // Number of generations
    int steps = 300; // Each step = 5-10 years with 0.10 % mortality rate
    float mortality_rate = 0.25; // Proportion of landscape dying per step
    
    // Species parameters
    int n_sp_init = 50;
    int n_alleles_init = 2;
    float seed_disp_dist = 3; // In units of cells
    int neighbor_radius = 5; // Miranda = 5; Banitz = 20

    int seeds_per_adult = 500; // Equal to fecundity..
    
    // NDD parameters
    float max_cndd = 1.;
    float min_cndd = 1.; // Min must be greater numerically than max, but means weaker NDD
   
    float max_gndd = 1.0;
    float min_gndd = 1.0; // Min must be greater numerically than max, but means weaker NDD
    
    
    // Landscape parameters
    int width  = 100;
    int height = 100;
    int area = width * height;
    
    float migration_rate = 0.0001; // Immigrant per recruit (~1 in 10,000) is from BCI paper; 1 in 9000 used in Muller Landau 2007
    
    int dispersal_mode = 0; // 1 == global; 0 == local
    
    int n_dead_per_step = mortality_rate * area;
    std::vector<int> empty_cell_indices(n_dead_per_step); // Initialize
    
    
    std::vector<Summary_step> summary_over_time; // Initialize
 
    
    
    
    
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
    
    // RNG for Migration
    std::uniform_real_distribution<> migration_rng(0, 1);

    
    // RNG for species - ranging from 0 to # of species - 1 (species ID used for indexing)
        std::uniform_int_distribution<int> species_rng(0, n_sp_init - 1);
        
        std::cout << "Individuals: " << area << " .. " << area/n_sp_init << " per species \n";
        
    // Initialize species vector and place species randomly throughout landscape
        std::vector<int> sp(area); // Init vector that will hold species
        
        for (auto& iter : sp){
            iter = species_rng(generator);
            assert(iter <= n_sp_init);
        }

    
    // Add genotypes randomly to species
    
        // Initialize genotype RNG
        std::uniform_int_distribution<int> gen_rng(0, (n_alleles_init-1));

        std::vector<int> gen(area); // Init vector that will hold species
    
        for (auto& iter : gen){
            iter = gen_rng(generator);
        }
    
    
    
    
    /////////////////////////////////
    // Assign CNDD values to species - conspecific density dependence
    // Lower species id == Stronger NDD
    
    std::vector<float> cndd_sp(n_sp_init);
    
    float cndd_increment = std::abs((max_cndd-min_cndd))/(float)n_sp_init;
    int i = 0;
    
    // If min and max are the same...
    if(min_cndd == max_cndd){
        for (auto& iter : cndd_sp){
            iter = max_cndd;
        }
    } else {
    
    // If variation in NDD (min_ndd != max_ndd)...
    for(float cndd_val = max_cndd; cndd_val <= min_cndd; cndd_val += cndd_increment){
        cndd_sp[i] = cndd_val;
        i++;
        if(i == n_sp_init) break; // Make sure it doensn't go out of bounds
        }
    }
    
    
    /////////////////////////////////
    // GENOTYPE DEPENDENT NDD (GNDD)
    
    std::vector<float> gndd_sp(n_sp_init);
    
    float gndd_increment = std::abs((max_gndd-min_gndd))/(float)n_sp_init;
    i = 0;
    
    // If min and max are the same...
    if(min_gndd == max_gndd){
        for (auto& iter : gndd_sp){
            iter = max_gndd;
        }
    } else {
        
        // If variation in NDD (min_ndd != max_ndd)...
        for(float gndd_val = max_gndd; gndd_val <= min_gndd; gndd_val += gndd_increment){
            gndd_sp[i] = gndd_val;
            i++;
            if(i == n_sp_init) break; // Make sure it doensn't go out of bounds
        }
    }
    

    
    // Initialize Neighbors object
    
    Neighbors neighbors(neighbor_radius,
                        n_sp_init, n_alleles_init);
    
    neighbors.initSeedRNG(neighbor_radius,
                           seed_disp_dist,
                           seeds_per_adult);
    

////////////////////////
////////////////////////
    
// Begin simulation
    
    std::cout << "\n\n";
    
    for (int step = 1; step <= steps; step++){
        
        
    // Death of adults
    // A set proportion dies every step, dependent on mortality rate
        
        // Choose which cells will die
        for(int i = 0; i < empty_cell_indices.size(); i++) {
            empty_cell_indices[i] = cell_rng(generator);
        }
        
        // Loop through empty cells and assign new species
 
        for(auto& empty_cell_iter : empty_cell_indices){

            // MIGRATION -
            
            if(migration_rng(generator) <= migration_rate){
            
                // Replace with random species with random genotype
                
                sp[empty_cell_iter] = species_rng(generator);
                gen[empty_cell_iter] = gen_rng(generator);
            
                
                continue; // Jump to next empty cell
            }
            
            // GLOBAL DISPERAL
            // Need to change so that some proportion of adults contribute seedlings, but these are still evaluated for NDD effects.. just that dispersal is not limited by proximity - randomly choose 8 adults to contribute seeds
            
                // Assigns species to empty cell based on relative frequency in population
                // May need to exclude dead trees from being chosen..
                if(dispersal_mode == 1){
                    sp[empty_cell_iter] = sp[cell_rng(generator)]; // Random species
                    continue;
                }
        
            
            
            // LOCAL DISPERSAL
            if(dispersal_mode == 0){

                // Reset neighbors object
                neighbors.reset();
                
                // Find neighbors
                neighbors.updateNeighbors(empty_cell_iter, height,
                                          width, area,
                                          neighbor_radius,
                                          sp, gen, n_alleles_init,
                                          n_sp_init);
       
                // Disperse seeds into cell
                neighbors.disperseSeeds(generator);
                
                // GNDD process
                neighbors.GNDD(gndd_sp);
                
                
                // Print
//                neighbors.printStatus(neighbor_radius, n_sp_init, n_alleles_init);

                // CNDD process
                neighbors.CNDD(cndd_sp);
                
                // Print
               if(verbose) neighbors.printStatus(neighbor_radius, n_sp_init, n_alleles_init);
                
                // Count up seeds
                neighbors.totalSeeds();
                
                // Choose winner and assign to cell
                neighbors.chooseWinner(generator,
                                       empty_cell_iter,
                                       sp, gen,
                                       n_sp_init,
                                       n_alleles_init);
         
            } // End looping over empty cells
            
        } // End local dispersal if
        

        if(step % save_every_n_steps == 0 | step == 1 ){
           // std::cout << "Step:" << step << " | " << calcSpeciesRichness(sp) ;
            
            Summary_step summary(sp, gen, n_sp_init, step);
            
            summary_over_time.push_back(summary);
            
            
            // Write landscape of species to tab delimited .txt file
            std::string species_filename = "landscape_species_step_";
            std::string suffix = ".txt";
            
            // Write to buffer to add leading 0s
            std::stringstream buffer;
            
            buffer << species_filename <<  std::setw(3) << std::setfill('0') << step << suffix;
        
            write_landscape(sp, gen, height, width,  buffer.str());
            
            
     
        }
                       
    } // End step loop

    
    // Writing things to file
    
    // Write landscape of species to tab delimited .txt file
   std::string species_filename = "landscape_species_FINAL.txt";
   write_landscape(sp, gen, height, width, species_filename);

    
    // Write summary to file
    
    write_summary(summary_over_time);

    
    return 0;
}
