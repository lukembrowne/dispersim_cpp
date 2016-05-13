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
#include "params.hpp"
#include <iomanip>
#include <sstream>



/// MAIN

int main(int argc, const char * argv[]) {


    // Initialize parameter list
    Params params;
    
    // Initialize Simulation
    Sim sim(params);

    // Initialize vector of summary object
    std::vector<Summary_step> summary_over_time;
    
    
    // Initialize Neighbors object
    Neighbors neighbors(params);
    
    neighbors.initSeedRNG(params);
    

////////////////////////
////////////////////////
    
// Begin simulation
    
    std::cout << "\n\n";
    
    for (int step = 1; step <= params.steps; step++){
        
        
    // Death of adults
    // A set proportion dies every step, dependent on mortality rate
        
        // Choose which cells will die
        for(int i = 0; i < sim.empty_cell_indices.size(); i++) {
            sim.empty_cell_indices[i] = sim.cell_rng(sim.generator);
        }
        
        // Loop through empty cells and assign new species
 
        for(auto& empty_cell_iter : sim.empty_cell_indices){

            // MIGRATION -
            
            if(sim.migration_rng(sim.generator) <= params.migration_rate){
            
                // Replace with random species with random genotype
                
                sim.sp[empty_cell_iter] = sim.species_rng(sim.generator);
                sim.gen[empty_cell_iter] = sim.gen_rng(sim.generator);
            
                
                continue; // Jump to next empty cell
            }
            
            // GLOBAL DISPERAL
            // Need to change so that some proportion of adults contribute seedlings, but these are still evaluated for NDD effects.. just that dispersal is not limited by proximity - randomly choose 8 adults to contribute seeds
            
                // Assigns species to empty cell based on relative frequency in population
                // May need to exclude dead trees from being chosen..
                if(params.dispersal_mode == 1){
                    sim.sp[empty_cell_iter] = sim.sp[sim.cell_rng(sim.generator)]; // Random species
                    continue;
                }
        
            
            
            // LOCAL DISPERSAL
            if(params.dispersal_mode == 0){

                // Reset neighbors object
                neighbors.reset();
                
                // Find neighbors
                neighbors.updateNeighbors(empty_cell_iter, params, sim);
       
                // Disperse seeds into cell
                neighbors.disperseSeeds(sim.generator);
                
                // GNDD process
                neighbors.GNDD(sim.gndd_sp);
                
                // Print
//                neighbors.printStatus(neighbor_radius, n_sp_init, n_alleles_init);

                // CNDD process
                neighbors.CNDD(sim.cndd_sp);
                
                // Print
               if(params.verbose) neighbors.printStatus(params);
                
                // Count up seeds
                neighbors.totalSeeds();
                
                // Choose winner and assign to cell
                neighbors.chooseWinner(empty_cell_iter, params, sim);
         
            } // End looping over empty cells
            
        } // End local dispersal if
        

        if(step % params.save_every_n_steps == 0 | step == 1 ){
           // std::cout << "Step:" << step << " | " << calcSpeciesRichness(sp) ;
            
            Summary_step summary(sim, params, step);
            
            summary_over_time.push_back(summary);
            
            
            // Write landscape of species to tab delimited .txt file
            std::string species_filename = "landscape_species_step_";
            std::string suffix = ".txt";
            
            // Write to buffer to add leading 0s
            std::stringstream buffer;
            
            buffer << species_filename <<  std::setw(3) << std::setfill('0') << step << suffix;
        
            write_landscape(buffer.str(), sim, params);
            
            
     
        }
                       
    } // End step loop

    
    // Writing things to file
    
    // Write landscape of species to tab delimited .txt file
   std::string species_filename = "landscape_species_FINAL.txt";
   write_landscape(species_filename, sim, params);

    
    // Write summary to file
    
    write_summary(summary_over_time);

    
    return 0;
}
