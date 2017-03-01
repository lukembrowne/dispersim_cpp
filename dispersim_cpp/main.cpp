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
#include <boost/filesystem.hpp>
#include "utils.hpp"
#include "summary.hpp"
#include "neighbors.hpp"
#include "sim.hpp"
#include "params.hpp"
#include <iomanip>
#include <sstream>

#include <boost/random/binomial_distribution.hpp>


/// MAIN

int main(int argc, const char *argv[]) {

    // Parsing command line arguments
    
    // If no arguments, exit
    if(argc == 1)
    {
        std::cout << "No parameter values given. Exiting program... \n";
        exit(1);
    }
    
    // If incorrect number of parameters, exit
    if(argc != 19)
    {
        std::cout << "Incorrect number of parameters given.. Should be 18.. try again. \n";
        exit(1);
    }
    
    
    // Print out command line arguments
     std::cout << argc << "\n";
    for (int count=0; count < argc; ++count) {
        std::cout << count << " " << argv[count] << '\n';

    }
    

    // Initialize parameter list - pass command line arguments
    Params params(argv);
    
    
    // Write params list to file
    boost::filesystem::create_directories("./params_out");

    std::stringstream param_buffer;
    
    param_buffer << "./params_out/params_out_sim_" << std::setw(3) << std::setfill('0') << params.sim_id << ".txt";


    write_params(param_buffer.str(), params);
    
    // Calculate number of dead individuals per step
    // Can't do this within Params.hpp w/o getting compiler errors
    // "cannot appear in a constant-expression"
    int n_dead_per_step = params.mortality_rate * params.area;
    

    
    // Initialize Simulation
    Sim sim(params, n_dead_per_step);

    // Initialize vector of summary object
    std::vector<Summary_step> summary_over_time;
    
    // Initialize Neighbors object
    Neighbors neighbors(params);
    
    
    // Initialize directories
    //boost::filesystem::remove_all("./landscape_out"); // This would erase the folder
    boost::filesystem::create_directories("./landscape_out");
    
    //boost::filesystem::remove_all("./summary_out"); // This would erase the folder
    boost::filesystem::create_directories("./summary_out");
    
    
 
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
        
        
        
        // Loop through empty cells
        for(auto& empty_cell_iter : sim.empty_cell_indices){

            
            /////////////
            // MIGRATION -
            
            if(sim.migration_rng(sim.generator) <= params.migration_rate){
            
                // Replace with random species with random genotype
                
                sim.sp[empty_cell_iter] = sim.species_rng(sim.generator);
                sim.gen[empty_cell_iter] = sim.gen_rng(sim.generator);
            
                
                continue; // Jump to next empty cell
                
            } // End migration if
           

            // Reset neighbors object
            neighbors.reset();
            
            // Find neighbors
            neighbors.updateNeighbors(empty_cell_iter, params, sim);

            // Disperse seeds into cell
            neighbors.disperseSeeds(sim.generator);
            
            // Print
            if(params.verbose){
                std::cout << "\n\n--------------------------- \n Before NDD \n";
                neighbors.printStatus(params);
            }
            
            // NDD process 
            neighbors.NDD(sim.gndd_sp, sim.cndd_sp);
            
            // Print
            if(params.verbose){
                std::cout << "\n\n--------------------------- \n After NDD \n";
                // Seeds by species is not updated here
                neighbors.printStatus(params);
            }
            
            // Count up seeds
            neighbors.totalSeeds();
            
            // Choose winner and assign to cell
            neighbors.chooseWinner(empty_cell_iter, params, sim);
     
        } // End looping over empty cells


        // Saving summaries
        if(step % params.save_every_n_steps == 0 | step == 1 ){
            
            Summary_step summary(sim, params, step);
            
            summary_over_time.push_back(summary);
            
            
//            // Write landscape of species to tab delimited .txt file
//            std::string species_filename = "./landscape_out/landscape_species_";
//            std::string suffix = ".txt";
//            
//            // Write to buffer to add leading 0s
//            std::stringstream buffer;
//            
//            // Adding sim_id and step number to file name, setw and setfill add leading for sorting
//            buffer << species_filename << "sim_" << std::setw(3) << std::setfill('0') << params.sim_id << "_step_" << std::setw(3) << std::setfill('0') << step << suffix;
//        
//            write_species_landscape(buffer.str(), sim, params);
            
        }
                       
    } // End step loop
    
    // Write landscape of species to tab delimited .txt file
    std::string species_filename = "./landscape_out/landscape_species_";
    std::string suffix = "_FINAL.txt";
    
    // Write to buffer to add leading 0s
    std::stringstream buffer;
    
    // Adding sim_id and step number to file name, setw and setfill add leading for sorting
    buffer << species_filename << "sim_" << std::setw(3) << std::setfill('0') << params.sim_id << suffix;
    
    write_species_landscape(buffer.str(), sim, params);

   // write_genotype_landscape(sim, params);
    
    // Write summary overall and summary by sp to file
    write_summary(summary_over_time, params, sim);

    
    return 0;
}
