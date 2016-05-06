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
#include "summary.hpp"

// Holds information used in NN and GNDD and CNDD calculation
struct key {
    int sp;
    int gen;
    int gen_1d_index;
};

bool SortKey( const key& elem1, const key& elem2 )
{
    return elem1.gen_1d_index < elem2.gen_1d_index;
}

bool CompKey(const key& elem1, const key& elem2){
    return elem1.gen_1d_index == elem2.gen_1d_index;
}


/// MAIN

int main(int argc, const char * argv[]) {

////////////////////////
////////////////////////
// PARAMETERS
    
    // Number of generations
    int steps = 100; // Each step = 5-10 years with 0.10 % mortality rate
    float mortality_rate = 0.1; // Proportion of landscape dying per step
    
    // Species parameters
    int n_sp_init = 50;
    int n_alleles_init = 15;
    float seed_disp_dist = 5; // In units of cells
    int seeds_per_adult = 500; // Equal to fecundity..
    
    // NDD parameters
    float max_cndd = 0.1; // Lowering this produced more clustered patterns... WHY??
    float min_cndd = 0.9; // Min must be greater numerically than max, but means weaker NDD
   
    float max_gndd = 0.5;
    float min_gndd = 0.5; // Min must be greater numerically than max, but means weaker NDD
    
    
    // Landscape parameters
    int width  = 100;
    int height = 100;
    int area = width * height;
    
    float migration_rate = 0.0001; // Immigrant per recruit (~1 in 10,000) is from BCI paper
    
    int dispersal_mode = 0; // 1 == global; 0 == local
    
    int n_dead_per_step = mortality_rate * area;
    int empty_cell_indices[n_dead_per_step]; // Initialize
    
    
    // Lagniappe parameters
    int save_every_n_steps = 10;
    
    std::vector<Summary_step> summary_over_time;
    
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
        }

    
    // Add genotypes randomly to species
    
        // Initialize genotype RNG
        std::uniform_int_distribution<int> gen_rng(0, (n_alleles_init-1));

        std::vector<int> gen(area); // Init vector that will hold species
    
        for (auto& iter : gen){
            iter = gen_rng(generator);
        }
    
    // Assign NDD values to species
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
    std::vector<float> seeds_by_sp(n_sp_init, 0.0); // Initialize to 0
    std::vector<float> seeds_by_gen(n_sp_init * n_alleles_init, 0.0);
    float seeds_total = {0.0};
    
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
                
                // Initialize variables
                std::fill(seeds_by_sp.begin(), seeds_by_sp.end(), 0);
                seeds_total = 0; // Reset seeds total to 0
                float seeds_to_add;
                std::vector<key> nn_keys(8);
                
                
                // Find surrounding 8 neighbors
                // Returns position in 1d array
                 neighbors = findNN(empty_cell_iter, height, width, area);
                
                
                // Disperse seeds from 8 neighbors into empty cell

                // Loop over neighbors (1d array index for 8 NN)
                i = 0;
                for(auto iter : neighbors){
                    
                    // Generate species key
                    // A vector of length 8 that holds species names
                    nn_keys[i].sp = sp[iter]; // SPecies id
                    nn_keys[i].gen = gen[iter]; // Genotype id
                    nn_keys[i].gen_1d_index = sp[iter] * n_alleles_init + gen[iter]; // Position in 1d array

                    // Vector holding number of seeds contributed by each NN
                    seeds_to_add = seed_rng(generator);
                    
                    seeds_by_gen[nn_keys[i].gen_1d_index] += seeds_to_add; // Add seeds for that specific genotype
                    
//                    // Printing for error checking
//                    std::cout << "Species: " << nn_keys[i].sp << " | Gen:" << nn_keys[i].gen << " | Gen_index: " << nn_keys[i].gen_1d_index << " | Seeds: " << seeds_by_gen[nn_keys[i].gen_1d_index] << "\n" ;

                     i++;
                    
                }
                
//                std::cout << "\n| ------------- | \n";
                
                
                // Need to sort and remove duplicates of species - genotype pairs to avoid double counting
                // Relies on comparison functions in 'key' structure
                
                    // Sort based on gen_1d_index
                     std::sort(nn_keys.begin(), nn_keys.end(), SortKey);
                    // Delete duplicates based on gen_1d_index
                    nn_keys.erase(std::unique(nn_keys.begin(), nn_keys.end(), CompKey), nn_keys.end());
                
//                
//                // Printing for error checking
//                    std::cout << "\n| ----- Post sort and delete dupes -------- | \n";
//                    i = 0;
//                    for(auto iter : nn_keys){
//                    // Printing for error checking
//                    std::cout << "Species: " << nn_keys[i].sp << " | Gen:" << nn_keys[i].gen << " | Gen_index: " << nn_keys[i].gen_1d_index << " | Seeds: " << seeds_by_gen[nn_keys[i].gen_1d_index] << "\n" ;
//                    i++;
//                    }
//                std::cout << "\n| ------------- | \n";
                
                
                
                
                // GNDD
                // Reduce densities based on Genotype - GNDD
                
                std::vector<float> n_gens_per_sp(n_sp_init, 0.0); // Move outside loop?
                
//                std::cout << "\n| ----- After GNDD -------- | \n";

                for(int i = 0; i < nn_keys.size(); i++){
                    
                    // Following equation of Harms et al. 2000
                    // Where ## of recruits is a function of beta exponent and log + 1 density of seeds
                    // Here, converting predicted log +1 density of recruits back to actual number of recruits with std::exp
                    // Adding + 1 to recruit density makes it so theres always at least 2.7 possible recruits
                    // Assume intercept (alpha) = 0
                    
                    seeds_by_gen[nn_keys[i].gen_1d_index] = std::exp(gndd_sp[nn_keys[i].sp] *
                                                                                   std::log(seeds_by_gen[nn_keys[i].gen_1d_index]));
                    
                    seeds_by_sp[nn_keys[i].sp] += seeds_by_gen[nn_keys[i].gen_1d_index]; // Add to total seeds
                    
                    n_gens_per_sp[nn_keys[i].sp] += 1;
                    
                    
//                    // Printing for error checking
//                    std::cout << "Species: " << nn_keys[i].sp << " | Gen:" << nn_keys[i].gen << " | Gen_index: " << nn_keys[i].gen_1d_index << " | Seeds: " << seeds_by_gen[nn_keys[i].gen_1d_index] << "\n" ;
                }
                
                
                
                /// CNDD
                // Reduce densities overall based on number of conspecifics, without regard to genotype
                
//                std::cout << "\n| ----- After CNDD -------- | \n";
               
                
                for(int i = 0; i < nn_keys.size(); i++){ // Loop across species and gen keys
                    
                    seeds_by_gen[nn_keys[i].gen_1d_index] -= std::exp((cndd_sp[nn_keys[i].sp] *
                                                                    std::log(seeds_by_sp[nn_keys[i].sp]))) / n_gens_per_sp[nn_keys[i].sp];
                    
                    seeds_total += seeds_by_gen[nn_keys[i].gen_1d_index];
                    
                 
                }
                


                
                // Assign probs based on relative frequency
                    std::vector<float> probabilities(nn_keys.size());
                    i = 0;
                    
                    for(auto& iter : nn_keys){
                        
                        probabilities[i] = seeds_by_gen[nn_keys[i].gen_1d_index]/seeds_total;
                        
                        
//                        // Printing for error checking
//                        std::cout << "Species: " << nn_keys[i].sp << " | Gen:" << nn_keys[i].gen << " | Gen_index: " << nn_keys[i].gen_1d_index << " | Seeds: " << seeds_by_gen[nn_keys[i].gen_1d_index] << " Probability: " << probabilities[i] << "\n" ;
//                        
                        // Reset seeds count
                        seeds_by_gen[nn_keys[i].gen_1d_index] = 0;
                        
                        i++;
                    }
                    
                
                
                    // Choose species to establish in empty cell
                    // Choose a winner / seed to establish based on relative frequency and weighted probability
                       boost::random::discrete_distribution<> seed_winner_rng(probabilities);
                
                    int winner_index; // Move up and out of loop?
                    winner_index = seed_winner_rng(generator);
                
                // Reassign Species ID and genotype for 'winner' then continue looping through empty cells
                
                    sp[empty_cell_iter] = nn_keys[winner_index].sp; // Reassign species
                    gen[empty_cell_iter] = nn_keys[winner_index].gen;

                
            } // End looping over empty cells
            
        } // End local dispersal if
        

        if(step % save_every_n_steps == 0 | step == 1 ){
           // std::cout << "Step:" << step << " | " << calcSpeciesRichness(sp) ;
            
            Summary_step summary(sp, gen, n_sp_init, step);
            
            summary_over_time.push_back(summary);
     
        }
                       
    } // End step loop

    
    // Writing things to file
    
    // Write landscape of species to tab delimited .txt file
   // write_landscape(sp, gen, height, width);

    
    // Write summary to file
    
    write_summary(summary_over_time);

    
    return 0;
}
