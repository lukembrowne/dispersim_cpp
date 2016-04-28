//
//  landscape.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/28/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "landscape.hpp"
#include <iostream>
#include <random>



// Initialize Landscape based on height and width - fill with 0s
// Stored as a 1d map data structure with key - value pairs
// Length of map is = length * width
void Landscape::initLandscape(){
 
    std::cout << "Initializing map with.. " << width * height << " cells ... \n";
    
    for(int i = 0; i < (width * height); i++){
    
        species_list.insert(std::make_pair(i, 0));
    }
}





// Add species to landscape
// Added uniformly randomly across the landscape
void Landscape::initSpecies(){
    
    
    // Good info on generating random numbers in C++
    // http://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c
    
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<int> rand_species(1, max_species);
    
    // Use this code to loop over map list
    for (auto& iter : species_list){
        
        // iter.second refers to value in key-value map pair
        // iter.first would refer to first value
        iter.second = rand_species(generator);

    }
}



// Get species

int Landscape::getSpecies(int x, int y){
    
    return species_list[x + (y * width)];

}



// Calculate species richness - maybe move to it's own file later on

int Landscape::calcSpeciesRichness(){
    
    // Make a vector that will store species ID
    std::vector<int> species_vector;
    species_vector.reserve(species_list.size()); // Reserve memory
    
    // Loop over map list and save species into vector
    for (auto& iter : species_list){
        species_vector.push_back(iter.second);
    }
    
    // Sort vector before eliminating duplicates
    std::sort(species_vector.begin(), species_vector.end());
  
    species_vector.erase(std::unique(species_vector.begin(), species_vector.end()), species_vector.end());
    
//    for( int i = 0; i < 20; i++){
//        std::cout << species_vector[i] << "\n";
//    }

    int richness = species_vector.size();
    
    std::cout << "Richness.. " << richness << " .. \n";
    
    
    return richness;
}
