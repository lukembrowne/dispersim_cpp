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



// Initialize Landscape - fill with 0s
void Landscape::initLandscape(){
 
    std::cout << "Initializing map with.. " << width * height << " cells ... \n";
    
    for(int i = 0; i < (width * height); i++){
    
        species_list.insert(std::make_pair(i, 0));
    }
}


// Add species to landscape

void Landscape::initSpecies(){
    
    
    // Good info on generating random numbers in C++
    // http://diego.assencio.com/?index=6890b8c50169ef45b74db135063c227c
    
    std::random_device device;
    std::mt19937 generator(device());
    std::uniform_int_distribution<int> rand_species(0, max_species);
    
    for (auto& iter : species_list){
        
        iter.second = rand_species(generator);
        
    //    std::cout << iter.second << ' ';
      //  std::cout << '\n';
    }
}



// Get species

int Landscape::getSpecies(int x, int y){
    
    return species_list[x + (y * width)];

}