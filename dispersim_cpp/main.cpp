//
//  main.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/28/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include <iostream>
#include <random>
#include "landscape.hpp"

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    
    
    // Initialize Landscape
    
    Landscape landscape = {50,50, 100};
    
    std::cout << landscape.height << "\n";
  
    landscape.initLandscape();

    std::cout << "Before adding species.. " << landscape.getSpecies(0, 0) << "\n";
    
    
    landscape.initSpecies();
    
    std::cout << "Adding adding species.. " << landscape.getSpecies(0, 0) << "\n";
    

    
    // Playground
    
    
    return 0;
}
