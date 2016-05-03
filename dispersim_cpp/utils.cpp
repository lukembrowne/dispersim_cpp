//
//  utils.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/29/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "utils.hpp"
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>


// Returns probability of dispersal for negative exponential, discretized
// Following Banitz et al. 2008
// alpha = 1 / mean dispersal distance
// dij = distance between source and endpoint
// R = Radius of cells from source cell
float neg_expo_discrete(float alpha, float dij, int R){
    
    float numerator = exp(-alpha * dij);
    
    float denominator{0.0};
    
    for(int i = -R; i < R; i++){
         for(int j = -R; j < R; j++){
             denominator += exp(-alpha * sqrt(pow(i, 2) + pow(j, 2)));
         }
    }
    
    float prob;
    prob = numerator / denominator;
    
    return(prob);
    
}




// Calculate species richness

int calcSpeciesRichness(std::vector<int>& sp){
    
    // Copy species vector
    
    
    std::vector<int> sp_copy(sp);
    
    
    // Sort vector before eliminating duplicates
    std::sort(sp_copy.begin(), sp_copy.end());
    
    sp_copy.erase(std::unique(sp_copy.begin(), sp_copy.end()), sp_copy.end());
    
    int richness = sp_copy.size();
    
    std::cout << "Richness.. " << richness << " .. \n";
    
    
    return richness;
}


// Write landscape to file
// Export as tab delimited file

void writeLandscape(std::vector<int>& sp, int height, int width){
    
    std::ofstream out("landscape_species.txt");
    
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            out << sp[i*width + j] << "\t";
        } // End column loop
        
        out << "\n";
        
    } // End row loop
    
    out.close();
    
    std::cout << "Successfully wrote to file... \n";
    
}

