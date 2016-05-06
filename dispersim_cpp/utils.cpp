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





// Write landscape for both species and genotypes to file
// Export as tab delimited file

void writeLandscape(std::vector<int>& sp, std::vector<int>& gen, int height, int width){
    
    std::ofstream out_species("landscape_species.txt");
    std::ofstream out_genotypes("landscape_genotypes.txt");
    
    for(int i = 0; i < width; i++){
        for(int j = 0; j < height; j++){
            out_species << sp[i*width + j] << "\t";
            out_genotypes<< gen[i*width + j] << "\t";

        } // End column loop
        
        out_species << "\n";
        out_genotypes << "\n";
        
    } // End row loop
    
    out_species.close();
    out_genotypes.close();
    
    std::cout << "Successfully wrote to file... \n";
    
}

