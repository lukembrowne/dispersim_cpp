//
//  params.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/13/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef params_hpp
#define params_hpp

#include <stdio.h>


////////////////////////
////////////////////////
// PARAMETERS

class Params {
    
public:
    // Lagniappe parameters
    int save_every_n_steps;
    bool verbose;

    // Number of generations
    int steps; // Each step = 5-10 years with 0.10 % mortality rate
    float mortality_rate; // Proportion of landscape dying per step

    // Species parameters
    int n_sp_init;
    int n_alleles_init;
    float seed_disp_dist; // In units of cells
    int neighbor_radius; // Miranda = 5; Banitz = 20

    int seeds_per_adult; // Equal to fecundity..

    // NDD parameters
    float mean_cndd; // Closer to 1 = weaker NDD
    float range_cndd;

    float mean_gndd; // Closer to 1 = weaker NDD
    float range_gndd;


    // Landscape parameters
    int width;
    int height;
    int area;

    float migration_rate; // Immigrant per recruit (~1 in 10,000) is from BCI paper; 1 in 9000 used in Muller Landau 2007

    int dispersal_mode ; // 1 == global; 0 == local
    
    int sim_id ; // Unique ID for the sim
        
    Params(const char * argv[]);

};

#endif /* params_hpp */
