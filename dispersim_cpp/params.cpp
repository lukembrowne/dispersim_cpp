//
//  params.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/13/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "params.hpp"
#include <sstream> // for std::stringstream





////////////////////////
////////////////////////
// PARAMETERS

// Lagniappe parameters
Params::Params(const char *argv[]) {

   // Parameters must be entered in the command line in this order
   // With first argument starting at index 1
  
    
    // 1 - How often to save out data save_every_n_steps

        // Convert to int
        std::stringstream convert(argv[1]); // set up a stringstream variable named convert, initialized with the input from argv[1]
        convert >> save_every_n_steps;
    
    
    // 2 - Whether printing is verbose or not
    
    int verbose_temp = 2;
    std::stringstream convert2(argv[2]);
    convert2 >> verbose_temp;
    
    if(verbose_temp == 1){
        verbose =  true;
    } else {
        verbose = false;
    }
    
    
    
    // 3 -  Number of steps
    
        // Each step = 5-10 years with 0.10 % mortality rate
        std::stringstream convert3(argv[3]);
        convert3 >> steps;
    
    // 4 - Mortality rate
    
        // Proportion of landscape dying per step
        std::stringstream convert4(argv[4]);
        convert4 >> mortality_rate;

    
    // 5 - Number of initial species
    
        std::stringstream convert5(argv[5]);
        convert5 >> n_sp_init;
    
    
    // 6 - Number of initial alleles per species
    
        std::stringstream convert6(argv[6]);
        convert6 >> n_alleles_init;


    // 7 - Seed dispersal distance
    
        // In units of cells
        std::stringstream convert7(argv[7]);
        convert7 >> seed_disp_dist;
    
    
    // 8 - Neighbor radius for dispersal
    
        // In units of cells
        // Miranda = 5; Banitz = 20
        std::stringstream convert8(argv[8]);
        convert8 >> neighbor_radius;
    
    // 9 - Seeds per adult
    
        // Equal to fecundity..
        std::stringstream convert9(argv[9]);
        convert9 >> seeds_per_adult;
    
    // 10, 11 - CNDD parameters

        // Min must be greater numerically than max, but means weaker NDD
        std::stringstream convert10(argv[10]);
        convert10 >> max_cndd;
    
        std::stringstream convert11(argv[11]);
        convert11 >> min_cndd;
    
    // 12, 13 - GNDD parameters
    
        // Min must be greater numerically than max, but means weaker NDD
        std::stringstream convert12(argv[12]);
        convert12 >> max_gndd;
        
        std::stringstream convert13(argv[13]);
        convert13 >> min_gndd;

    
    // 14, 15  - Landscape parameters - height, width, area
    
        std::stringstream convert14(argv[14]);
        convert14 >> width;
        
        std::stringstream convert15(argv[15]);
        convert15 >> height;
    
        area = width * height;
    
    // 16 - Migration rate
    
        // Immigrant per recruit (~1 in 10,000) is from BCI paper; 1 in 9000 used in Muller Landau 2007
        std::stringstream convert16(argv[16]);
        convert16 >> migration_rate;
    
    // 17 - Dispersal mode, 1 == global; 0 == local
    
        std::stringstream convert17(argv[17]);
        convert17 >> dispersal_mode;
}