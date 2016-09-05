//
//  params.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/13/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "params.hpp"





////////////////////////
////////////////////////
// PARAMETERS

// Lagniappe parameters
Params::Params() :
save_every_n_steps{5},
verbose{false},

// Number of generations
steps{50}, // Each step = 5-10 years with 0.10 % mortality rate
mortality_rate{0.25}, // Proportion of landscape dying per step

// Species parameters
n_sp_init{50},
n_alleles_init{5},
seed_disp_dist{1}, // In units of cells
neighbor_radius{3}, // Miranda = 5; Banitz = 20

seeds_per_adult{500}, // Equal to fecundity..

// NDD parameters
max_cndd{1.0},
min_cndd{1.0}, // Min must be greater numerically than max, but means weaker NDD

max_gndd{1.0},
min_gndd{1.0}, // Min must be greater numerically than max, but means weaker NDD


// Landscape parameters
width{100},
height{100},
area{width * height},

migration_rate{0.0001}, // Immigrant per recruit (~1 in 10,000) is from BCI paper; 1 in 9000 used in Muller Landau 2007

dispersal_mode{0}// 1 == global; 0 == local

{}