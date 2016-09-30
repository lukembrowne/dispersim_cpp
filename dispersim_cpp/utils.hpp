//
//  utils.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/29/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <vector>
#include "summary.hpp"


float neg_expo_discrete(float alpha, float dij, int R);

void write_landscape(std::string species_filename, Sim& sim, Params& params);

void write_summary(std::vector<Summary_step>& summary_over_time, Params& params);

void write_params(std::string params_filename, Params& params);


#endif /* utils_hpp */
