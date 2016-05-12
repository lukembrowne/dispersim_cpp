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

void write_landscape(std::vector<int>& sp, std::vector<int>& gen,
                     int height, int width,
                     std::string species_filename);

void write_summary(std::vector<Summary_step>& summary_over_time);


#endif /* utils_hpp */
