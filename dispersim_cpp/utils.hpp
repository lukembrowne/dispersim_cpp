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


float neg_expo_discrete(float alpha, float dij, int R);

void writeLandscape(std::vector<int>& sp, std::vector<int>& gen, int height, int width);


#endif /* utils_hpp */
