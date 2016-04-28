//
//  landscape.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/28/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef landscape_hpp
#define landscape_hpp

#include <stdio.h>
#include <map>


class Landscape {
public:
    int width;
    int height;
    int max_species;
    std::map <int, int> species_list;
    
    //  Create map data structure
    void initLandscape();
    
    //  Add species to landscape
    void initSpecies();
    
    // Find species given x,y coords
    int getSpecies(int x, int y);
    
};


#endif /* landscape_hpp */
