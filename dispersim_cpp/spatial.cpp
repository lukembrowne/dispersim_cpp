//
//  spatial.cpp
//  dispersim_cpp
//
//  Created by Luke Browne on 4/29/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#include "spatial.hpp"
#include <iostream>
#include <assert.h>

// Find 8 nearest neighbors of a cell
void findNN(int focal_cell, int height, int width, int area){
    
    std::vector<int> NN(8);
    
    // Returns a vector with the index to 8 closest neighbors
    // In this order: N NE E SE S SW W NW
    // With toroidal boundary wrapping..
    // https://community.oracle.com/thread/1257021?start=0&tstart=0
    // http://stackoverflow.com/questions/3775905/programming-logic-how-to-check-for-neighbors-in-a-grid
    
    assert(focal_cell >= 0 && focal_cell < area);

    // N
    NN[0] = focal_cell - width;
    if(NN[0] < 0) NN[0] = area + NN[0]; // focal_cell on top row
    
    // S
    NN[4] = focal_cell + width;
    if(NN[4] > (area - 1)) NN[4] = NN[4] - area; // focal_cell on btm row

    // E
    NN[2] = focal_cell + 1;
    
    // W
    NN[6] = focal_cell - 1;
    
    // NE
    NN[1] = NN[0] + 1;
    
    // NW
    NN[7] = NN[0] - 1;
    
    // SE
    NN[3] = NN[4] + 1;
    
    // SW
    NN[5] = NN[4] - 1;
    

    // Toroidal boundary checks
    
    // If on left column
    if((focal_cell) % width == 0 ){// If on left column
        
        // W
        NN[6] = focal_cell + width - 1;
        
        // SW
        NN[5] = NN[6] + width;
        
        // NW
        NN[7] = NN[6] - width;
    }
    
    // If on right column
    if((focal_cell + 1 ) % width == 0 ){
        // E
        NN[2] = focal_cell - width + 1;
        
        // SE
        NN[3] = NN[2] + width;
        
        // NE
        NN[1] = NN[2] - width;
    }
    
    // If on top left corner
    if(focal_cell == 0) {
        NN[7] = area - 1; // NW
    }
    
    // If bottom right corner
    if(focal_cell == (area - 1)){
        NN[3] = 0; // SE
    }

    // If top right corner
    if(focal_cell == (width - 1)){
        NN[1] = (height - 1) * width; // NE
    }
    
    // If bottom left corner
    if(focal_cell == (height - 1) * width){
        NN[5] = width - 1; // SW
    }
    
    
//    // Print 8 NN
//    for(auto& iter : NN){
//        std::cout << iter << " | ";
//    }
}
