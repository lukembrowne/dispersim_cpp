//
//  summary.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/6/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include <stdio.h>
#include <vector>



// Summary class object that calculates summary statistics of the simulation and
// Writes them to file, etc

class Summary {
    
    
public:
    
    // Fields
    int step;
    
    // Species diversity
    int sp_richness;
    
    // Genetic diversity
    std::vector<int> allelic_richness_by_sp;
    float allelic_richness_avg;
    
    //Constructor function
    Summary(std::vector<int>& sp, std::vector<int>& gen, int n_sp_init, int step);
    
    // Member functions
    void calc_sp_richness(std::vector<int>& sp);
    void calc_allelic_richness(int n_sp_init, std::vector<int>& sp, std::vector<int>& gen);
    float calc_avg(std::vector<int> metric); // Maybe need to overload another for float?
    void print();
};


#endif /* summary_hpp */
