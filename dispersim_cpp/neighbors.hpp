//
//  neighbors.hpp
//  dispersim_cpp
//
//  Created by Luke Browne on 5/6/16.
//  Copyright © 2016 Luke Browne. All rights reserved.
//

#ifndef neighbors_hpp
#define neighbors_hpp

#include <stdio.h>
#include <random>
#include <vector>


// Create a class that holds information about neighbors around an empty cell
// Holds information on:

// Who the neighbors are in a certain Radius of cells R
// Species of neighbors
// Genotypes of neighbors
// probability of dispersal (neg exponential RNG) for each neighbor


//This does bring us to an important topic though: What should you store in a vector? Built-ins like integers can be stored easily, of course. But what about objects of class type? Suppose you are programming a particle system and have a particle class, and want to store a long list of particles, what is the best way to do it? One way is to store a vector<particle>, another is to store a vector<particle*>
//— that is, don’t store the particles themselves directly in the vector, but store pointers to particles.

//How do you choose which? Generally, the larger and bulkier the object is, the more likely it is you want to store pointers to it, rather than the object itself. Storing a vector<int*> would be very inefficient, since the pointers would be as large or larger than the integers and you’d have to have the overhead of the memory allocations too. But for a large object, like Frogatto’s custom_object class, a vector<custom_object*> is probably what we want. Note that to store an object directly, it must be copyable, i.e. have accessible copy constructors and assignment operators.
//
//Note also that if you store a vector of pointers, the vector will not manage the memory pointed to by the pointers. If you want the object’s memory to be managed for you, you could use a vector
//<boost::shared_ptr<particle> > to have a vector of ‘smart pointers’ that manage the memory they point to.


//It means something like this:

//std::vector<Movie *> movies;
//Then you add to the vector as you read lines:
//
//movies.push_back(new Movie(...));



// Create a neighbors class that is re-used during demographic stage that holds
// Species ID of neighbors
// Genotypes of neighbors
// Gen 1d index
// Seeds by genotype
// RNG for dispersal for each neighbor that's called each iteration
// size will be dependent on R = or how many cells in radius are searched for potential dispersal

// Can have a printing function that prints out ## of seeds per species per genotype, etc


class Neighbors {
    
public:
    
    // Fields
    int n_neighbors;
    std::vector<int> nn_sp;
    std::vector<int> nn_index;
    std::vector<int> nn_gen;
    std::vector<int> nn_gen_1d_index;
    std::vector<bool> nn_gen_1d_index_dupe; // TRUE if is a duplicate
    std::vector<float> probabilities; // Recruitment probabilities
    
    std::vector<std::binomial_distribution<int>> seed_rng;
    std::vector<float> seeds_by_sp; // Initialize to 0
    std::vector<float> dead_seeds_sp;
    std::vector<float> seeds_by_gen;
    float seeds_total;    

    //Constructor function
    Neighbors(int radius, int n_sp_init, int n_alleles_init);
    
    //Initialize seed RNGS based on distance
    void initSeedRNG(int neighbor_radius, int seed_disp_dist,
                     int seeds_per_adult);
    
    // Find 1d indices of neighhbors
    void getNeighborIndex(int focal_cell, int height,
                          int width, int area,
                          int neighbor_radius);
    
    // Updates nn_sp, nn_gen - calls getNeighborIndex
    void updateNeighbors(int focal_cell, int height,
                         int width, int area,
                         int neighbor_radius,
                         std::vector<int>& sp,
                         std::vector<int>& gen,
                         int n_alleles_init,
                         int n_sp_init);
    
    // Generate numeber of seeds contributed by each neighbor
    // Based on distance and neg exponential dispersal kernel
    // Updates seeds_by_sp and seeds_by_gen
    void disperseSeeds(std::mt19937& generator);
    
    
    // Genotype dependent density dependence
    void GNDD(std::vector<float>& gndd_sp);
    
    // Conspecific NDD
    void CNDD(std::vector<float>& cndd_sp);
    
    // Count total number of seeds
    void totalSeeds();
    
    // Choose recruit based on relative frequency in cell after GNDD and CNDD
    void chooseWinner(std::mt19937& generator,
                      int focal_cell,
                      std::vector<int>& sp,
                      std::vector<int>& gen,
                      int n_sp_init,
                      int n_alleles_init);
    
    // Print status
    void printStatus(int neighbor_radius, int n_sp_init, int n_alleles_init);
    
    // Reset - set values back to 0
    void reset();
    
};



#endif /* neighbors_hpp */
