#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "mcpele/mc.h"
#include "mcpele/random_coords_displacement.h"
#include "mcpele/particle_pair_swap.h"

using pele::Array;

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

class TakeStepTest: public ::testing::Test{
public:
    typedef pele::Array<double> arr_t;

    size_t seed;
    size_t ndof;
    arr_t coor;
    arr_t reference;
    double stepsize;
    size_t niterations;
    double f;
    mcpele::MC* mc;

    virtual void SetUp(){
        seed = 42;
        ndof = 33;
        coor = Array<double>(ndof);
        for (size_t i = 0; i < ndof; ++i) {
            coor[i] = 4242 + i;
        }
        reference = coor.copy();
        stepsize = 0.1;
        niterations = 10000;
        f = 0.1;
        mc = NULL;
    }
};

TEST_F(TakeStepTest, BasicFunctionalityAveragingErasing_OneIteration){
    //one iteration gives expected variation
    mcpele::RandomCoordsDisplacement displ(seed, stepsize);
    displ.displace(coor, mc);
    for (size_t i = 0; i < ndof; ++i){
        EXPECT_NEAR( reference[i], coor[i], stepsize*0.5 );
    }
}

TEST_F(TakeStepTest, BasicFunctionalityAveragingErasing_NIterations){
    //n iterations give expected variation
    mcpele::RandomCoordsDisplacement displ(seed, stepsize);
    for (size_t i = 0; i < niterations; ++i){
    displ.displace(coor, mc);
    }
    for (size_t i = 0; i < ndof; ++i){
    EXPECT_NEAR( reference[i], coor[i], f*sqrt(niterations) );
    }
}

TEST_F(TakeStepTest, BasicFunctionalityAveragingErasing_NIterationsReAllocate){
    // n iterations give expected vairation even if step generator is deleted and re-allocated
    mcpele::RandomCoordsDisplacement displ(seed, stepsize); // this constructor re-seeds the rng
    for (size_t i = 0; i < niterations; ++i){
        displ.displace(coor, mc);
    }
    for (size_t i = 0; i < ndof; ++i){
        EXPECT_NEAR( reference[i], coor[i], f*sqrt(niterations) );
    }
}

TEST_F(TakeStepTest, PairSwapWorks){
    const size_t box_dimension = 3;
    const size_t nr_particles = ndof / box_dimension;
    mcpele::ParticlePairSwap swap(42, nr_particles, 1);
    auto coor1 = coor.copy();
    auto coor2 = coor.copy();
    const size_t a = 1;
    const size_t b = 4;
    swap.swap_coordinates(a, b, coor1);
    for (size_t i = 0; i < ndof; ++i) {
        if (i >= a * box_dimension && i < (a + 1) * box_dimension) {
            EXPECT_DOUBLE_EQ(coor1[i], coor[(i + b * box_dimension) - a * box_dimension]);
        }
        else if (i >= b * box_dimension && i < (b + 1) * box_dimension) {
            EXPECT_DOUBLE_EQ(coor1[i], coor[(i + a * box_dimension) - b * box_dimension]);
        }
        else {
            EXPECT_DOUBLE_EQ(coor1[i], coor[i]);
        }
    }
    swap.swap_coordinates(8, 8, coor2);
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_DOUBLE_EQ(coor2[i], coor[i]);
    }
}
