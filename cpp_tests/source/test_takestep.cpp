#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "pele/array.h"

#include "mcpele/takestep.h"

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

    virtual void SetUp(){
    seed = 42;
    ndof = 33;
    coor = Array<double>(ndof);
    for (size_t i = 0; i < ndof; ++i){ coor[i] = 4242; }
    reference = coor.copy();
    stepsize = 0.1;
    niterations = 10000;
    f = 0.1;
    }
    virtual void TearDown() {
    /////
    }
};

TEST_F(TakeStepTest, BasicFunctionalityAveragingErasing_OneIteration){
    //one iteration gives expected variation
    mcpele::RandomCoordsDisplacement displ(seed);
    displ.takestep(coor,stepsize);
    for (size_t i = 0; i < ndof; ++i){
    EXPECT_NEAR( reference[i], coor[i], stepsize*0.5 );
    }
}

TEST_F(TakeStepTest, BasicFunctionalityAveragingErasing_NIterations){
    //n iterations give expected variation
    mcpele::RandomCoordsDisplacement displ(seed);
    for (size_t i = 0; i < niterations; ++i){
    displ.takestep(coor,stepsize);
    }
    for (size_t i = 0; i < ndof; ++i){
    EXPECT_NEAR( reference[i], coor[i], f*sqrt(niterations) );
    }
}

TEST_F(TakeStepTest, BasicFunctionalityAveragingErasing_NIterationsReAllocate){
    // n iterations give expected vairation even if step generator is deleted and re-allocated
    for (size_t i = 0; i < niterations; ++i){
    //mcpele::RandomCoordsDisplacement displ(seed); // this constructor re-seeds the rng
    mcpele::RandomCoordsDisplacement displ; // this constructor does not re-seed the rng
    displ.takestep(coor,stepsize);
    }
    for (size_t i = 0; i < ndof; ++i){
    EXPECT_NEAR( reference[i], coor[i], f*sqrt(niterations) );
    }
}
