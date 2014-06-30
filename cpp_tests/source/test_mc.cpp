#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include "pele/harmonic.h"

#include "takestep.h"

#include "mc.h"

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

class TestMC: public ::testing::Test{
public:

    typedef pele::Array<double> arr_t;
    size_t boxdim;
    size_t nparticles;
    size_t ndof;
    arr_t origin;
    arr_t x;
    double stepsize;
    double k;
    pele::Harmonic* potential;
    size_t max_iter;

    virtual void SetUp(){
	boxdim = 3;
	nparticles = 1e2;
	ndof = boxdim*nparticles;
	origin.resize(ndof);
	std::fill(origin.data(),origin.data()+ndof,0);
	x.resize(ndof);
	std::fill(x.data(),x.data()+ndof,0);
	stepsize = 1e-2;
	k = 400;
	potential = new pele::Harmonic(origin, k, boxdim);
	max_iter = 1e3;
    }

    virtual void TearDown() {
	delete potential;
    }
};

TEST_F(TestMC, BasicFunctionalityAddingModules){
    mcpele::MC mc(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc.get_potential_ptr())->get_k() );
    mcpele::RandomCoordsDisplacement sampler_uniform(42);
    mc.set_takestep(&sampler_uniform);
    EXPECT_TRUE( mc.take_step_specified() );
    mc.run(max_iter);
}
