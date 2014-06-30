#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include "pele/harmonic.h"

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
    /*

    size_t ndof;
    size_t nsteps;
    size_t ntot;
    arr_t displ_uniform;
    arr_t displ_gaussian;
    arr_t baseline;
    double ss;
    */
    virtual void SetUp(){
	boxdim = 3;
	nparticles = 1000;
	ndof = boxdim*nparticles;
	origin.resize(ndof);
	std::fill(origin.data(),origin.data()+ndof,0);
	x.resize(ndof);
	std::fill(x.data(),x.data()+ndof,0);
	stepsize = 1e-2;
	k = 400;
	potential = new pele::Harmonic(origin, k, boxdim);
	max_iter = 1e5;
	/*
	displ_uniform.resize(ndof);
	displ_gaussian.resize(ndof);
	baseline.resize(ndof);
	for (size_t i = 0; i < ndof; ++i){
	    displ_uniform[i] = 0;
	    displ_gaussian[i] = 0;
	    baseline[i] = 0;
	}
	ss = 4.2;
	*/
    }
    virtual void TearDown() {
	delete potential;
    }
};

TEST_F(TestMC, BasicFunctionalityAddingModules){
    //mcpele::RandomCoordsDisplacement sampler_uniform(42);
    //mcpele::GaussianCoordsDisplacement sampler_gaussian(42);
    mcpele::MC mc(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc.get_potential_ptr())->get_k() );
    mc.run(max_iter);
}
