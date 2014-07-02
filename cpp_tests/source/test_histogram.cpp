#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include "pele/array.h"

#include "mcpele/takestep.h"
#include "mcpele/histogram.h"

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

class TestHistogram: public ::testing::Test{
public:
    typedef pele::Array<double> arr_t;
    size_t ndof;
    size_t nsteps;
    size_t ntot;
    arr_t displ_uniform;
    arr_t displ_gaussian;
    arr_t baseline;
    double ss;
    virtual void SetUp(){
	ndof = 1000;
	nsteps = 100;
	ntot = nsteps*ndof;
	displ_uniform.resize(ndof);
	displ_gaussian.resize(ndof);
	baseline.resize(ndof);
	for (size_t i = 0; i < ndof; ++i){
	    displ_uniform[i] = 0;
	    displ_gaussian[i] = 0;
	    baseline[i] = 0;
	}
	ss = 4.2;
    }
    virtual void TearDown() {
	;
    }
};

TEST_F(TestHistogram, TestMoments){
    mcpele::RandomCoordsDisplacement sampler_uniform(42);
    mcpele::GaussianCoordsDisplacement sampler_gaussian(42);
    mcpele::Histogram hist_uniform(-0.5*ss,0.5*ss,0.01);
    mcpele::Histogram hist_gaussian(-5,5,0.1);
    for (size_t step = 0; step < nsteps; ++step){
	std::fill(displ_uniform.data(),displ_uniform.data()+ndof,0);
	std::fill(displ_gaussian.data(),displ_gaussian.data()+ndof,0);
	sampler_uniform.takestep(displ_uniform,ss);
	sampler_gaussian.takestep(displ_gaussian,ss);
	for (size_t dof = 0; dof < ndof; ++dof){
	    hist_uniform.add_entry(displ_uniform[dof]);
	    hist_gaussian.add_entry(displ_gaussian[dof]);
	}
    }
    EXPECT_TRUE(static_cast<size_t>(hist_uniform.entries())==ntot);
    EXPECT_TRUE(static_cast<size_t>(hist_gaussian.entries())==ntot);
    EXPECT_NEAR_RELATIVE(hist_uniform.get_mean(),sampler_uniform.expected_mean(),2*ss/sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_gaussian.get_mean(),sampler_gaussian.expected_mean(),2*ss/sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_uniform.get_variance(),sampler_uniform.expected_variance(ss),2*ss/sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_gaussian.get_variance(),sampler_gaussian.expected_variance(ss),2*ss/sqrt(ntot));
}
