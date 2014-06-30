#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include "pele/harmonic.h"

#include "takestep.h"
#include "accept_test.h"
#include "actions.h"

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
	nparticles = 1e1;
	ndof = boxdim*nparticles;
	origin.resize(ndof);
	std::fill(origin.data(),origin.data()+ndof,0);
	x.resize(ndof);
	std::fill(x.data(),x.data()+ndof,0);
	stepsize = 1e-2;
	k = 400;
	potential = new pele::Harmonic(origin, k, boxdim);
	max_iter = 1e4;
    }

    virtual void TearDown() {
	delete potential;
    }
};

TEST_F(TestMC, BasicFunctionalityAddingModulesStatic){
    mcpele::MC mc(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc.get_potential_ptr())->get_k() );
    mcpele::RandomCoordsDisplacement sampler_uniform;
    EXPECT_TRUE( !mc.take_step_specified() );
    mc.set_takestep(&sampler_uniform);
    EXPECT_TRUE( mc.take_step_specified() );
    mcpele::MetropolisTest metropolis(42);
    mc.add_accept_test(&metropolis);
    //mc.set_print_progress(true);
    mc.run(max_iter);
}

TEST_F(TestMC, BasicFunctionalityAddingModulesDynamic){
    mcpele::MC* mc = new mcpele::MC(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc->get_potential_ptr())->get_k() );
    mcpele::RandomCoordsDisplacement* sampler_uniform = new mcpele::RandomCoordsDisplacement;
    EXPECT_TRUE( !mc->take_step_specified() );
    mc->set_takestep(sampler_uniform);
    EXPECT_TRUE( mc->take_step_specified() );
    mcpele::MetropolisTest* metropolis = new mcpele::MetropolisTest(42);
    mc->add_accept_test(metropolis);
    //mc->set_print_progress(true);
    mc->run(max_iter);
    delete mc;
    delete sampler_uniform;
    delete metropolis;
}

TEST_F(TestMC, BasicFunctionalityAddingModulesPoly){
    mcpele::MC* mc = new mcpele::MC(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc->get_potential_ptr())->get_k() );
    std::cout << "initial step size: " << mc->get_stepsize() << std::endl;
    mcpele::TakeStep* sampler_uniform = new mcpele::RandomCoordsDisplacement;
    EXPECT_TRUE( !mc->take_step_specified() );
    mc->set_takestep(sampler_uniform);
    EXPECT_TRUE( mc->take_step_specified() );
    mcpele::AcceptTest* metropolis = new mcpele::MetropolisTest(42);
    mc->add_accept_test(metropolis);
    //AdjustStep(double target, double factor, size_t niter, size_t navg)
    mcpele::Action* adjust_step = new mcpele::AdjustStep(0.2, 0.5, max_iter/1e1, max_iter/1e2);
    mc->add_action(adjust_step);
    //mc->set_print_progress(true);
    mc->run(max_iter);
    EXPECT_TRUE( mc->get_iterations_count() == max_iter );
    //std::cout << "acc frac: " << mc->get_accepted_fraction() << std::endl;
    std::cout << "final step size: " << mc->get_stepsize() << std::endl;
    delete mc;
    delete sampler_uniform;
    delete metropolis;
    delete adjust_step;
}
