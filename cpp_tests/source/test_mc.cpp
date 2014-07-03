#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>
#include <memory>

#include "pele/harmonic.h"

#include "mcpele/takestep.h"
#include "mcpele/accept_test.h"
#include "mcpele/actions.h"

#include "mcpele/mc.h"

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

using std::shared_ptr;

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
	origin = Array<double>(ndof);
	std::fill(origin.data(),origin.data()+ndof,0);
	x = Array<double>(ndof);
	std::fill(x.data(),x.data()+ndof,0);
	stepsize = 1e-2;
	k = 400;
	potential = new pele::Harmonic(origin, k, boxdim);
	max_iter = 1e5;
    }

    virtual void TearDown() {
	delete potential;
    }
};

TEST_F(TestMC, BasicFunctionalityAddingModulesStatic){
    mcpele::MC mc(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc.get_potential_ptr())->get_k() );
    shared_ptr<mcpele::RandomCoordsDisplacement> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacement>();
    EXPECT_TRUE( !mc.take_step_specified() );
    mc.set_takestep(sampler_uniform);
    EXPECT_TRUE( mc.take_step_specified() );
    shared_ptr<mcpele::MetropolisTest> metropolis = std::make_shared<mcpele::MetropolisTest>(42);
    mc.add_accept_test(metropolis);
    //mc.set_print_progress(true);
    mc.run(max_iter);
}

TEST_F(TestMC, BasicFunctionalityAddingModulesDynamic){
    mcpele::MC* mc = new mcpele::MC(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc->get_potential_ptr())->get_k() );
    shared_ptr<mcpele::RandomCoordsDisplacement> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacement>();
    EXPECT_TRUE( !mc->take_step_specified() );
    mc->set_takestep(sampler_uniform);
    EXPECT_TRUE( mc->take_step_specified() );
    shared_ptr<mcpele::MetropolisTest> metropolis = std::make_shared<mcpele::MetropolisTest>(42);
    mc->add_accept_test(metropolis);
    //mc->set_print_progress(true);
    mc->run(max_iter);
    delete mc;
}

TEST_F(TestMC, BasicFunctionalityPolyHarmonic){
    mcpele::MC* mc = new mcpele::MC(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == static_cast<pele::Harmonic*>(mc->get_potential_ptr())->get_k() );
    //std::cout << "initial step size: " << mc->get_stepsize() << std::endl;
    shared_ptr<mcpele::TakeStep> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacement>();
    EXPECT_TRUE( !mc->take_step_specified() );
    mc->set_takestep(sampler_uniform);
    EXPECT_TRUE( mc->take_step_specified() );
    shared_ptr<mcpele::AcceptTest> metropolis = std::make_shared<mcpele::MetropolisTest>(42);
    mc->add_accept_test(metropolis);
    //AdjustStep(double target, double factor, size_t niter, size_t navg)
    const size_t adj_iter(max_iter/1e1);
    shared_ptr<mcpele::Action> adjust_step = std::make_shared<mcpele::AdjustStep>(0.2, 0.5, adj_iter, adj_iter/1e1);
    mc->add_action(adjust_step);
    //(double min, double max, double bin, size_t eqsteps)
    shared_ptr<mcpele::Action> record_histogram = std::make_shared<mcpele::RecordEnergyHistogram>(0,10,1,adj_iter);
    mc->add_action(record_histogram);
    //mc->set_print_progress(true);
    mc->run(max_iter);
    EXPECT_TRUE( mc->get_iterations_count() == max_iter );
    //std::cout << "mean energy: " << std::endl;
    //std::cout << static_cast<mcpele::RecordEnergyHistogram*>(record_histogram)->get_mean() << std::endl;
    //std::cout << sqrt(static_cast<mcpele::RecordEnergyHistogram*>(record_histogram)->get_variance()) << std::endl;
    //std::cout << "prediced mean energy: " << std::endl;
    //std::cout << 0.5*ndof << std::endl;
    const double computed_mean = std::static_pointer_cast<mcpele::RecordEnergyHistogram>(record_histogram)->get_mean();
    const double computed_std = sqrt(std::static_pointer_cast<mcpele::RecordEnergyHistogram>(record_histogram)->get_variance());
    const double expected_mean = 0.5*ndof;
    EXPECT_NEAR(computed_mean, expected_mean, computed_std);
    //std::cout << "acc frac: " << mc->get_accepted_fraction() << std::endl;
    //std::cout << "final step size: " << mc->get_stepsize() << std::endl;
    delete mc;
}
