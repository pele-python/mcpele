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
using mcpele::ConfTest;
using mcpele::MC;

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

struct TrivialConfTest : public mcpele::ConfTest{
    bool result;
    size_t call_count;
    TrivialConfTest(bool return_val)
        : result(return_val), call_count(0)
    {}

    virtual bool test(Array<double> &trial_coords, mcpele::MC * mc)
    {
        call_count++;
        return result;
    }

};

struct TrivialAcceptTest : public mcpele::AcceptTest{
    bool result;
    size_t call_count;
    TrivialAcceptTest(bool return_val)
        : result(return_val), call_count(0)
    {}

    virtual bool test(Array<double> &trial_coords, double trial_energy,
            Array<double> & old_coords, double old_energy, double temperature,
            MC * mc)
    {
        call_count++;
        return result;
    }

};

struct TrivialTakestep : public mcpele::TakeStep{
    size_t call_count;
    TrivialTakestep()
        : call_count(0)
    {}
    virtual void takestep(Array<double> &coords, double stepsize, MC * mc=NULL)
    {
        call_count++;
    }
};

//class TestMCMock: public ::testing::Test{

TEST_F(TestMC, ConfTest_Fails){
    TrivialConfTest * ct = new TrivialConfTest(false);
    std::shared_ptr<ConfTest> ctest( ct );
    TrivialTakestep * ts = new TrivialTakestep();
    mcpele::MC mc(potential, x, 1, stepsize);
    mc.set_takestep(std::shared_ptr<mcpele::TakeStep>(ts));
    mc.add_conf_test(ctest);
    mc.run(10);
    EXPECT_EQ(ct->call_count, 10);
    EXPECT_EQ(ts->call_count, 10);
    EXPECT_EQ(mc.get_neval(), 1);
}

TEST_F(TestMC, ConfTest_Passes){
    mcpele::MC mc(potential, x, 1, stepsize);

    TrivialConfTest * ct = new TrivialConfTest(true);
    mc.add_conf_test(std::shared_ptr<ConfTest>(ct));

    TrivialTakestep * ts = new TrivialTakestep();
    mc.set_takestep(std::shared_ptr<mcpele::TakeStep>(ts));

    mc.run(10);
    EXPECT_EQ(ct->call_count, 10);
    EXPECT_EQ(ts->call_count, 10);
    EXPECT_EQ(mc.get_neval(), 11);
}

TEST_F(TestMC, AcceptTest_Passes){
    mcpele::MC mc(potential, x, 1, stepsize);

    TrivialConfTest * ct = new TrivialConfTest(true);
    mc.add_conf_test(std::shared_ptr<ConfTest>(ct));

    TrivialTakestep * ts = new TrivialTakestep();
    mc.set_takestep(std::shared_ptr<mcpele::TakeStep>(ts));

    TrivialAcceptTest * at = new TrivialAcceptTest(true);
    mc.add_accept_test(std::shared_ptr<mcpele::AcceptTest>(at));

    TrivialConfTest * lct = new TrivialConfTest(true);
    mc.add_late_conf_test(std::shared_ptr<ConfTest>(lct));


    mc.run(10);
    EXPECT_EQ(ct->call_count, 10);
    EXPECT_EQ(ts->call_count, 10);
    EXPECT_EQ(at->call_count, 10);
    EXPECT_EQ(lct->call_count, 10);
    EXPECT_EQ(mc.get_neval(), 11);
    EXPECT_EQ(mc.get_iterations_count(), 10);
}

TEST_F(TestMC, AcceptTest_Fails){
    mcpele::MC mc(potential, x, 1, stepsize);

    TrivialConfTest * ct = new TrivialConfTest(true);
    mc.add_conf_test(std::shared_ptr<ConfTest>(ct));

    TrivialTakestep * ts = new TrivialTakestep();
    mc.set_takestep(std::shared_ptr<mcpele::TakeStep>(ts));

    TrivialAcceptTest * at = new TrivialAcceptTest(false);
    mc.add_accept_test(std::shared_ptr<mcpele::AcceptTest>(at));

    TrivialConfTest * lct = new TrivialConfTest(true);
    mc.add_late_conf_test(std::shared_ptr<ConfTest>(lct));


    mc.run(10);
    EXPECT_EQ(ct->call_count, 10);
    EXPECT_EQ(ts->call_count, 10);
    EXPECT_EQ(at->call_count, 10);
    EXPECT_EQ(lct->call_count, 0);
    EXPECT_EQ(mc.get_neval(), 11);
    EXPECT_EQ(mc.get_iterations_count(), 10);
}
