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
    std::shared_ptr<pele::Harmonic> potential;
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
        potential = std::make_shared<pele::Harmonic>(origin, k, boxdim);
        max_iter = 1e5;
    }
};

TEST_F(TestMC, BasicFunctionalityAddingModulesStatic){
    mcpele::MC mc(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc.get_potential_ptr().get())->get_k() );
    shared_ptr<mcpele::RandomCoordsDisplacement> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacement>(42);
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
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc->get_potential_ptr().get())->get_k() );
    shared_ptr<mcpele::RandomCoordsDisplacement> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacement>(42);
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
    max_iter *= 10;
    mcpele::MC* mc = new mcpele::MC(potential, x, 1, stepsize);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == static_cast<pele::Harmonic*>(mc->get_potential_ptr().get())->get_k() );
    //std::cout << "initial step size: " << mc->get_stepsize() << std::endl;
    shared_ptr<mcpele::TakeStep> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacement>(42);
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
    mc->set_print_progress(true);
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
    bool return_val;
    size_t call_count;
    TrivialConfTest(bool return_val1)
        : return_val(return_val1), call_count(0)
    {}

    virtual bool test(Array<double> &trial_coords, mcpele::MC * mc)
    {
        call_count++;
        return return_val;
    }

};

struct TrivialAcceptTest : public mcpele::AcceptTest{
    bool return_val;
    size_t call_count;
    TrivialAcceptTest(bool return_val1)
        : return_val(return_val1), call_count(0)
    {}

    virtual bool test(Array<double> &trial_coords, double trial_energy,
            Array<double> & old_coords, double old_energy, double temperature,
            MC * mc)
    {
        call_count++;
        return return_val;
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

struct TrivialAction : public mcpele::Action{
    size_t call_count;
    TrivialAction()
        : call_count(0)
    {}
    virtual void action(Array<double> &coords, double energy, bool accepted,
            MC* mc)
    {
        call_count++;
    }
};

struct TrivialPotential : public pele::BasePotential{
    size_t call_count;
    virtual double get_energy(Array<double> coords)
    {
        call_count++;
        return 0.;
    }
};

class TestMCMock: public ::testing::Test{
public:
    MC * mc;
    Array<double> x0;
    std::shared_ptr<TrivialPotential> pot;
    TrivialConfTest * ct;
    TrivialConfTest * lct;
    TrivialTakestep * ts;
    TrivialAcceptTest * at;
    TrivialAction * a;

    virtual void SetUp()
    {
        x0 = Array<double>(2, 0.);
        pot = std::make_shared<TrivialPotential>();
        ct = new TrivialConfTest(true);
        lct = new TrivialConfTest(true);
        at = new TrivialAcceptTest(true);
        ts = new TrivialTakestep();
        a = new TrivialAction();
        mc = new mcpele::MC(pot, x0, 1, 1);

        mc->add_conf_test(std::shared_ptr<ConfTest>(ct));
        mc->set_takestep(std::shared_ptr<mcpele::TakeStep>(ts));
        mc->add_accept_test(std::shared_ptr<mcpele::AcceptTest>(at));
        mc->add_late_conf_test(std::shared_ptr<ConfTest>(lct));
        mc->add_action(std::shared_ptr<mcpele::Action>(a));
    }
    virtual void TearDown()
    {
        delete mc;
    }
};

TEST_F(TestMCMock, AllPass_AllCalled){
    mc->run(10);
    EXPECT_EQ(ct->call_count, size_t(10));
    EXPECT_EQ(ts->call_count, size_t(10));
    EXPECT_EQ(at->call_count, size_t(10));
    EXPECT_EQ(lct->call_count, size_t(10));
    EXPECT_EQ(a->call_count, size_t(10));
    EXPECT_EQ(mc->get_neval(), size_t(11));
    EXPECT_EQ(mc->get_iterations_count(), size_t(10));
    EXPECT_EQ(mc->get_naccept(), size_t(10));
    EXPECT_EQ(mc->get_nreject(), size_t(0));
}

TEST_F(TestMCMock, ConfTest_Fails){
    ct->return_val = false;
    mc->run(10);
    EXPECT_EQ(ct->call_count, size_t(10));
    EXPECT_EQ(ts->call_count, size_t(10));
    EXPECT_EQ(at->call_count, size_t(0));
    EXPECT_EQ(lct->call_count, size_t(0));
    EXPECT_EQ(a->call_count, size_t(10));
    EXPECT_EQ(mc->get_neval(), size_t(1));
}

TEST_F(TestMCMock, AcceptTest_Fails){
    at->return_val = false;
    mc->run(10);
    EXPECT_EQ(ct->call_count, size_t(10));
    EXPECT_EQ(ts->call_count, size_t(10));
    EXPECT_EQ(at->call_count, size_t(10));
    EXPECT_EQ(lct->call_count, size_t(0));
    EXPECT_EQ(a->call_count, size_t(10));
    EXPECT_EQ(mc->get_neval(), size_t(11));
    EXPECT_EQ(mc->get_iterations_count(), size_t(10));
    EXPECT_EQ(mc->get_naccept(), size_t(0));
    EXPECT_EQ(mc->get_nreject(), size_t(10));
}

TEST_F(TestMCMock, LateConfTest_Fails){
    lct->return_val = false;
    mc->run(10);
    EXPECT_EQ(ct->call_count, size_t(10));
    EXPECT_EQ(ts->call_count, size_t(10));
    EXPECT_EQ(at->call_count, size_t(10));
    EXPECT_EQ(lct->call_count, size_t(10));
    EXPECT_EQ(a->call_count, size_t(10));
    EXPECT_EQ(mc->get_neval(), size_t(11));
    EXPECT_EQ(mc->get_iterations_count(), size_t(10));
    EXPECT_EQ(mc->get_naccept(), size_t(0));
    EXPECT_EQ(mc->get_nreject(), size_t(10));
}
