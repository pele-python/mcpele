#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>
#include <string>
#include <gtest/gtest.h>

#include "pele/harmonic.h"

#include "mcpele/random_coords_displacement.h"
#include "mcpele/metropolis_test.h"
#include "mcpele/energy_window_test.h"
#include "mcpele/record_energy_histogram.h"
#include "mcpele/adaptive_takestep.h"
#include "mcpele/take_step_pattern.h"
#include "mcpele/take_step_probabilities.h"
#include "mcpele/progress.h"

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

using std::shared_ptr;
using mcpele::ConfTest;
using mcpele::MC;
using pele::Array;

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
    mcpele::MC mc(potential, x, 1);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc.get_potential_ptr().get())->get_k() );
    shared_ptr<mcpele::RandomCoordsDisplacementAll> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacementAll>(42, stepsize);
    EXPECT_TRUE( !mc.take_step_specified() );
    mc.set_takestep(sampler_uniform);
    EXPECT_TRUE( mc.take_step_specified() );
    shared_ptr<mcpele::MetropolisTest> metropolis = std::make_shared<mcpele::MetropolisTest>(42);
    mc.add_accept_test(metropolis);
    //mc.set_print_progress(true);
    mc.run(max_iter);
}

TEST_F(TestMC, BasicFunctionalityAddingModulesDynamic){
    mcpele::MC* mc = new mcpele::MC(potential, x, 1);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == reinterpret_cast<pele::Harmonic*>(mc->get_potential_ptr().get())->get_k() );
    shared_ptr<mcpele::RandomCoordsDisplacementAll> sampler_uniform = std::make_shared<mcpele::RandomCoordsDisplacementAll>(42, stepsize);
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
    //max_iter *= 10;
    mcpele::MC* mc = new mcpele::MC(potential, x, 1);
    EXPECT_TRUE( k == potential->get_k() );
    EXPECT_TRUE( k == static_cast<pele::Harmonic*>(mc->get_potential_ptr().get())->get_k() );
    //std::cout << "initial step size: " << mc->get_stepsize() << std::endl;
    std::shared_ptr<mcpele::TakeStep> sampler_uniform_adaptive = std::make_shared<mcpele::AdaptiveTakeStep>(std::shared_ptr<mcpele::TakeStep>(new mcpele::RandomCoordsDisplacementAll(42)), 50);
    EXPECT_TRUE( !mc->take_step_specified() );
    mc->set_takestep(sampler_uniform_adaptive);
    EXPECT_TRUE( mc->take_step_specified() );
    shared_ptr<mcpele::AcceptTest> metropolis = std::make_shared<mcpele::MetropolisTest>(42);
    mc->add_accept_test(metropolis);
    const size_t adj_iter(max_iter/1e1);
    mc->set_report_steps(adj_iter);
    shared_ptr<mcpele::Action> record_histogram = std::make_shared<mcpele::RecordEnergyHistogram>(0,10,1,adj_iter);
    mc->add_action(record_histogram);
    mc->run(max_iter);
    EXPECT_TRUE( mc->get_iterations_count() == max_iter );
    const double computed_mean = std::static_pointer_cast<mcpele::RecordEnergyHistogram>(record_histogram)->get_mean();
    const double computed_std = sqrt(std::static_pointer_cast<mcpele::RecordEnergyHistogram>(record_histogram)->get_variance());
    const double expected_mean = 0.5*ndof;
    EXPECT_NEAR(computed_mean, expected_mean, computed_std);
    delete mc;
}

struct TrivialConfTest : public mcpele::ConfTest{
    bool return_val;
    size_t call_count;
    TrivialConfTest(bool return_val1)
        : return_val(return_val1), call_count(0)
    {}

    virtual bool conf_test(Array<double> &trial_coords, mcpele::MC * mc)
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
    virtual void displace(Array<double> &coords, MC * mc=NULL)
    {
        call_count++;
    }
};

struct TrivialTakestepMessage : public mcpele::TakeStep{
    size_t call_count;
    std::string message;
    TrivialTakestepMessage(const std::string message_)
        : call_count(0),
          message(message_)
    {}
    virtual void displace(Array<double> &coords, MC * mc=NULL)
    {
        call_count++;
        std::cout << message << "\n";
        std::cout << "call_count: " << call_count << "\n";
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
    virtual ~TrivialPotential() {}
    TrivialPotential()
        : call_count(0)
    {}
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
        mc = new mcpele::MC(pot, x0, 1);

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

TEST_F(TestMCMock, Progress_Works){
    mcpele::progress stat(100);
    for (size_t i = 0; i < 100; ++i) {
        stat.next(i + 1);
        EXPECT_EQ(i + 1, stat.get_current_percentage());
    }
}

TEST_F(TestMCMock, PatternStep_Works){
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, x0, 1);
    auto step_pattern = std::make_shared<mcpele::TakeStepPattern>();
    //auto ts0 = std::make_shared<TrivialTakestepMessage>("hello 00");
    //auto ts1 = std::make_shared<TrivialTakestepMessage>("hello 01");
    //auto ts2 = std::make_shared<TrivialTakestepMessage>("hello 02");
    auto ts0 = std::make_shared<TrivialTakestep>();
    auto ts1 = std::make_shared<TrivialTakestep>();
    auto ts2 = std::make_shared<TrivialTakestep>();
    const size_t repetitions0 = 1;
    const size_t repetitions1 = 2;
    const size_t repetitions2 = 42;
    step_pattern->add_step(ts0);
    step_pattern->add_step(ts1, repetitions1);
    step_pattern->add_step(ts2, repetitions2);
    mc->set_takestep(step_pattern);
    mc->run(100);
    std::vector<size_t> repetitions;
    repetitions.push_back(repetitions0);
    repetitions.push_back(repetitions1);
    repetitions.push_back(repetitions2);
    std::vector<size_t> expected_pattern(std::accumulate(repetitions.begin(), repetitions.end(), 0));
    std::fill(expected_pattern.begin(), expected_pattern.begin() + repetitions0, 0);
    std::fill(expected_pattern.begin() + repetitions0, expected_pattern.begin() + repetitions0 + repetitions1, 1);
    std::fill(expected_pattern.begin() + repetitions0 + repetitions1, expected_pattern.end(), 2);
    const auto actual_pattern = step_pattern->get_pattern();
    EXPECT_EQ(actual_pattern.size(), expected_pattern.size());
    EXPECT_EQ(actual_pattern.size(), std::accumulate(repetitions.begin(), repetitions.end(), 0u));
    for (size_t i = 0; i < actual_pattern.size(); ++i) {
        EXPECT_EQ(actual_pattern.at(i), expected_pattern.at(i));
    }
}

TEST_F(TestMCMock, ProbabilityTakeStep_BasicWorks){
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, x0, 1);
    auto step = std::make_shared<mcpele::TakeStepProbabilities>(42);
    auto ts0 = std::make_shared<TrivialTakestep>();
    auto ts1 = std::make_shared<TrivialTakestep>();
    auto ts2 = std::make_shared<TrivialTakestep>();
    const size_t weight0 = 1;
    const size_t weight1 = 2;
    const size_t weight2 = 42;
    std::vector<double> weights;
    weights.push_back(weight0);
    weights.push_back(weight1);
    weights.push_back(weight2);
    step->add_step(ts0, weight0);
    step->add_step(ts1, weight1);
    step->add_step(ts2, weight2);
    mc->set_takestep(step);
    mc->run(1e2);
    const auto sw = step->get_weights();
    for (size_t i = 0; i < sw.size(); ++i) {
        EXPECT_DOUBLE_EQ(sw.at(i), weights.at(i));
    }
}

TEST_F(TestMCMock, AdaptiveTakeStep_WorksDown){
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, x0, 1);
    mc->add_accept_test(std::make_shared<TrivialAcceptTest>(false));
    mcpele::RandomCoordsDisplacementAll* rs = new  mcpele::RandomCoordsDisplacementAll(42, 1);
    auto astep = std::make_shared<mcpele::AdaptiveTakeStep>(std::shared_ptr<mcpele::RandomCoordsDisplacementAll>(rs), 100, 0.8);
    mc->set_takestep(astep);
    const size_t nr_iterations = 1e3;
    mc->set_report_steps(nr_iterations);
    mc->run(nr_iterations);
    EXPECT_DOUBLE_EQ(rs->get_stepsize(), pow(0.8, 10));
}

TEST_F(TestMCMock, AdaptiveTakeStep_WorksUp){
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, x0, 1);
    mc->add_accept_test(std::make_shared<TrivialAcceptTest>(true));
    mcpele::RandomCoordsDisplacementAll* rs = new  mcpele::RandomCoordsDisplacementAll(42, 1);
    auto astep = std::make_shared<mcpele::AdaptiveTakeStep>(std::shared_ptr<mcpele::RandomCoordsDisplacementAll>(rs), 100, 0.8);
    mc->set_takestep(astep);
    const size_t nr_iterations = 1e3;
    mc->set_report_steps(nr_iterations);
    mc->run(nr_iterations);
    EXPECT_DOUBLE_EQ(rs->get_stepsize(), 1/pow(0.8, 10));
}

TEST_F(TestMCMock, RecordEnergyHistogram_Works) {
    const size_t report_steps = 1e2;
    mcpele::RecordEnergyHistogram* hist_ = new mcpele::RecordEnergyHistogram(0, 42, 2, report_steps);
    auto hist = std::shared_ptr<mcpele::Action>(hist_);
    mc->add_action(hist);
    mc->set_report_steps(report_steps);
    mc->run(10 * report_steps);
    pele::Array<double> hd = hist_->get_histogram();
    EXPECT_EQ(hd.size(), 22);
    EXPECT_GE(hd[0], 0);
    for (size_t i = 1; i < hd.size(); ++i) {
        EXPECT_DOUBLE_EQ(hd[i], 0);
    }
}

