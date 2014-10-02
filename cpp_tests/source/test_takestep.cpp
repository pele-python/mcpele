#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "pele/harmonic.h"

#include "mcpele/random_coords_displacement_adaptive.h"
#include "mcpele/particle_pair_swap.h"
#include "mcpele/metropolis_test.h"
#include "mcpele/adaptive_takestep.h"

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
    mcpele::ParticlePairSwap swap(42, nr_particles);
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

struct TrivialTakestep : public mcpele::TakeStep{
    size_t call_count;
    virtual ~TrivialTakestep() {}
    TrivialTakestep()
        : call_count(0)
    {}
    virtual void displace(Array<double> &coords, MC * mc=NULL)
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

TEST_F(TakeStepTest, TakeStepProbabilities_Correct){
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, x0, 1);
    auto step = std::make_shared<mcpele::TakeStepProbabilities>(42);
    auto ts0 = std::make_shared<TrivialTakestep>();
    auto ts1 = std::make_shared<TrivialTakestep>();
    auto ts2 = std::make_shared<TrivialTakestep>();
    const size_t weight0 = 1;
    const size_t weight1 = 2;
    const size_t weight2 = 42;
    step->add_step(ts0, weight0);
    step->add_step(ts1, weight1);
    step->add_step(ts2, weight2);
    mc->set_takestep(step);
    const size_t total_iterations = 1e4;
    mc->run(total_iterations);

}

class AdaptiveTakeStepTest: public ::testing::Test{
public:
    size_t seed;
    size_t nr_particles;
    size_t box_dimension;
    size_t ndof;
    Array<double> coords;
    Array<double> reference;
    double stepsize;
    size_t niterations;
    double k;
    std::shared_ptr<pele::Harmonic> potential;
    double temperature;
    std::shared_ptr<mcpele::MC> mc;
    std::shared_ptr<mcpele::MetropolisTest> mrt2;
    virtual void SetUp(){
        seed = 42;
        box_dimension = 3;
        nr_particles = 10;
        ndof = box_dimension * nr_particles;
        coords = Array<double>(ndof);
        for (size_t i = 0; i < ndof; ++i) {
            coords[i] = 4242 + i;
        }
        reference = coords.copy();
        stepsize = 0.1;
        niterations = 10000;
        k = 100;
        potential = std::make_shared<pele::Harmonic>(coords, k, box_dimension);
        temperature = 1;
        mc = std::make_shared<mcpele::MC>(potential, coords, temperature);
        mrt2 = std::make_shared<mcpele::MetropolisTest>(42);
        mc->add_accept_test(mrt2);
    }
};

TEST_F(AdaptiveTakeStepTest, UniformAdaptive_Works) {
    const double initial_stepsize = 1e-2;
    auto uni = std::make_shared<mcpele::RandomCoordsDisplacementAdaptive>(42, initial_stepsize);
    const size_t total_iterations(2e5);
    const size_t equilibration_report_iterations(total_iterations / 2);
    mc->set_report_steps(equilibration_report_iterations);
    const size_t interval(equilibration_report_iterations / 1e2);
    const double factor(0.9);
    const double min_acc(0.2);
    const double max_acc(0.3);
    auto uni_adaptive = std::make_shared<mcpele::AdaptiveTakeStep>(uni, interval, factor, min_acc, max_acc);
    mc->set_takestep(uni_adaptive);
    mc->run(equilibration_report_iterations);
    const size_t total_eq = total_iterations - equilibration_report_iterations;
    size_t eq_acc = 0;
    for (size_t i = 0; i < total_eq; ++i) {
        mc->one_iteration();
        eq_acc += mc->get_success();
    }
    const double eq_acc_ratio = static_cast<double>(eq_acc) / static_cast<double>(total_eq);
    EXPECT_EQ(mc->get_iterations_count(), total_iterations);
    EXPECT_LE(initial_stepsize, uni->get_stepsize());
    EXPECT_DOUBLE_EQ(mc->get_nreject(), mc->get_E_rejection_fraction() * mc->get_iterations_count());
    EXPECT_LE(eq_acc_ratio, max_acc);
    EXPECT_LE(min_acc, eq_acc_ratio);
}

