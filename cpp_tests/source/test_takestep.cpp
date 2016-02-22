#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "pele/harmonic.h"
#include "pele/meta_pow.h"

#include "mcpele/adaptive_takestep.h"
#include "mcpele/histogram.h"
#include "mcpele/metropolis_test.h"
#include "mcpele/particle_pair_swap.h"
#include "mcpele/random_coords_displacement.h"
#include "mcpele/take_step_pattern.h"
#include "mcpele/take_step_probabilities.h"
#include "mcpele/uniform_rectangular_sampling.h"
#include "mcpele/uniform_spherical_sampling.h"

using pele::Array;
using mcpele::MC;
using mcpele::TakeStep;

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

class TakeStepTest: public ::testing::Test {
public:
    typedef pele::Array<double> arr_t;
    size_t seed;
    size_t nparticles;
    size_t ndim;
    size_t ndof;
    arr_t coor;
    arr_t reference;
    double stepsize;
    size_t niterations;
    double f;
    mcpele::MC* mc;
    virtual void SetUp(){
        seed = 42;
        nparticles = 10;
        ndim = 3;
        ndof = nparticles * ndim;
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

TEST_F(TakeStepTest, Global_BasicFunctionalityAveragingErasing_OneIteration) {
    //one iteration gives expected variation
    mcpele::RandomCoordsDisplacementAll displ(seed, stepsize);
    displ.displace(coor, mc);
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_NEAR(reference[i], coor[i], stepsize * 0.5);
    }
}

TEST_F(TakeStepTest, Global_BasicFunctionalityAveragingErasing_NIterations){
    //n iterations give expected variation
    mcpele::RandomCoordsDisplacementAll displ(seed, stepsize);
    for (size_t i = 0; i < niterations; ++i) {
        displ.displace(coor, mc);
    }
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_NEAR(reference[i], coor[i], f * sqrt(niterations));
    }
}

TEST_F(TakeStepTest, Single_BasicFunctionalityAveragingErasing_OneIteration){
    //one iteration gives expected variation
    mcpele::RandomCoordsDisplacementSingle displ(seed, nparticles, ndim, stepsize);
    displ.displace(coor, mc);
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_NEAR(reference[i], coor[i], stepsize * 0.5);
    }
}

TEST_F(TakeStepTest, Single_BasicFunctionality_OneParticleMoves){
    //test that only 1 particle moves
    mcpele::RandomCoordsDisplacementSingle displ(seed, nparticles, ndim, stepsize);
    displ.displace(coor, mc);
    size_t part = displ.get_rand_particle();
    //check that particle i has moved
    for (size_t j = part * ndim; j < part * ndim + ndim; ++j) {
        EXPECT_NE(reference[j], coor[j]);
    }
    //check that all other particles have not moved
    for (size_t i = 0; i < nparticles; ++i) {
        if (i != part) {
            for (size_t j = 0; j < ndim; ++j) {
                EXPECT_EQ(reference[i * ndim + j], coor[i * ndim + j]);
            }
        }
    }
}

TEST_F(TakeStepTest, Single_BasicFunctionality_AllParticlesMove){
    //test that all particles move
    mcpele::RandomCoordsDisplacementSingle displ(seed, nparticles, ndim, stepsize);
    for (size_t i = 0; i < 1000; ++i) {
        displ.displace(coor, mc);
    }
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_NE(reference[i], coor[i]);
    }
}

TEST_F(TakeStepTest, Single_BasicFunctionality_AllParticlesSampledUniformly){
    //test that particles are sampled uniformly
    mcpele::RandomCoordsDisplacementSingle displ(seed, nparticles, ndim, stepsize);
    mcpele::Histogram hist_uniform_single(0, nparticles - 1, 1);
    const size_t ntot = 5000;
    for (size_t i = 0; i < ntot; ++i) {
        displ.displace(coor, mc);
        hist_uniform_single.add_entry(displ.get_rand_particle());
    }
    EXPECT_NEAR_RELATIVE(hist_uniform_single.get_mean(), (nparticles - 1) / 2, nparticles / sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_uniform_single.get_variance(), (ntot * ntot - 1) / 12, (ntot * ntot - 1) / (12 * sqrt(ntot)));
}

TEST_F(TakeStepTest, UniformCubic_CorrectMoments){
    // test that first two moments of uniform cubic are correct
    // http://mathworld.wolfram.com/UniformDistribution.html
    const size_t ndim = 20;
    const size_t nsamples = 1e4;
    const double delta = 42.42;
    std::vector<std::shared_ptr<mcpele::Histogram> > hist(100, std::make_shared<mcpele::Histogram>(0, nparticles - 1, 1));
    mcpele::UniformRectangularSampling sampler(42, {2 * delta});
    pele::Array<double> x(ndim);
    for (size_t i = 0; i < nsamples; ++i) {
        sampler.displace(x, NULL);
        for (size_t k = 0; k < ndim; ++k) {
            hist[k]->add_entry(x[k]);
        }
    }
    for (size_t i = 0; i < ndim; ++i) {
        EXPECT_NEAR_RELATIVE(hist[i]->get_mean(), 0, 2 / sqrt(nsamples));
        EXPECT_NEAR_RELATIVE(hist[i]->get_variance(), delta * delta / 3, 2 / sqrt(nsamples));
    }
}

TEST_F(TakeStepTest, UniformSpherical_CorrectMoments){
    // test uniform spherical is OK in 2d
    const size_t nsamples = 1e5;
    const double radius = 42.42;
    mcpele::Histogram hist(0, nparticles - 1, 1);
    mcpele::Histogram hist3(0, nparticles - 1, 1);
    mcpele::Histogram hist_rectangle_2d(0, nparticles -1, 1);
    mcpele::UniformSphericalSampling sampler(42, radius);
    const double lx = 42.42;
    const double ly = 3;
    mcpele::UniformRectangularSampling sampler_rectangle(42, {lx, ly});
    pele::Array<double> x2(2);
    pele::Array<double> x3(3);
    for (size_t i = 0; i < nsamples; ++i) {
        sampler.displace(x2, NULL);
        hist.add_entry(pele::dot(x2, x2));
        sampler.displace(x3, NULL);
        hist3.add_entry(pele::dot(x3, x3));
        sampler_rectangle.displace(x2, NULL);
        hist_rectangle_2d.add_entry(pele::dot(x2, x2));
    }
    /**For a random walk in a disk of radius R, the mean of the squared
     * displacement from the origin is R^/2 and the variance of the
     * squared displacement is R^4/12.
     */
    EXPECT_NEAR_RELATIVE(radius * radius / 2, hist.get_mean(), 1e-3);
    EXPECT_NEAR_RELATIVE(pele::pos_int_pow<4>(radius) / 12, hist.get_variance(), 1e-2);
    /**For a random walk in a 3d sphere of radius R, the mean of the
     * squared displacement from the origin is 3 * R^2 / 5 and the
     * variance of the squared displacement from the origin is
     * 12 * R^4 / 175.
     */
     EXPECT_NEAR_RELATIVE(3 * pele::pos_int_pow<2>(radius) / 5, hist3.get_mean(), 1e-3);
     EXPECT_NEAR_RELATIVE(12 * pele::pos_int_pow<4>(radius) / 175, hist3.get_variance(), 1e-2);
    /**Random walk in a rectangle.
     * Mean of squared displacement is 1/12 * (lx^2 + ly^2).
     * Variance of the squared displacement is 1/180 * (lx^4 + ly^4).
     */
     const double displ2_rw = 1. / 12. * (lx * lx + ly * ly);
     const double var_displ2_rw = 1. / 180. * (pele::pos_int_pow<4>(lx) + pele::pos_int_pow<4>(ly));
     EXPECT_NEAR_RELATIVE(displ2_rw, hist_rectangle_2d.get_mean(), 1e-2);
     EXPECT_NEAR_RELATIVE(var_displ2_rw, hist_rectangle_2d.get_variance(), 1e-2);
}

TEST_F(TakeStepTest, Single_BasicFunctionalityAveragingErasing_NIterations){
    //n iterations give expected variation
    mcpele::RandomCoordsDisplacementSingle displ(seed, nparticles, ndim, stepsize);
    for (size_t i = 0; i < niterations; ++i) {
        displ.displace(coor, mc);
    }
    for (size_t i = 0; i < ndof; ++i) {
        EXPECT_NEAR( reference[i], coor[i], f * sqrt(niterations) );
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

TEST_F(TakeStepTest, SwapDisplace_Works) {
    bool throws = false;
    try {
        mcpele::ParticlePairSwap(42, 0);
    }
    catch (...) {
        throws = true;
    }
    EXPECT_TRUE(throws);
    const size_t box_dimension = 3;
    const size_t nr_particles = ndof / box_dimension;
    mcpele::ParticlePairSwap swap(42, nr_particles);
    const size_t new_seed = 44;
    EXPECT_EQ(swap.get_seed(), 42u);
    swap.set_generator_seed(new_seed);
    EXPECT_EQ(swap.get_seed(), new_seed);
    auto coor1 = coor.copy();
    auto coor2 = coor.copy();
    swap.displace(coor1, NULL);
    size_t nr_different_elements = 0;
    size_t nr_identical_elements = 0;
    for (size_t i = 0; i < coor1.size(); ++i) {
       const bool id = coor1[i] == coor2[i];
       nr_different_elements += !id;
       nr_identical_elements += id;
    }
    EXPECT_EQ(nr_different_elements + nr_identical_elements, coor1.size());
    EXPECT_EQ(nr_different_elements, 2 * box_dimension);
    EXPECT_EQ(nr_identical_elements, (nr_particles - 2) * box_dimension);
}

class TrivialTakestep : public mcpele::TakeStep {
private:
    size_t call_count;
    size_t report_count;
public:
    virtual ~TrivialTakestep() {}
    TrivialTakestep()
        : call_count(0),
          report_count(0)
    {}
    virtual void displace(Array<double> &coords, MC * mc=NULL)
    {
        call_count++;
    }
    size_t get_call_count() const { return call_count; }
    size_t get_report_count() const { return report_count; }
    virtual void report(pele::Array<double>&, const double, pele::Array<double>&, const double, const bool, MC*)
    {
        ++report_count;
    }
};

struct TrivialTakestep2 : public mcpele::TakeStep {
    size_t call_counter;
    size_t report_counter;
    size_t get_call_count() const { return call_counter; }
    size_t get_report_count() const { return report_counter; } 
    virtual ~TrivialTakestep2() {}
    TrivialTakestep2()
        : call_counter(0),
          report_counter(0)
    {}
    virtual void displace(Array<double>&, MC*) 
    {
        ++call_counter;
    }
    virtual void report(pele::Array<double>&, const double, pele::Array<double>&, const double, const bool, MC*)
    {
        ++report_counter;
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
    typedef TrivialTakestep2 tt;
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, coor, 1);
    auto step = std::make_shared<mcpele::TakeStepProbabilities>(42);
    std::shared_ptr<TakeStep> ts0(new tt());
    std::shared_ptr<TakeStep> ts1(new tt());
    std::shared_ptr<TakeStep> ts2(new tt());
    const size_t weight0 = 1;
    const size_t weight1 = 2;
    const size_t weight2 = 42;
    step->add_step(ts0, weight0);
    step->add_step(ts1, weight1);
    step->add_step(ts2, weight2);
    mc->set_takestep(step);
    const size_t total_iterations = 1e4;
    mc->run(total_iterations);
    const double total_input_weight = weight0 + weight1 + weight2;
    const double freq0 = weight0 / total_input_weight;
    const double freq1 = weight1 / total_input_weight;
    const double freq2 = weight2 / total_input_weight;
    EXPECT_NEAR(freq0, static_cast<double>(static_cast<tt*>(ts0.get())->get_call_count()) / static_cast<double>(total_iterations), 2e-3);
    EXPECT_NEAR(freq1, static_cast<double>(static_cast<tt*>(ts1.get())->get_call_count()) / static_cast<double>(total_iterations), 2e-3);
    EXPECT_NEAR(freq2, static_cast<double>(static_cast<tt*>(ts2.get())->get_call_count()) / static_cast<double>(total_iterations), 2e-3);
}

TEST_F(TakeStepTest, TakeStepProbabilities_Throws){
    bool threw = false;
    try {
        auto step = std::make_shared<mcpele::TakeStepProbabilities>(42);
        step->displace(coor, NULL);
    }
    catch (...) {
        threw = true;
    }
    EXPECT_TRUE(threw);
}

TEST_F(TakeStepTest, TakeStepPattern_Correct){
    typedef TrivialTakestep2 tt;
    auto pot = std::make_shared<TrivialPotential>();
    auto mc = std::make_shared<mcpele::MC>(pot, coor, 1);
    auto step = std::make_shared<mcpele::TakeStepPattern>();
    std::shared_ptr<TakeStep> ts0(new tt());
    std::shared_ptr<TakeStep> ts1(new tt());
    std::shared_ptr<TakeStep> ts2(new tt());
    const size_t weight0 = 1;
    const size_t weight1 = 2;
    const size_t weight2 = 42;
    std::vector<size_t> weights = {weight0, weight1, weight2};
    step->add_step(ts0, weight0);
    step->add_step(ts1, weight1);
    step->add_step(ts2, weight2);
    mc->set_takestep(step);
    const size_t total_iterations = 1e4;
    const size_t report_iterations(total_iterations / 2);
    mc->set_report_steps(report_iterations);
    EXPECT_EQ(report_iterations, mc->get_report_steps());
    mc->run(total_iterations);
    const double total_input_weight = weight0 + weight1 + weight2;
    const double freq0 = weight0 / total_input_weight;
    const double freq1 = weight1 / total_input_weight;
    const double freq2 = weight2 / total_input_weight;
    const std::vector<size_t> pattern = step->get_pattern();
    const std::vector<size_t> pattern_direct = step->get_pattern_direct();
    EXPECT_EQ(pattern.size(), pattern_direct.size());
    EXPECT_EQ(0u, pattern.front());
    EXPECT_EQ(0u, pattern_direct.front());
    EXPECT_EQ(weights.size() - 1, pattern.back());
    EXPECT_EQ(weights.size() - 1, pattern_direct.back());
    std::vector<size_t>::const_iterator i = pattern.begin();
    for (size_t w = 0; w < weights.size(); ++w) {
       for (size_t b = 0; b < weights.at(w); ++b) {
           EXPECT_EQ(*i++, w);
       } 
    }
    for (size_t i = 0; i < pattern.size(); ++i) {
        EXPECT_EQ(pattern.at(i), pattern_direct.at(i));
    }
    std::cout << "ts0->get_call_count(): " << static_cast<tt*>(ts0.get())->get_call_count() << "\n";
    std::cout << "ts1->get_call_count(): " << static_cast<tt*>(ts1.get())->get_call_count() << "\n";
    std::cout << "ts2->get_call_count(): " << static_cast<tt*>(ts2.get())->get_call_count() << "\n";
    std::cout << "total_iterations: " << total_iterations << "\n";
    EXPECT_NEAR(freq0, static_cast<double>(static_cast<tt*>(ts0.get())->get_call_count()) / static_cast<double>(total_iterations), 2e-3);
    EXPECT_NEAR(freq1, static_cast<double>(static_cast<tt*>(ts1.get())->get_call_count()) / static_cast<double>(total_iterations), 2e-3);
    EXPECT_NEAR(freq2, static_cast<double>(static_cast<tt*>(ts2.get())->get_call_count()) / static_cast<double>(total_iterations), 2e-3);
    std::cout << "ts0->get_report_count(): " << static_cast<tt*>(ts0.get())->get_report_count() << "\n";
    std::cout << "ts1->get_report_count(): " << static_cast<tt*>(ts1.get())->get_report_count() << "\n";
    std::cout << "ts2->get_report_count(): " << static_cast<tt*>(ts2.get())->get_report_count() << "\n";
    std::cout << "report_iterations: " << report_iterations << "\n";
    EXPECT_NEAR(freq0, static_cast<double>(static_cast<tt*>(ts0.get())->get_report_count()) / static_cast<double>(report_iterations), 1e-2);
    EXPECT_NEAR(freq1, static_cast<double>(static_cast<tt*>(ts1.get())->get_report_count()) / static_cast<double>(report_iterations), 1e-2);
    EXPECT_NEAR(freq2, static_cast<double>(static_cast<tt*>(ts2.get())->get_report_count()) / static_cast<double>(report_iterations), 1e-2);
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

TEST_F(AdaptiveTakeStepTest, Global_UniformAdaptive_Works) {
    const double initial_stepsize = 1e-2;
    auto uni = std::make_shared<mcpele::RandomCoordsDisplacementAll>(42, initial_stepsize);
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

TEST_F(AdaptiveTakeStepTest, Single_UniformAdaptive_Works) {
    const double initial_stepsize = 1e-2;
    auto uni = std::make_shared<mcpele::RandomCoordsDisplacementSingle>(42, nr_particles, box_dimension, initial_stepsize);
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
