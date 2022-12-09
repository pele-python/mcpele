/*
 * Tests for Monte Carlo Swap
*/
#include <gtest/gtest.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>




#include <algorithm>
#include <ctime>
#include <iostream>
#include <memory>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <vector>

#include "pele/inversepower.hpp"
#include "pele/meta_pow.hpp"
#include "pele/utils.hpp"

#include "mcpele/adaptive_takestep.h"
#include "mcpele/histogram.h"
#include "mcpele/metropolis_test.h"
#include "mcpele/particle_pair_swap.h"
#include "mcpele/random_coords_displacement.h"
#include "mcpele/take_step_pattern.h"
#include "mcpele/uniform_rectangular_sampling.h"
#include "mcpele/uniform_spherical_sampling.h"



// import make_shared library

using namespace pele;
using mcpele::MC;
using mcpele::TakeStep;


class TestParticlePairSwap : public::testing::Test {
    // Set up a physical system for swap tests
    std::shared_ptr<BasePotential> potential;
    MC mc;
    void SetUp() {
    constexpr size_t dim = 2;
    constexpr int pow = 2;

    size_t n_particles = 32;
    double eps = 1.0;
    double non_additivity = 0.0;
    bool exact_sum = false;

    double cutoff_factor = 1.25;

    double dmin_by_dmax = 0.449;
    double d_mean = 1.0;
    BerthierDistribution3d berthier_dist =
        BerthierDistribution3d(dmin_by_dmax, d_mean);
    Array<double> radii = berthier_dist.sample(n_particles);

    double phi = 1.2;
    double box_length = get_box_length(radii, dim, phi);
    Array<double> boxvec = {box_length, box_length};

    Array<double> coordinates =
        generate_random_coordinates(box_length, n_particles, dim);
    potential =
        std::make_shared<InverseIntPowerPeriodic<dim, pow>>(
            eps, radii, boxvec, exact_sum, non_additivity);

    double temperature = 1.0;
    mc = MC(potential, coordinates, temperature);
    int seed = 0;
    mc.add_accept_test(std::make_shared<mcpele::MetropolisTest>(seed));
    }
};







