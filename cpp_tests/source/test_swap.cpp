/*
 * Tests for Monte Carlo Swap
 */
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
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

class TestParticlePairSwap : public ::testing::Test {
public:
  // Set up a physical system for swap tests
  std::shared_ptr<BasePotential> potential;
  Array<double> coordinates;
  int runtime_dim;
  int runtime_pow;
  size_t n_particles;
  void SetUp() {
    // need to be known at compile time for the potential
    constexpr size_t dim = 2;
    constexpr int pow = 2;
    runtime_dim = dim;
    runtime_pow = pow;

    n_particles = 32;
    double eps = 1.0;
    double non_additivity = 0.0;
    bool exact_sum = false;

    double dmin_by_dmax = 0.449;
    double d_mean = 1.0;
    BerthierDistribution3d berthier_dist =
        BerthierDistribution3d(dmin_by_dmax, d_mean);
    Array<double> radii = berthier_dist.sample(n_particles);

    double phi = 1.2;
    double box_length = get_box_length(radii, dim, phi);
    Array<double> boxvec = {box_length, box_length};

    coordinates = generate_random_coordinates(box_length, n_particles, dim);
    potential = std::make_shared<InverseIntPowerPeriodic<dim, pow>>(
        eps, radii, boxvec, exact_sum, non_additivity);
  }
};

TEST_F(TestParticlePairSwap, test_swap) {

  double temperature = 1.0;
  MC mc = MC(potential, coordinates, temperature);
  int seed = 0;

  mc.add_accept_test(std::make_shared<mcpele::MetropolisTest>(seed));
  mcpele::ParticlePairSwap swap(42, n_particles, runtime_dim);

  const size_t new_seed = 44;
  EXPECT_EQ(swap.get_seed(), 42u);
  swap.set_generator_seed(new_seed);
  EXPECT_EQ(swap.get_seed(), new_seed);
  auto coordinates_1 = coordinates.copy();
  auto coordinates_2 = coordinates.copy();

  // No need to add swap to mc
  swap.displace(coordinates_1, &mc);
  size_t nr_different_elements = 0;
  size_t nr_identical_elements = 0;
  for (size_t i = 0; i < coordinates_1.size(); ++i) {
    const bool id = coordinates_1[i] == coordinates_2[i];
    nr_different_elements += !id;
    nr_identical_elements += id;
  }
  EXPECT_EQ(nr_different_elements + nr_identical_elements, coordinates_1.size());
  EXPECT_EQ(nr_different_elements, 2 * runtime_dim);
  EXPECT_EQ(nr_identical_elements, (n_particles - 2) * runtime_dim);
}
