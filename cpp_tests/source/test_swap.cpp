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
#include <numeric>
#include <omp.h>
#include <random>
#include <stdexcept>
#include <unordered_map>
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
  EXPECT_EQ(nr_different_elements + nr_identical_elements,
            coordinates_1.size());
  EXPECT_EQ(nr_different_elements, 2 * runtime_dim);
  EXPECT_EQ(nr_identical_elements, (n_particles - 2) * runtime_dim);
}

/// Tests on random number generator

class TestDiscreteUniform : public ::testing::Test {
public:
  mcpele::DiscreteUniformDistribution dist;
  TestDiscreteUniform() : dist(0) {}
};

TEST_F(TestDiscreteUniform, Sample) {
  const size_t n_samples = 100000;
  std::vector<size_t> samples(n_samples);

  size_t lower = 0;
  size_t upper = 10;
  std::vector<size_t> histogram(upper - lower, 0);
  double probability;

  for (size_t i = 0; i < n_samples; ++i) {
    EXPECT_GE(samples[i], lower);
    EXPECT_LT(samples[i], upper);
    samples[i] = dist.sample(lower, upper);
    ++histogram[samples[i]];
  }
  for (size_t i = 0; i < histogram.size(); ++i) {
    probability = static_cast<double>(histogram[i]) / n_samples;
    EXPECT_NEAR(probability, 1. / (upper - lower), 1e-2);
  }
}

TEST_F(TestDiscreteUniform, SampleIgnoringValue) {
  const size_t n_samples = 100000;
  std::vector<size_t> samples(n_samples);

  size_t lower = 0;
  size_t upper = 10;
  size_t ignore = 5;
  std::vector<size_t> histogram(upper - lower, 0);
  double probability;

  for (size_t i = 0; i < n_samples; ++i) {
    EXPECT_GE(samples[i], lower);
    EXPECT_LT(samples[i], upper);
    samples[i] = dist.sample_ignoring_value(lower, upper, ignore);
    ++histogram[samples[i]];
  }
  for (size_t i = 0; i < histogram.size(); ++i) {
    probability = static_cast<double>(histogram[i]) / n_samples;
    if (i == ignore) {
      EXPECT_EQ(histogram[i], 0u);
      EXPECT_EQ(probability, 0.0);
    } else {
      EXPECT_NEAR(probability, 1.0 / (upper - lower - 1), 1e-2);
    }
  }
}

TEST_F(TestDiscreteUniform, SelectRandom) {
  const size_t n_samples = 100000;
  std::vector<size_t> samples(n_samples);

  std::vector<size_t> vec = {5, 6, 1, 2};

  std::unordered_map<size_t, size_t> histogram;
  for (size_t i = 0; i < vec.size(); ++i) {
    histogram[vec[i]] = 0u;
  }
  double probability;

  for (size_t i = 0; i < n_samples; ++i) {
    samples[i] = dist.select_random(vec);
    ++histogram[samples[i]];
  }
  for (size_t i = 0; i < histogram.size(); ++i) {
    probability = static_cast<double>(histogram[vec[i]]) / n_samples;
    EXPECT_NEAR(probability, 1.0 / vec.size(), 1e-2);
  }
}

// Test to check whether swapping works with a specific preset window

// Test to check whether the swap dictionary is being built correctly

// End to end test to check that the swap works