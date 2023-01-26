#ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
#define _MCPELE_PARTICLE_PAIR_SWAP_H__

#include <cassert>
#include <cstddef>
#include <random>
#include <unordered_map>
#include <vector>

#include "mc.h"
#include "pele/array.hpp"

namespace mcpele {

/*
 * Discrete uniform distribution, generates a random integer between 0 and N-1
 * where N is the number of particles in the system.
 * This is designed not to regenerate the random sample
 *
 * While uniform_int_distribution exists in the standard library,
 * it requires respecifying the distribution if the limits are changed
 */
class DiscreteUniformDistribution {
 private:
  std::mt19937_64 m_generator;
  std::uniform_real_distribution<double> m_real_distribution;

 public:
  DiscreteUniformDistribution(const size_t rseed)
      : m_generator(rseed), m_real_distribution(0, 1){};

  /*
   *  Gets a random integer between lower and upper, inclusive of lower. does
   * not include upper ex: sample(0, 3) will return 0, 1, 2. Upper HAS
   * to be greater than lower otherwise it will return garbage
   * Also lower has to be greater than 0
   */
  size_t sample(const size_t lower, const size_t upper) {
    return static_cast<size_t>(
        std::floor(m_real_distribution(m_generator) * (upper - lower) + lower));
  }

  /*
   * Samples but ignores a specific value.  Upper HAS to be greater than lower
   * otherwise it will return garbage
   * Also lower has to be greater than 0
   */
  size_t sample_ignoring_value(const size_t lower, const size_t upper,
                               const size_t ignore) {
    size_t uniform_sample = sample(lower, upper - 1);
    if (uniform_sample == ignore) {
      uniform_sample = upper - 1;
    }
    return uniform_sample;
  }

  size_t select_random(const std::vector<size_t> &vec) {
    return vec[sample(0, vec.size())];
  }

  void set_seed(const size_t rseed) { m_generator.seed(rseed); }
};

class ParticlePairSwap : public TakeStep {
 private:
  size_t m_seed;
  DiscreteUniformDistribution m_uniform_distribution;
  const size_t m_nr_particles;
  const size_t m_ndim;
  std::vector<size_t> m_changed_atoms = std::vector<size_t>(2);

  std::vector<double> m_changed_coords_old;
  pele::Array<double> m_radii;
  bool m_radii_set = false;
  double max_radii_diff = 0;
  double max_diff_to_swap_radii =
      0;  // maximum difference in radii allowed between swaps
  std::unordered_map<size_t, std::vector<size_t>> m_particle_to_allowed_swaps;

 public:
  virtual ~ParticlePairSwap() {}
  ParticlePairSwap(const size_t seed, const size_t nr_particles,
                   const size_t ndim);
  void displace(pele::Array<double> &coords, MC *mc);
  void swap_coordinates(const size_t particle_a, const size_t particle_b,
                        pele::Array<double> &coords);
  size_t get_seed() const { return m_seed; }
  void set_generator_seed(const size_t inp);
  const std::vector<size_t> get_changed_atoms() const {
    return m_changed_atoms;
  }
  const std::vector<double> get_changed_coords_old() const {
    return m_changed_coords_old;
  }

  /*
   * Decreases the maximum difference in radii allowed between swaps
   * window_decreasing_factor must be between 0 and 1
   * The larger the difference the less likely a move will be accepted
   */
  void increase_acceptance(const double window_decreasing_factor) {
    assert(m_radii_set);  // Radii needs to be set for max_radii_diff to be
                          // initialized
    max_diff_to_swap_radii *= window_decreasing_factor;
    m_particle_to_allowed_swaps.clear();
  }
  /*
   * Increases the maximum difference in radii allowed between swaps
   * window_decreasing_factor must be between 0 and 1
   * The larger the difference the less likely a move will be accepted
   * If the maximum difference is already at the minimum value, it will not
   * decrease it further
   */
  void decrease_acceptance(const double window_decreasing_factor) {
    assert(m_radii_set);  // Radii needs to be set for max_radii_diff to be
                          // initialized
    max_diff_to_swap_radii /= window_decreasing_factor;
    m_particle_to_allowed_swaps.clear();
  }
  double get_max_diff_to_swap_radii() const { return max_diff_to_swap_radii; }

  // Finds allowed radii and stores them to m_particle_to_allowed_swaps if
  // necessary
  void add_allowed_radii(std::vector<size_t> &allowed_radii, double min_radius,
                         double max_radius, pele::Array<double> const &radii,
                         const size_t particle_a);

  void add_closest_radii(std::vector<size_t> &allowed_radii,
                         pele::Array<double> const &radii,
                         const size_t particle_a);

  std::vector<size_t> find_allowed_radii_to_swap(
      const pele::Array<double> &radii, const double max_diff,
      const size_t particle_a);
  size_t find_swap_partner(size_t particle_a);
  std::unordered_map<size_t, std::vector<size_t>>
  get_particle_to_allowed_swaps() const {
    return m_particle_to_allowed_swaps;
  }
  void print_swap_map() {
    for (size_t i = 0; i < m_radii.size(); ++i) {
      const auto allowed_swaps = m_particle_to_allowed_swaps[i];
      std::cout << "Particle " << i << " can swap with: ";
      for (size_t j = 0; j < allowed_swaps.size(); ++j) {
        std::cout << allowed_swaps[j] << " ";
      }
      std::cout << std::endl;
    }
  }
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
