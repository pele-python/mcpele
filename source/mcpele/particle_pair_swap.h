#ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
#define _MCPELE_PARTICLE_PAIR_SWAP_H__

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
  DiscreteUniformDistribution(const size_t rseed) : m_generator(rseed), m_real_distribution(0, 1) {};

  /*
  *  Gets a random integer between lower and upper, inclusive of lower. does not include upper
  *  ex: get_random_int(0, 3) will return 0, 1, 2
  */
  int sample(const int lower, const int upper) {
    return static_cast<int>(std::floor(m_real_distribution(m_generator) * (upper - lower) + lower));
  }


  /*
  * Samples but ignores a specific value. 
  */
  int sample_ignoring_value(const int lower, const int upper, const int ignore) {
    int uniform_sample = sample(lower, upper - 1);
    if (uniform_sample == ignore) {
      uniform_sample = upper - 1;
    }
    return uniform_sample;
  }

};

class ParticlePairSwap : public TakeStep {
private:
  size_t m_seed;
  std::mt19937_64 m_generator;
  std::uniform_int_distribution<size_t> m_distribution;
  const size_t m_nr_particles;
  const size_t m_ndim;
  std::vector<long> m_changed_atoms = std::vector<long>(2);

  std::vector<double> m_changed_coords_old;
  pele::Array<double> m_radii;
  bool m_radii_set = false;
  double max_radii_diff = 0;
  double max_diff_to_swap_radii =
      0; // maximum difference in radii allowed between swaps
  std::unordered_map<int, std::vector<int>> m_particle_to_allowed_swaps;

public:
  virtual ~ParticlePairSwap() {}
  ParticlePairSwap(const size_t seed, const size_t nr_particles,
                   const size_t ndim);
  void displace(pele::Array<double> &coords, MC *mc);
  void swap_coordinates(const size_t particle_a, const size_t particle_b,
                        pele::Array<double> &coords);
  size_t get_seed() const { return m_seed; }
  void set_generator_seed(const size_t inp);
  const std::vector<long> get_changed_atoms() const { return m_changed_atoms; }
  const std::vector<double> get_changed_coords_old() const {
    return m_changed_coords_old;
  }
  void increase_acceptance(const double factor) {
    if (max_diff_to_swap_radii < max_radii_diff) {
      max_diff_to_swap_radii *= factor;
      m_particle_to_allowed_swaps.clear();
    } else {
      std::cout << "max_diff_to_swap_radii is already at maximum value of "
                << max_radii_diff << std::endl;
      std::cout << "not increasing it further" << std::endl;
    }
  }
  void decrease_acceptance(const double factor) {
    max_diff_to_swap_radii /= factor;
    m_particle_to_allowed_swaps.clear();
  }

  // Finds allowed radii and stores them to m_particle_to_allowed_swaps if
  // necessary
  void find_allowed_radii_to_swap(const pele::Array<double> &radii,
                                  const double max_diff);
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
