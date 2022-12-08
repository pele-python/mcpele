#ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
#define _MCPELE_PARTICLE_PAIR_SWAP_H__

#include <random>

#include "mc.h"
#include "pele/array.hpp"

namespace mcpele {

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
    } else {
      std::cout << "max_diff_to_swap_radii is already at maximum value of "
                << max_radii_diff << std::endl;
    }
  }
  void decrease_acceptance(const double factor) {
    max_diff_to_swap_radii /= factor;
  }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
