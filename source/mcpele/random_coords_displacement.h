#ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_H__
#define _MCPELE_RANDOM_COORDS_DISPLACEMENT_H__

#include <cstddef>
#include <random>

#include "mc.h"

namespace mcpele {

/**
 * Random coords displacement, generates a random displacement for a N
 * dimensional system sampling from a N-dimensional sphere.
 * The stepsize is defined per coordinates, that's why the maximum stepsize is
 * sqrt(N) * stepsize.
 */
class RandomCoordsDisplacement : public TakeStep {
 protected:
  size_t m_seed;
  std::mt19937_64 m_generator;
  std::uniform_real_distribution<double> m_real_distribution;
  double m_stepsize;
  size_t m_count;
  // if max stepsize is 0, then no max stepsize is set
  double m_max_stepsize;

 public:
  RandomCoordsDisplacement(const size_t rseed, const double stepsize = 1,
                           const double max_stepsize = 0);
  virtual ~RandomCoordsDisplacement() {}
  virtual void displace(pele::Array<double> &coords, MC *mc) = 0;
  size_t get_seed() const { return m_seed; }
  void set_generator_seed(const size_t inp) { m_generator.seed(inp); }
  double expected_mean() const { return 0; }
  double get_stepsize() const { return m_stepsize; }
  void set_stepsize(const double input) { m_stepsize = input; }
  /**
   * Reference: http://mathworld.wolfram.com/UniformDistribution.html
   */
  double expected_variance(const double ss) const {
    return ss * ss / static_cast<double>(12);
  }
  void increase_acceptance(const double factor) { m_stepsize *= factor; }
  void decrease_acceptance(const double factor) {
    if (m_max_stepsize == 0 or m_stepsize < m_max_stepsize) {
      // 0 means the max stepsize has not been set
      m_stepsize /= factor;
    }
  }
  size_t get_count() const { return m_count; }
  void set_count(const size_t input) { m_count = input; }
};

class

    RandomCoordsDisplacementAll : public RandomCoordsDisplacement {
  std::vector<size_t> m_changed_atoms;
  std::vector<double> m_changed_coords_old;

 public:
  RandomCoordsDisplacementAll(const size_t rseed, const double stepsize = 1,
                              double max_stepsize = 0);
  // The overloaded constructor helps return changed atoms and particles for
  RandomCoordsDisplacementAll(const size_t rseed, const size_t nparticles,
                              const size_t ndim, const double stepsize = 1,
                              double max_stepsize = 0);
  virtual ~RandomCoordsDisplacementAll() {}
  virtual void displace(pele::Array<double> &coords, MC *mc);
  const std::vector<size_t> get_changed_atoms() const {
    return m_changed_atoms;
  }
  const std::vector<double> get_changed_coords_old() const {
    return m_changed_coords_old;
  }
};

class RandomCoordsDisplacementSingle : public RandomCoordsDisplacement {
  size_t m_nparticles, m_ndim;
  std::uniform_int_distribution<size_t> m_int_distribution;
  std::vector<size_t> m_changed_atoms = std::vector<size_t>(1);
  std::vector<double> m_changed_coords_old;
  double max_stepsize;

 public:
  RandomCoordsDisplacementSingle(const size_t rseed, const size_t nparticles,
                                 const size_t ndim, const double stepsize = 1,
                                 const double max_stepsize = 0.0);
  virtual ~RandomCoordsDisplacementSingle() {}
  virtual void displace(pele::Array<double> &coords, MC *mc);
  size_t get_rand_particle() {
    return m_changed_atoms[0];
  }  // dangerous function, should be used only for testing purposes
  const std::vector<size_t> get_changed_atoms() const {
    return m_changed_atoms;
  }
  const std::vector<double> get_changed_coords_old() const {
    return m_changed_coords_old;
  }
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_H__
