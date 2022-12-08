#include <algorithm>
#include <cstddef>
#include <limits>
#include <memory>
#include <vector>

#include "mcpele/particle_pair_swap.h"
#include "pele/array.hpp"
#include "pele/pairwise_potential_interface.hpp"

namespace mcpele {

ParticlePairSwap::ParticlePairSwap(const size_t seed, const size_t nr_particles,
                                   const size_t ndim)
    : m_seed(seed), m_uniform_distribution(0), m_nr_particles(nr_particles),
      m_ndim(ndim), m_changed_coords_old(2 * ndim), m_radii(nr_particles) {
  if (nr_particles == 0) {
    throw std::runtime_error("ParticlePairSwap: illegal input");
  }
}

void ParticlePairSwap::displace(pele::Array<double> &coords, MC *mc) {
  if (!m_radii_set) {
    auto pairwise_ptr = static_pointer_cast<pele::PairwisePotentialInterface>(
        mc->get_potential_ptr());
    if (pairwise_ptr == nullptr) {
      throw std::runtime_error("ParticlePairSwap: potential must be pairwise");
    }
    m_radii = pairwise_ptr->get_radii();

    auto [min, max] = std::minmax_element(m_radii.begin(), m_radii.end());
    max_radii_diff = *max - *min;
    max_diff_to_swap_radii =
        max_radii_diff + std::numeric_limits<double>::epsilon();

    m_radii_set = true;
  }

  size_t first_particle = 42;
  size_t second_particle = 42;
  first_particle = m_uniform_distribution.sample(0, m_nr_particles);

  second_particle = m_uniform_distribution.sample_ignoring_value(
      0, m_nr_particles, first_particle);
  m_changed_atoms[0] = first_particle;
  m_changed_atoms[1] = second_particle;
  swap_coordinates(first_particle, second_particle, coords);
}

void ParticlePairSwap::swap_coordinates(const size_t particle_a,
                                        const size_t particle_b,
                                        pele::Array<double> &coords) {
  const size_t index_a = particle_a * m_ndim;
  const size_t index_b = particle_b * m_ndim;
  double *const &x = coords.data();
  double *const &xa = x + index_a;
  double *const &xb = x + index_b;
  std::copy(xa, xa + m_ndim, m_changed_coords_old.begin());
  std::copy(xb, xb + m_ndim, m_changed_coords_old.begin() + m_ndim);
  if (particle_a != particle_b) {
    std::swap_ranges(xa, xa + m_ndim, xb);
  }
}

void ParticlePairSwap::set_generator_seed(const size_t inp) {
  m_uniform_distribution.set_seed(inp);
  m_seed = inp;
}

inline void ParticlePairSwap::add_allowed_radii(
    std::vector<size_t> &allowed_radii, double min_radius, double max_radius,
    pele::Array<double> const &radii, const size_t particle_a) {
  for (size_t i = 0; i < radii.size(); ++i) {
    if (i != particle_a) {
      const double radius_b = radii[i];
      if (radius_b >= min_radius && radius_b <= max_radius) {
        allowed_radii.push_back(i);
      }
    }
  }
}

inline void
ParticlePairSwap::add_closest_radii(std::vector<size_t> &allowed_radii,
                                    const pele::Array<double> &radii,
                                    const size_t particle_a) {
  double min_diff = std::numeric_limits<double>::max();
  for (size_t i = 0; i < radii.size(); ++i) {
    if (i != particle_a) {
      const double radius_b = radii[i];
      const double diff = std::abs(radius_b - radii[particle_a]);
      if (diff < min_diff) {
        min_diff = diff;
        allowed_radii.clear();
        allowed_radii.push_back(i);
      } else if (diff == min_diff) {
        allowed_radii.push_back(i);
      }
    }
  }
}

/*
 * Find all radii that are within max_diff of the radius of particle_a.
 * If no radii are found, find the closest radius. and add all indices
 * of radii that have the same distance to particle_a.
 */
std::vector<size_t>
ParticlePairSwap::find_allowed_radii_to_swap(const pele::Array<double> &radii,
                                             const double max_diff,
                                             const size_t particle_a) {

  const double radius_a = radii[particle_a];
  const double min_radius = radius_a - max_diff;
  const double max_radius = radius_a + max_diff;

  std::vector<size_t> allowed_radii_to_swap;

  add_allowed_radii(allowed_radii_to_swap, min_radius, max_radius, radii,
                    particle_a);

  if (allowed_radii_to_swap.size() == 0) {
    add_closest_radii(allowed_radii_to_swap, radii, particle_a);
  }

  return allowed_radii_to_swap;
}

size_t ParticlePairSwap::find_swap_partner(const size_t particle_a) {
  if (max_diff_to_swap_radii > max_radii_diff) {
    return m_uniform_distribution.sample_ignoring_value(0, m_radii.size(),
                                                        particle_a);
  }

  if (m_particle_to_allowed_swaps.find(particle_a) ==
      m_particle_to_allowed_swaps.end()) {
    const auto allowed_radii_to_swap =
        find_allowed_radii_to_swap(m_radii, max_diff_to_swap_radii, particle_a);
    m_particle_to_allowed_swaps[particle_a] = allowed_radii_to_swap;
  }
  const auto &allowed_radii_to_swap = m_particle_to_allowed_swaps[particle_a];

  return m_uniform_distribution.select_random(allowed_radii_to_swap);
}

} // namespace mcpele
