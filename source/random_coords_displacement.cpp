#include "mcpele/random_coords_displacement.h"
#include <ostream>

namespace mcpele {

/*RandomCoordsDisplacement*/

RandomCoordsDisplacement::RandomCoordsDisplacement(const size_t rseed,
                                                   const double stepsize,
                                                   const double max_stepsize)
    : m_seed(rseed),
      m_generator(rseed),
      m_real_distribution(0.0, 1.0),
      m_stepsize(stepsize),
      m_count(0),
      m_max_stepsize(max_stepsize) {
#ifdef DEBUG
  std::cout << "seed TakeStep:" << _seed << "\n";
#endif
}

/*RandomCoordsDisplacementAll*/

RandomCoordsDisplacementAll::RandomCoordsDisplacementAll(
    const size_t rseed, const double stepsize, const double max_stepsize)
    : RandomCoordsDisplacement(rseed, stepsize, max_stepsize) {}

RandomCoordsDisplacementAll::RandomCoordsDisplacementAll(
    const size_t rseed, const size_t nparticles, const size_t ndim,
    const double stepsize, const double max_stepsize)
    : RandomCoordsDisplacement(rseed, stepsize, max_stepsize),
      m_changed_atoms(nparticles),
      m_changed_coords_old(nparticles * ndim) {}

void RandomCoordsDisplacementAll::displace(pele::Array<double> &coords,
                                           MC *mc) {
  std::iota(m_changed_atoms.begin(), m_changed_atoms.end(), 0);
  if (m_changed_coords_old.size() != 0) {
    std::copy(coords.begin(), coords.end(), m_changed_coords_old.begin());
  }

  for (size_t i = 0; i < coords.size(); ++i) {
    double rand = m_real_distribution(m_generator);
    coords[i] += (0.5 - rand) * m_stepsize;
  }
  ++m_count;
}

/*RandomCoordsDisplacementSingle*/

RandomCoordsDisplacementSingle::RandomCoordsDisplacementSingle(
    const size_t rseed, const size_t nparticles, const size_t ndim,
    const double stepsize, const double max_stepsize)
    : RandomCoordsDisplacement(rseed, stepsize, max_stepsize),
      m_nparticles(nparticles),
      m_ndim(ndim),
      m_int_distribution(0, m_nparticles - 1),
      m_changed_coords_old(ndim) {}

void RandomCoordsDisplacementSingle::displace(pele::Array<double> &coords,
                                              MC *mc) {
  m_changed_atoms[0] = m_int_distribution(m_generator);
  size_t offset = m_changed_atoms[0] * m_ndim;
  for (size_t i = 0; i < m_ndim; ++i) {
    double rand = m_real_distribution(m_generator);
    m_changed_coords_old[i] = coords[offset + i];
    coords[offset + i] += (0.5 - rand) * m_stepsize;
  }
  ++m_count;
}

}  // namespace mcpele
