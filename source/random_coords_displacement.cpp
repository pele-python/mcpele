#include "mcpele/random_coords_displacement.h"

namespace mcpele {

/*RandomCoordsDisplacement*/

  RandomCoordsDisplacement::RandomCoordsDisplacement(const size_t rseed, const double stepsize)
    : m_seed(rseed),
      m_generator(rseed),
      m_real_distribution(0.0, 1.0),
      m_stepsize(stepsize),
      m_count(0)
  {
#ifdef DEBUG
    std::cout<<"seed TakeStep:"<<_seed<< "\n";
#endif
  }

  /*RandomCoordsDisplacementAll*/

  RandomCoordsDisplacementAll::RandomCoordsDisplacementAll(const size_t rseed, const double stepsize)
    : RandomCoordsDisplacement(rseed, stepsize)
  {}

  RandomCoordsDisplacementAll::RandomCoordsDisplacementAll(const size_t rseed, const size_t nparticles, const size_t ndim, const double stepsize)
    : RandomCoordsDisplacement(rseed, stepsize),
      m_changed_atoms(nparticles),
      m_changed_coords_old(nparticles*ndim)
  {}
  
  

  void RandomCoordsDisplacementAll::displace(pele::Array<double>& coords, MC* mc)
  { std::iota(m_changed_atoms.begin(), m_changed_atoms.end(), 0);
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

  RandomCoordsDisplacementSingle::RandomCoordsDisplacementSingle(const size_t rseed, const size_t nparticles, const size_t ndim, const double stepsize)
    : RandomCoordsDisplacement(rseed, stepsize),
      m_nparticles(nparticles),
      m_ndim(ndim),
      m_int_distribution(0, m_nparticles-1),
      m_changed_coords_old(ndim)
    {}

  
  void RandomCoordsDisplacementSingle::displace(pele::Array<double>& coords, MC* mc)
  {
    m_changed_atoms[0] = m_int_distribution(m_generator);
    // std::cout << m_changed_atoms[0] << "m_changed_atoms \n";
    // std::cout << m_generator << " M generaaator \n";
    size_t offset = m_changed_atoms[0] * m_ndim;
    for (size_t i = 0; i < m_ndim; ++i) {
      double rand = m_real_distribution(m_generator);
      m_changed_coords_old[i] = coords[offset + i];
      coords[offset + i] += (0.5 - rand) * m_stepsize;
    }
    ++m_count;
  }

} // namespace mcpele
