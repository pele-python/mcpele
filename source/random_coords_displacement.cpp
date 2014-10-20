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
    : RandomCoordsDisplacement(rseed, stepsize){}

void RandomCoordsDisplacementAll::displace(pele::Array<double>& coords, MC* mc)
{
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
      m_rand_particle(0),
      m_int_distribution(0, m_nparticles-1){}

void RandomCoordsDisplacementSingle::displace(pele::Array<double>& coords, MC* mc)
{
    m_rand_particle = m_int_distribution(m_generator);
    size_t rand_particle_dof = m_rand_particle * m_ndim;
    for (size_t i = rand_particle_dof; i < rand_particle_dof + m_ndim; ++i) {
        double rand = m_real_distribution(m_generator);
        coords[i] += (0.5 - rand) * m_stepsize;
    }
    ++m_count;
}

} // namespace mcpele
