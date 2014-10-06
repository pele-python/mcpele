#include "mcpele/random_coords_displacement.h"

namespace mcpele {

/*RandomCoordsDisplacement*/

RandomCoordsDisplacement::RandomCoordsDisplacement(const size_t rseed, const double stepsize)
    : m_seed(rseed),
      m_generator(rseed),
      m_real_distribution(0.0, 1.0),
      m_stepsize(stepsize)
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
}

/*RandomCoordsDisplacementSingle*/

RandomCoordsDisplacementSingle::RandomCoordsDisplacementSingle(const size_t rseed, const size_t nparticles, const size_t ndim, const double stepsize)
    : RandomCoordsDisplacement(rseed, stepsize),
      m_nparticles(nparticles),
      m_ndim(ndim),
      m_int_distribution(0, m_nparticles-1){}

void RandomCoordsDisplacementSingle::displace(pele::Array<double>& coords, MC* mc)
{
    size_t rand_particle = m_int_distribution(m_generator);
    for (size_t i = rand_particle; i < rand_particle+m_ndim; ++i) {
        double rand = m_real_distribution(m_generator);
        coords[i] += (0.5 - rand) * m_stepsize;
    }
}

} // namespace mcpele
