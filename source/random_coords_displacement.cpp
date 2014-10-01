#include "mcpele/random_coords_displacement.h"

namespace mcpele {

RandomCoordsDisplacement::RandomCoordsDisplacement(const size_t rseed, const double stepsize)
    : m_seed(rseed),
      m_generator(rseed),
      m_distribution(0.0, 1.0),
      m_stepsize(stepsize)
{
    #ifdef DEBUG
        std::cout<<"seed TakeStep:"<<_seed<< "\n";
    #endif
}

void RandomCoordsDisplacement::displace(pele::Array<double>& coords, MC* mc)
{
    double rand;
    for (size_t i = 0; i < coords.size(); ++i) {
        rand = m_distribution(m_generator);
        coords[i] += (0.5 - rand) * m_stepsize;
    }
}

} // namespace mcpele
