#include "gaussian_coords_displacement.h"

namespace mcpele {

GaussianCoordsDisplacement::GaussianCoordsDisplacement(const size_t rseed, const double stepsize)
    : m_seed(rseed),
      m_mean(0.0),
      m_stdev(1.0),
      m_generator(rseed),
      m_distribution(m_mean, m_stdev),
      m_stepsize(stepsize)
{
    #ifdef DEBUG
        std::cout<<"seed TakeStep:"<<_seed<< "\n";
    #endif
}

void GaussianCoordsDisplacement::displace(pele::Array<double>& coords, MC* mc)
{
    //assert(coords.size() == _N);
    for(size_t i = 0; i < coords.size(); ++i){
        double randz = m_distribution(m_generator); //this is sample from N(0,1)
        coords[i] += randz * m_stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
}

} // namespace mcpele
