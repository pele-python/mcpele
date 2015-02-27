#include "mcpele/gaussian_coords_displacement.h"

namespace mcpele {

GaussianTakeStep::GaussianTakeStep(const size_t rseed, const double stepsize, const size_t ndim)
    : m_seed(rseed),
      m_mean(0.0),
      m_stdev(1.0),
      m_generator(rseed),
      m_distribution(m_mean, m_stdev),
      m_stepsize(stepsize),
      m_count(0),
      m_ndim(ndim),
      m_normal_vec(ndim)
{
#ifdef DEBUG
    std::cout<<"seed TakeStep:"<<_seed<< "\n";
#endif
}

/*draw ndim random variates from N(0,1) and fill up the m_normal_vec array with them*/
inline void GaussianTakeStep::m_sample_normal_vec(){
    for(size_t i = 0; i < m_ndim; ++i){
        double randz = m_distribution(m_generator); //this is sample from N(0,1)
        m_normal_vec[i] = randz;                    //This is sampled from N(0,stepsize)
    }
}

GaussianCoordsDisplacement::GaussianCoordsDisplacement(const size_t rseed, const double stepsize, const size_t ndim)
    : GaussianTakeStep(rseed, stepsize, ndim){}

/*see https://en.wikipedia.org/wiki/Multivariate_normal_distribution*/
void GaussianCoordsDisplacement::displace(pele::Array<double>& coords, MC* mc)
{
    this->m_sample_normal_vec();
    m_normal_vec /= norm(m_normal_vec);          //sample from surface of unit hypersphere
    for(size_t i = 0; i < m_ndim; ++i){
        coords[i] += m_normal_vec[i] * m_stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
    ++m_count;
}

SampleGaussian::SampleGaussian(const size_t rseed, const double stepsize, const pele::Array<double> origin)
    : GaussianTakeStep(rseed, stepsize, origin.size()),
      m_origin(origin.copy())
    {}

/*see https://en.wikipedia.org/wiki/Multivariate_normal_distribution
 * SampleGaussian::displace ignore the coords passed to it and it sample
 * with mean centered at the origin and stdev defined by the stepsize
 * */
void SampleGaussian::displace(pele::Array<double>& coords, MC* mc)
{
    this->m_sample_normal_vec();
    for(size_t i = 0; i < m_ndim; ++i){
         coords[i] = m_origin[i] + m_normal_vec[i] * m_stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
    ++m_count;
}

} // namespace mcpele
