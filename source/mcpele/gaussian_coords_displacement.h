#ifndef _MCPELE_GAUSSIAN_COORDS_DISPLACEMENT_H__
#define _MCPELE_GAUSSIAN_COORDS_DISPLACEMENT_H__

#include <random>

#include "mc.h"

namespace mcpele {

/*GaussianTakeStep is a base class from which all the gaussian based take step routines should derive
 * as it implements all its basic components
 * */
class GaussianTakeStep : public TakeStep {
protected:
    size_t m_seed;
    double m_mean;
    double m_stdev;
    std::mt19937_64 m_generator;
    std::normal_distribution<double> m_distribution;
    double m_stepsize;
    size_t m_count, m_ndim;
    pele::Array<double> m_normal_vec;
    /*draw ndim random variates from N(0,1) and fill up the m_normal_vec array with them*/
    inline void m_sample_normal_vec(){
        for(size_t i = 0; i < m_ndim; ++i){
            double randz = m_distribution(m_generator); //this is sample from N(0,1)
            m_normal_vec[i] = randz;
        }
    }
public:
    GaussianTakeStep(const size_t rseed, const double stepsize, const size_t ndim);
    virtual ~GaussianTakeStep() {}
    virtual void displace(pele::Array<double>& coords, MC* mc)=0;
    size_t get_seed() const { return m_seed; }
    void set_generator_seed(const size_t inp) { m_generator.seed(inp); }
    double get_stepsize() const { return m_stepsize; }
    void set_stepsize(const double input) { m_stepsize = input; }
    size_t get_count() const { return m_count; }
    /*Reference: http://mathworld.wolfram.com/NormalDistribution.html*/
    double expected_mean() const { return 0; }
    double expected_variance(const double ss) const { return ss * ss; }
};

class GaussianCoordsDisplacement : public GaussianTakeStep {
public:
    GaussianCoordsDisplacement(const size_t rseed, const double stepsize, const size_t ndim);
    virtual ~GaussianCoordsDisplacement() {}
    virtual void displace(pele::Array<double>& coords, MC* mc);
};

/**
 * Sample a simple Gaussian distribution N(coords, stepsize)
 * this step samples first from the standard normal N(0, 1) and outputs a
 * random variate sampled from N(0, stepsize)
 */

class SampleGaussian : public GaussianTakeStep {
protected:
    pele::Array<double> m_origin;
public:
    SampleGaussian(const size_t rseed, const double stepsize, const pele::Array<double> origin);
    virtual ~SampleGaussian() {}
    virtual void displace(pele::Array<double>& coords, MC* mc);
};

} // namespace mcpele

#endif // #ifndef _MCPELE_GAUSSIAN_COORDS_DISPLACEMENT_H__
