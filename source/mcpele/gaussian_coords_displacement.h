#ifndef _MCPELE_GAUSSIAN_COORDS_DISPLACEMENT_H__
#define _MCPELE_GAUSSIAN_COORDS_DISPLACEMENT_H__

#include <random>

#include "mc.h"

namespace mcpele {

/**
 * Uniform Gaussian step
 * this step samples first from the standard normal N(0, 1) and outputs a
 * random variate sampled from N(0, stepsize)
 */
class GaussianCoordsDisplacement : public TakeStep {
protected:
    size_t m_seed;
    double m_mean;
    double m_stdev;
    std::mt19937_64 m_generator;
    std::normal_distribution<double> m_distribution;
    double m_stepsize;
public:
    GaussianCoordsDisplacement(const size_t rseed, const double stepsize);
    virtual ~GaussianCoordsDisplacement() {}
    virtual void displace(pele::Array<double>& coords, MC* mc);
    size_t get_seed() const { return m_seed; }
    void set_generator_seed(const size_t inp) { m_generator.seed(inp); }
    double expected_mean() const { return 0; }
    /**
     * Reference: http://mathworld.wolfram.com/NormalDistribution.html
     */
    double expected_variance(const double ss) const { return ss * ss; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_GAUSSIAN_COORDS_DISPLACEMENT_H__
