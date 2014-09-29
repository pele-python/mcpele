#ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_H__
#define _MCPELE_RANDOM_COORDS_DISPLACEMENT_H__

#include <random>

#include "mc.h"

namespace mcpele {

/**
 * Random coords displacement, generates a random displacement for a N
 * dimensional system sampling from a N-dimensional sphere.
 * The stepsize is defined per coordinates, that's why the maximum stepsize is
 * sqrt(N) * stepsize.
 */
class RandomCoordsDisplacement : public TakeStep {
protected:
    size_t m_seed;
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_distribution;
    double m_stepsize;
public:
    RandomCoordsDisplacement(const size_t rseed, const double stepsize);
    virtual ~RandomCoordsDisplacement() {}
    virtual void displace(pele::Array<double>& coords, MC* mc);
    size_t get_seed() const {return m_seed;}
    void set_generator_seed(const size_t inp) { m_generator.seed(inp); }
    double expected_mean() const { return 0; }
    /**
     * Reference: http://mathworld.wolfram.com/UniformDistribution.html
     */
    double expected_variance(const double ss) const { return ss * ss / static_cast<double>(12); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_H__
