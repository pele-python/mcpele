#ifndef _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__
#define _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__

#include <random>

namespace mcpele {

/**
 * Sample points uniformly at random within an N-ball.
 * See also: http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
 */

class UniformSphericalSampling : public TakeStep {
protected:
    std::mt19937_64 m_gen;
    const double m_radius;
    pele::Array<double> m_origin;
    std::normal_distribution<double> m_dist_normal;
    std::uniform_real_distribution<double> m_dist_uniform;
public:
    virtual ~UniformSphericalSampling() {}
    UniformSphericalSampling(const size_t seed=42, const double radius=1)
        : m_gen(seed),
          m_radius(radius),
          m_origin(),
          m_dist_normal(0, 1),
          m_dist_uniform(0, 1)
    {}
    void set_generator_seed(const size_t inp) { m_gen.seed(inp); }
    void set_origin(const pele::Array<double> origin){m_origin = origin.copy();}
    virtual void displace(pele::Array<double>& coords, MC* mc)
    {
        for (size_t i = 0; i < coords.size(); ++i) {
            coords[i] = m_dist_normal(m_gen);
        }
        /**
         * From Numerical Recipes:
         * Picking a random point on a sphere:
         * 1) generate n independent, identically distributed, normal random numbers, y_0, ..., y_{n-1}
         * 2) get point {x} on unit sphere in n dimensions by {x} = {y} / norm({y})
         * Picking a random point inside a sphere:
         * 3) generate an additional uniform random number u in [0,1]
         * 4) compute point x_i = y_i * u^{1/n} / norm({y})
         */
        // This computes 1 / norm({y}).
        double tmp = 1.0 / norm(coords);
        // This computes u^{1/n} / norm({y}) and rescales to a sphere of radius m_radius.
        tmp *= m_radius * std::pow(m_dist_uniform(m_gen), 1.0 / coords.size());
        // This computes the sampled random point in the sphere.
        if (m_origin.empty()){
            coords *= tmp;
        }
        else{
            for (size_t i = 0; i < coords.size(); ++i) {
                coords[i] = m_origin[i] + coords[i] * tmp;
            }
        }
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__
