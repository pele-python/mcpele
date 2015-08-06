#ifndef _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__
#define _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__

namespace mcpele {

/**
 * Sample points uniformly at random within an N-ball.
 * See also: http://mathworld.wolfram.com/BallPointPicking.html
 */    

class UniformSphericalSampling : public TakeStep {
protected:
    std::mt19937_64 m_gen;
    const double m_radius2;    
    std::uniform_real_distribution<double> m_dist;
public:
    virtual ~UniformSphericalSampling() {}
    UniformSphericalSampling(const size_t seed=42, const double radius=1)
        : m_gen(seed),
          m_radius2(radius * radius),
          m_dist(-radius, radius)
    {}
    void set_generator_seed(const size_t inp) { m_gen.seed(inp); }
    virtual void displace(pele::Array<double>& coords, MC* mc)
    {
        do {
            for (size_t i = 0; i < coords.size(); ++i) {
                coords[i] = m_dist(m_gen);
            }
        } while (pele::dot(coords, coords) > m_radius2);
    }
};
    
} // namespace mcpele

#endif // #ifndef _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__
