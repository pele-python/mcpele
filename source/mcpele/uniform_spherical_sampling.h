#ifndef _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__
#define _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__

namespace mcpele {

/**
 * Sample points uniformly at random within an N-ball.
 * See also: http://math.stackexchange.com/questions/87230/picking-random-points-in-the-volume-of-sphere-with-uniform-probability
 */    

class UniformSphericalSampling : public TakeStep {
protected:
    std::mt19937_64 m_gen;
    const double m_radius;    
    std::normal_distribution<double> m_dist_normal;
    std::uniform_real_distribution<double> m_dist_uniform;
public:
    virtual ~UniformSphericalSampling() {}
    UniformSphericalSampling(const size_t seed=42, const double radius=1)
        : m_gen(seed),
          m_radius(radius),
          m_dist_normal(0, 1),
          m_dist_uniform(0, 1)
    {}
    void set_generator_seed(const size_t inp) { m_gen.seed(inp); }
    virtual void displace(pele::Array<double>& coords, MC* mc)
    {
        for (size_t i = 0; i < coords.size(); ++i) {
            coords[i] = m_dist_normal(m_gen);
        }
        const double tmp = std::pow(pele::dot(coords, coords), -0.5) * m_radius * std::pow(m_dist_uniform(m_gen), 1 / static_cast<double>(coords.size()));
        for (size_t i = 0; i < coords.size(); ++i) {
            coords[i] *= tmp;
        }
        
    }
};
    
} // namespace mcpele

#endif // #ifndef _MCPELE_UNIFORM_SPHERICAL_SAMPLING_H__
