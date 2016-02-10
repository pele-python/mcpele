#ifndef _MCPELE_UNIFORM_CUBIC_SAMPLING_H__
#define _MCPELE_UNIFORM_CUBIC_SAMPLING_H__

namespace mcpele {
    
class UniformCubicSampling : public TakeStep {
protected:
    std::mt19937_64 m_gen;
    std::uniform_real_distribution<double> m_dist;
public:
    virtual ~UniformCubicSampling() {}
    UniformCubicSampling(const size_t seed=42, const double delta=1)
        : m_gen(seed),
          m_dist(-delta, delta)
    {}
    void set_generator_seed(const size_t inp) { m_gen.seed(inp); }
    virtual void displace(pele::Array<double>& coords, MC* mc)
    {
        for (size_t i = 0; i < coords.size(); ++i) {
            coords[i] = m_dist(m_gen);
        }
    }
};    
    
} // namespace mcpele

#endif // #ifndef _MCPELE_UNIFORM_CUBIC_SAMPLING_H__
