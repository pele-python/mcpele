#ifndef _MCPELE_UNIFORM_CUBIC_SAMPLING_H__
#define _MCPELE_UNIFORM_CUBIC_SAMPLING_H__

namespace mcpele {
    
class UniformCubicSampling : public TakeStep {
protected:
    std::mt19937_64 m_gen;
    std::uniform_real_distribution<double> m_dist;
    std::uniform_real_distribution<double> m_dist05;
    pele::Array<double> m_boxvec;
    bool m_cubic;
public:
    virtual ~UniformCubicSampling() {}
    UniformCubicSampling(const size_t seed=42, const double delta=1)
        : m_gen(seed),
          m_dist(-delta, delta),
          m_dist05(-0.5, 0.5),
          m_cubic(true)
    {}
    void set_generator_seed(const size_t inp) { m_gen.seed(inp); }
    void set_boxvec(const pele::Array<double>& boxvec) 
    {
        m_boxvec = boxvec;
        m_cubic = false;
    }
    virtual void displace(pele::Array<double>& coords, MC* mc)
    {
        if (m_cubic) {
            for (size_t i = 0; i < coords.size(); ++i) {
                coords[i] = m_dist(m_gen);            
            }
        }
        else {
            if (coords.size() % m_boxvec.size()) {
                throw std::runtime_error("UniformCubicSampling::displace: coods size incompatible with boxvec size");
            }
            const size_t nr_particles = coords.size() / m_boxvec.size();
            const size_t dim = m_boxvec.size();
            for (size_t i = 0; i < nr_particles; ++i) {
                for (size_t k = 0; k < dim; ++k) {
                    coords[i * dim + k] = m_boxvec[k] * m_dist05(m_gen);
                }
            }
        }
    }
};    
    
} // namespace mcpele

#endif // #ifndef _MCPELE_UNIFORM_CUBIC_SAMPLING_H__
