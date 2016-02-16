#ifndef _MCPELE_UNIFORM_CUBIC_SAMPLING_H__
#define _MCPELE_UNIFORM_CUBIC_SAMPLING_H__

namespace mcpele {
    
class UniformCubicSampling : public TakeStep {
protected:
    std::mt19937_64 m_gen;
    std::uniform_real_distribution<double> m_dist05;
    pele::Array<double> m_boxvec;
    bool m_cubic;
public:
    virtual ~UniformCubicSampling() {}
    UniformCubicSampling(const size_t seed=42, const pele::Array<double> boxvec={2})
        : m_gen(seed),
          m_dist05(-0.5, 0.5),
          m_boxvec(boxvec)
    {}
    void set_generator_seed(const size_t inp) { m_gen.seed(inp); }
    virtual void displace(pele::Array<double>& coords, MC* mc)
    {
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
};    
    
} // namespace mcpele

#endif // #ifndef _MCPELE_UNIFORM_CUBIC_SAMPLING_H__
