#ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
#define _MCPELE_PARTICLE_PAIR_SWAP_H__

#include <random>

#include "mc.h"

namespace mcpele {

class ParticlePairSwap : public TakeStep {
private:
    size_t m_seed;
    std::mt19937_64 m_generator;
    std::uniform_int_distribution<size_t> m_distribution;
    const size_t m_nr_particles;
    std::vector<size_t> changed_atoms = std::vector<size_t>(2);
public:
    virtual ~ParticlePairSwap() {}
    ParticlePairSwap(const size_t seed, const size_t nr_particles);
    void displace(pele::Array<double>& coords, MC* mc);
    void swap_coordinates(const size_t particle_a, const size_t particle_b, pele::Array<double>& coords);
    size_t get_seed() const { return m_seed; }
    void set_generator_seed(const size_t inp);
    const std::vector<size_t> get_changed_atoms() const { return changed_atoms; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
