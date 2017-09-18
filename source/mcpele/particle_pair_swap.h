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
    const size_t m_ndim;
    std::vector<long> m_changed_atoms = std::vector<long>(2);
    std::vector<double> m_changed_coords_old;
public:
    virtual ~ParticlePairSwap() {}
    ParticlePairSwap(const size_t seed, const size_t nr_particles, const size_t ndim);
    void displace(pele::Array<double>& coords, MC* mc);
    void swap_coordinates(const size_t particle_a, const size_t particle_b, pele::Array<double>& coords);
    size_t get_seed() const { return m_seed; }
    void set_generator_seed(const size_t inp);
    const std::vector<long> get_changed_atoms() const { return m_changed_atoms; }
    const std::vector<double> get_changed_coords_old() const { return m_changed_coords_old; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PARTICLE_PAIR_SWAP_H__
