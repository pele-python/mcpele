#include "mcpele/particle_pair_swap.h"

namespace mcpele {

ParticlePairSwap::ParticlePairSwap(const size_t seed, const size_t nr_particles)
    : m_seed(seed),
      m_generator(seed),
      m_distribution(0, nr_particles - 1),
      m_nr_particles(nr_particles)
{
    if (nr_particles == 0) {
        throw std::runtime_error("ParticlePairSwap: illegal input");
    }
}

void ParticlePairSwap::displace(pele::Array<double>& coords, MC* mc)
{
    size_t particle_a = 42;
    size_t particle_b = 42;
    while (particle_a == particle_b) {
        particle_a = m_distribution(m_generator);
        particle_b = m_distribution(m_generator);
    }
    assert(particle_a < m_nr_particles && particle_b < m_nr_particles);
    assert(particle_a != particle_b);
    swap_coordinates(particle_a, particle_b, coords);
}

void ParticlePairSwap::swap_coordinates(const size_t particle_a, const size_t particle_b, pele::Array<double>& coords)
{
    if (particle_a == particle_b) {
        return;
    }
    const size_t box_dimension = coords.size() / m_nr_particles;
    const size_t index_a = particle_a * box_dimension;
    const size_t index_b = particle_b * box_dimension;
    double*const& x = coords.data();
    double*const& xa = x + index_a;
    double*const& xb = x + index_b;
    std::swap_ranges(xa, xa + box_dimension, xb);
}

void ParticlePairSwap::set_generator_seed(const size_t inp)
{
    m_generator.seed(inp);
    m_seed = inp;
}

} // namespace mcpele
