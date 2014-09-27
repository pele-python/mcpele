#include <cassert>
#include <algorithm>

#include "takestep.h"

namespace mcpele{

RandomCoordsDisplacement::RandomCoordsDisplacement(const size_t rseed, const double stepsize)
    : m_seed(rseed),
      m_generator(rseed),
      m_distribution(0.0, 1.0),
      m_stepsize(stepsize)
{
    #ifdef DEBUG
        std::cout<<"seed TakeStep:"<<_seed<< "\n";
    #endif
}

void RandomCoordsDisplacement::displace(pele::Array<double>& coords, MC* mc)
{
    double rand;
    for (size_t i = 0; i < coords.size(); ++i) {
        rand = m_distribution(m_generator);
        coords[i] += (0.5 - rand) * m_stepsize;
    }
}

GaussianCoordsDisplacement::GaussianCoordsDisplacement(const size_t rseed, const double stepsize)
    : m_seed(rseed),
      m_mean(0.0),
      m_stdev(1.0),
      m_generator(rseed),
      m_distribution(m_mean, m_stdev),
      m_stepsize(stepsize)
{
    #ifdef DEBUG
        std::cout<<"seed TakeStep:"<<_seed<< "\n";
    #endif
}

void GaussianCoordsDisplacement::displace(pele::Array<double>& coords, MC* mc)
{
    //assert(coords.size() == _N);
    for(size_t i = 0; i < coords.size(); ++i){
        double randz = m_distribution(m_generator); //this is sample from N(0,1)
        coords[i] += randz * m_stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
}

ParticlePairSwap::ParticlePairSwap(const size_t seed, const size_t nr_particles, const size_t swap_every)
    : m_seed(seed),
      m_generator(seed),
      m_distribution(0, nr_particles - 1),
      m_nr_particles(nr_particles),
      m_swap_every(swap_every)
{
    if (nr_particles * swap_every == 0) {
        throw std::runtime_error("ParticlePairSwap: illegal input");
    }
}

void ParticlePairSwap::displace(pele::Array<double>& coords, MC* mc)
{
    if (mc->get_iterations_count() % m_swap_every != 0) {
        return;
    }
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

} // namespace mcpele
