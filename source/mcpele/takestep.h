#ifndef _MCPELE_TAKESTEP_H__
#define _MCPELE_TAKESTEP_H__

#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

#include "pele/array.h"
#include "mc.h"

namespace mcpele{

/**
 * Random coords displacement, generates a random displacement for a N dimensional system
 * sampling from a N-dimensional sphere
 * the stepsize is defined per coordinates, that's why the maximum stepsize is sqrt(N)*stepsize
 */

class RandomCoordsDisplacement : public TakeStep {
protected:
    size_t _seed;
    std::mt19937_64 _generator;
    std::uniform_real_distribution<double> _distribution;
public:
    RandomCoordsDisplacement(size_t rseed);
    virtual ~RandomCoordsDisplacement() {}
    virtual void takestep(pele::Array<double>& coords, double stepsize, MC * mc);
    size_t get_seed() const {return _seed;}
    void set_generator_seed(const size_t inp) { _generator.seed(inp); }
    double expected_mean() const { return 0; }
    /**
     * Reference: http://mathworld.wolfram.com/UniformDistribution.html
     */
    double expected_variance(const double ss) const { return ss * ss / static_cast<double>(12); }
};

class RandomCoordsDisplacementAdaptive : public RandomCoordsDisplacement {
protected:
    double m_stepsize;
    const double m_factor;
    const double m_min_acc_frac;
    const double m_max_acc_frac;
public:
    virtual ~RandomCoordsDisplacementAdaptive() {}
    RandomCoordsDisplacementAdaptive(const double stepsize=1, const double factor=0.9, const double min_acc_frac=0.2, const double max_acc_frac=0.5)
        : m_stepsize(stepsize),
          m_factor(factor),
          m_min_acc_frac(min_acc_frac),
          m_max_acc_frac(max_acc_frac)
    {}
    void displace(pele::Array<double> &coords, MC * mc)
    {
        takestep(coords, m_stepsize, mc);
    }
    void increase_acceptance()
    {
        m_stepsize *= factor;
    }
    void decrease_acceptance()
    {
        m_stepsize /= factor;
    }
    double get_min_acc_frac() const { return m_min_acc_frac; }
    double get_max_acc_frac() const { return m_max_acc_frac; }
};

/**
 * Uniform Gaussian step
 * this step samples first from the standard normal N(0,1) and outputs a random variate sampled from N(0,stepsize)
 */

class GaussianCoordsDisplacement : public TakeStep {
protected:
    size_t _seed;
    double _mean, _stdev;
    std::mt19937_64 _generator;
    std::normal_distribution<double> _distribution;
public:
    GaussianCoordsDisplacement(size_t rseed);
    virtual ~GaussianCoordsDisplacement() {}
    virtual void takestep(pele::Array<double>& coords, double stepsize, MC * mc);
    size_t get_seed() const { return _seed; }
    void set_generator_seed(const size_t inp) { _generator.seed(inp); }
    double expected_mean() const { return 0; }
    /**
     * Reference: http://mathworld.wolfram.com/NormalDistribution.html
     */
    double expected_variance(const double ss) const { return ss * ss; }
};

class ParticlePairSwap : public TakeStep {
private:
    size_t m_seed;
    std::mt19937_64 m_generator;
    std::uniform_int_distribution<size_t> m_distribution;
    const size_t m_nr_particles;
    const size_t m_swap_every;
public:
    virtual ~ParticlePairSwap() {}
    ParticlePairSwap(const size_t, const size_t, const size_t);
    virtual void takestep(pele::Array<double>& coords, double stepsize, MC * mc);
    void swap_coordinates(const size_t, const size_t, pele::Array<double>&);
    size_t get_seed() const { return m_seed; }
    void set_generator_seed(const size_t inp)
    {
        m_generator.seed(inp);
        m_seed = inp;
    }
};

}
#endif
