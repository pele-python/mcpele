#ifndef _MCPELE_ACCEPT_TEST_H__
#define _MCPELE_ACCEPT_TEST_H__

#include <limits>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

#include "pele/array.h"
#include "mc.h"

namespace mcpele{

/*Metropolis acceptance criterion
 * */

class MetropolisTest : public AcceptTest {
protected:
    size_t m_seed;
    std::mt19937_64 m_generator;
    std::uniform_real_distribution<double> m_distribution;
public:
    MetropolisTest(const size_t rseed);
    virtual ~MetropolisTest() {}
    virtual bool test(pele::Array<double> &trial_coords, double trial_energy,
            pele::Array<double> & old_coords, double old_energy, double temperature,
            MC * mc);
    size_t get_seed() const {return m_seed;}
    void set_generator_seed(const size_t inp) { m_generator.seed(inp); }
};

/*ENERGY WINDOW TEST
 * this test checks that the energy of the system stays within a certain energy window
 * */

class EnergyWindowTest : public AcceptTest {
protected:
    double m_min_energy;
    double m_max_energy;
public:
    EnergyWindowTest(const double min_energy, const double max_energy);
    virtual ~EnergyWindowTest() {}
    virtual bool test(pele::Array<double> &trial_coords, double trial_energy,
            pele::Array<double> & old_coords, double old_energy, double temperature,
            MC * mc);
};

} // namespace mcpele

#endif // #ifndef _MCPELE_ACCEPT_TEST_H__
