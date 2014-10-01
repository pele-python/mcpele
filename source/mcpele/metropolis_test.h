#ifndef _MCPELE_METROPOLIS_TEST_H__
#define _MCPELE_METROPOLIS_TEST_H__

#include <random>

#include "pele/array.h"
#include "mc.h"

namespace mcpele {

/**
 * Metropolis acceptance criterion
 */
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

} // namespace mcpele

#endif // #ifndef _MCPELE_METROPOLIS_TEST_H__
