#include "mcpele/metropolis_test.h"

#include <cmath>
//#include <chrono>

using pele::Array;

namespace mcpele {

MetropolisTest::MetropolisTest(const size_t rseed)
    : m_seed(rseed),
      m_generator(rseed),
      m_distribution(0.0, 1.0)
{
    #ifdef DEBUG
        std::cout << "seed Metropolis:" << _seed << "\n";
        //std::chrono::system_clock::now().time_since_epoch().count()
    #endif
}

bool MetropolisTest::test(Array<double> &trial_coords, double trial_energy,
        Array<double>& old_coords, double old_energy, double temperature,
        MC * mc)
{
    bool success = true;
    double dE = trial_energy - old_energy;
    if (dE > 0.){
        double w = exp(-dE / temperature);
        double rand = m_distribution(m_generator);
        if (rand > w) {
            success = false;
        }
    }
    return success;
}

} // namespace mcpele
