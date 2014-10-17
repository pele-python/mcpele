#include "pele/array.h"

#include "mcpele/energy_window_test.h"

using pele::Array;

namespace mcpele {

EnergyWindowTest::EnergyWindowTest(const double min_energy, const double max_energy)
    : m_min_energy(min_energy),
      m_max_energy(max_energy)
{}

bool EnergyWindowTest::test(Array<double>&, double trial_energy,
        Array<double>&, double, double, MC*)
{
    return ((trial_energy >= m_min_energy) and (trial_energy <= m_max_energy));
}

} // namespace mcpele
