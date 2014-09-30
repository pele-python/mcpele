#include "energy_window_test.h"

#include "pele/array.h"
#include "mc.h"

using pele::Array;

namespace mcpele {

EnergyWindowTest::EnergyWindowTest(const double min_energy, const double max_energy)
    : m_min_energy(min_energy),
      m_max_energy(max_energy)
{}

bool EnergyWindowTest::test(Array<double> &trial_coords, double trial_energy,
    Array<double> & old_coords, double old_energy, double temperature,
    MC * mc)
{
    return ((trial_energy >= m_min_energy) and (trial_energy <= m_max_energy));
}

} // namespace mcpele
