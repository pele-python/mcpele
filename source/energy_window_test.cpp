#include "pele/array.h"

#include "mcpele/energy_window_test.h"

using pele::Array;

namespace mcpele {

EnergyWindowTest::EnergyWindowTest(const double min_energy, const double max_energy)
    : m_min_energy(min_energy),
      m_max_energy(max_energy)
{}

bool EnergyWindowTest::test(__attribute__((unused)) Array<double> &trial_coords, double trial_energy,
        __attribute__((unused)) Array<double> & old_coords, __attribute__((unused)) double old_energy,
        __attribute__((unused)) double temperature, __attribute__((unused)) MC * mc)
{
    return ((trial_energy >= m_min_energy) and (trial_energy <= m_max_energy));
}

} // namespace mcpele
