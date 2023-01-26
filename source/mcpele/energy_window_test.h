#ifndef _MCPELE_ENERGY_WINDOW_TEST_H__
#define _MCPELE_ENERGY_WINDOW_TEST_H__

#include "mc.h"
#include "pele/array.hpp"

namespace mcpele {

/**
 * Energy window test
 * This test checks that the energy of the system stays within a certain energy
 * window
 */
class EnergyWindowTest : public AcceptTest {
 protected:
  double m_min_energy;
  double m_max_energy;

 public:
  EnergyWindowTest(const double min_energy, const double max_energy);
  virtual ~EnergyWindowTest() {}
  virtual bool test(pele::Array<double> &trial_coords, double trial_energy,
                    pele::Array<double> &old_coords, double old_energy,
                    double temperature, MC *mc);
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_ENERGY_WINDOW_TEST_H__
