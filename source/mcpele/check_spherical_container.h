#ifndef _MCPELE_CHECK_SPHERICAL_CONTAINER_H__
#define _MCPELE_CHECK_SPHERICAL_CONTAINER_H__

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>

#include "mc.h"
#include "pele/array.hpp"
#include "pele/distance.hpp"
#include "pele/optimizer.hpp"

namespace mcpele {

class CheckSphericalContainer : public ConfTest {
 protected:
  double m_radius2;
  size_t m_ndim;

 public:
  CheckSphericalContainer(const double radius, const size_t ndim);
  virtual bool conf_test(pele::Array<double> &trial_coords, MC *mc);
  virtual ~CheckSphericalContainer() {}
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_CHECK_SPHERICAL_CONTAINER_H__
