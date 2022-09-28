#ifndef _MCPELE_CHECK_SPHERICAL_CONTAINER_H__
#define _MCPELE_CHECK_SPHERICAL_CONTAINER_H__

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

#include "pele/array.hpp"
#include "pele/optimizer.hpp"
#include "pele/distance.hpp"

#include "mc.h"

namespace mcpele {

class CheckSphericalContainer : public ConfTest {
protected:
    double m_radius2;
    size_t m_ndim;
public:
    CheckSphericalContainer(const double radius, const size_t ndim);
    virtual bool conf_test(pele::Array<double> &trial_coords, MC * mc);
    virtual ~CheckSphericalContainer() {}
};

} // namespace mcpele

#endif // #ifndef _MCPELE_CHECK_SPHERICAL_CONTAINER_H__
