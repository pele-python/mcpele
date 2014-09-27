#ifndef _MCPELE_CONF_TEST_H__
#define _MCPELE_CONF_TEST_H__

#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <chrono>

#include "pele/array.h"
#include "pele/optimizer.h"
#include "pele/distance.h"
#include "mc.h"

namespace mcpele{

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

#endif // #ifndef _MCPELE_CONF_TEST_H__
