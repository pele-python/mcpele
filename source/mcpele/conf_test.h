#ifndef _MCPELE_CONF_TEST_H__
#define _MCPELE_CONF_TEST_H__

#include <iostream>
#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "pele/array.h"
#include "pele/optimizer.h"
#include "pele/distance.h"
#include "mc.h"

using std::runtime_error;
using pele::Array;

namespace mcpele{

class CheckSphericalContainer:public ConfTest{
protected:
    double _radius2;
    size_t _ndim;
public:

    CheckSphericalContainer(double radius, size_t ndim);
    virtual bool test(Array<double> &trial_coords, MC * mc);
    virtual ~CheckSphericalContainer(){}
};

}//namespace mcpele


#endif
