#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "pele/array.h"

#include "mcpele/conf_test.h"

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

TEST(CheckSphericalContainer, Works){
    size_t ndim = 4;
    pele::Array<double> x(2*ndim,0);
    double r = 2.;
    double eps = 1e-10;
    mcpele::MC * mcvoid = NULL;

    mcpele::CheckSphericalContainer check(r, ndim);
    x[0] = r - eps;
    EXPECT_TRUE(check.conf_test(x, mcvoid));
    x[0] = r + eps;
    EXPECT_FALSE(check.conf_test(x, mcvoid));
}
