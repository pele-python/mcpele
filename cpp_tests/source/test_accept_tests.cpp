#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>

#include "pele/array.h"

#include "mcpele/accept_test.h"

TEST(EnergyWindow, Works){
    double emin = 1.;
    double emax = 2.;
    pele::Array<double> xnew, xold;
    mcpele::MC * mc = NULL;
    double T = 0., eold = 0.;
    double eps = 1e-10;
    mcpele::EnergyWindowTest test(emin, emax);
    EXPECT_TRUE(test.test(xnew, emin+eps, xold, eold, T, mc));
    EXPECT_TRUE(test.test(xnew, emax-eps, xold, eold, T, mc));
    EXPECT_FALSE(test.test(xnew, emin-eps, xold, eold, T, mc));
    EXPECT_FALSE(test.test(xnew, emax+eps, xold, eold, T, mc));
}
