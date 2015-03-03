#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <gtest/gtest.h>

#include "pele/array.h"

#include "mcpele/random_coords_displacement.h"
#include "mcpele/gaussian_coords_displacement.h"
#include "mcpele/histogram.h"

using pele::Array;

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

class TestHistogram: public ::testing::Test{
public:
    typedef pele::Array<double> arr_t;
    size_t ndof;
    size_t nsteps;
    size_t ntot;
    arr_t displ_uniform_all, displ_uniform_single;
    arr_t displ_gaussian;
    arr_t baseline;
    double ss;
    mcpele::MC* mc;
    virtual void SetUp(){
        ndof = 1000;
        nsteps = 100;
        ntot = nsteps*ndof;
        displ_uniform_all = Array<double>(ndof);
        displ_uniform_single = Array<double>(ndof);
        displ_gaussian = Array<double>(ndof);
        baseline = Array<double>(ndof);
        for (size_t i = 0; i < ndof; ++i){
            displ_uniform_all[i] = 0;
            displ_uniform_single[i] = 0;
            displ_gaussian[i] = 0;
            baseline[i] = 0;
        }
        ss = 4.2;
        mc = NULL;
    }
};

TEST_F(TestHistogram, TestMomentsGlobalMoves){
    std::fill(displ_gaussian.data(), displ_gaussian.data() + ndof, 0);
    mcpele::SampleGaussian sampler_gaussian(42, ss, displ_gaussian);
    mcpele::RandomCoordsDisplacementAll sampler_uniform_all(42, ss);
    mcpele::Histogram hist_uniform_all(-0.5 * ss, 0.5 * ss, 0.01);
    mcpele::Histogram hist_gaussian(-5, 5, 0.1);
    for (size_t step = 0; step < nsteps; ++step) {
        std::fill(displ_uniform_all.data(), displ_uniform_all.data() + ndof, 0);
        std::fill(displ_gaussian.data(), displ_gaussian.data() + ndof, 0);
        sampler_uniform_all.displace(displ_uniform_all, mc);
        sampler_gaussian.displace(displ_gaussian, mc);
        for (size_t dof = 0; dof < ndof; ++dof) {
            hist_uniform_all.add_entry(displ_uniform_all[dof]);
            hist_gaussian.add_entry(displ_gaussian[dof]);
        }
    }
    EXPECT_EQ(static_cast<size_t>(hist_uniform_all.entries()), ntot);
    EXPECT_EQ(static_cast<size_t>(hist_gaussian.entries()), ntot);
    EXPECT_NEAR_RELATIVE(hist_uniform_all.get_mean(), sampler_uniform_all.expected_mean(), 2 * ss / sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_gaussian.get_mean(), sampler_gaussian.expected_mean(), 2 * ss / sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_uniform_all.get_variance(), sampler_uniform_all.expected_variance(ss), 2 * ss / sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_gaussian.get_variance(), sampler_gaussian.expected_variance(ss), 2 * ss / sqrt(ntot));
}

TEST_F(TestHistogram, TestMomentsSingleMoves){
    mcpele::RandomCoordsDisplacementSingle sampler_uniform_single(42, 500, 2, ss);
    mcpele::Histogram hist_uniform_single(-0.5 * ss, 0.5 * ss, 0.01);
    for (size_t step = 0; step < nsteps*500; ++step) {
        std::fill(displ_uniform_single.data(), displ_uniform_single.data() + ndof, 0);
        sampler_uniform_single.displace(displ_uniform_single, mc);
        for (size_t dof = 0; dof < ndof; ++dof) {
            if (displ_uniform_single[dof]!= 0){
                hist_uniform_single.add_entry(displ_uniform_single[dof]);
            }
        }
    }

    EXPECT_EQ(static_cast<size_t>(hist_uniform_single.entries()), ntot);
    EXPECT_NEAR_RELATIVE(hist_uniform_single.get_mean(), sampler_uniform_single.expected_mean(), 2 * ss / sqrt(ntot));
    EXPECT_NEAR_RELATIVE(hist_uniform_single.get_variance(), sampler_uniform_single.expected_variance(ss), 2 * ss / sqrt(ntot));
}

TEST_F(TestHistogram, TestBinning){
    std::fill(displ_gaussian.data(), displ_gaussian.data() + ndof, 0);
    mcpele::SampleGaussian sampler(42, ss, displ_gaussian);
    const double min = -42;
    const double max = 42;
    const double bin = 2;
    mcpele::Histogram hist(min, max, bin);
    for (size_t step = 0; step < nsteps; ++step) {
        sampler.displace(displ_gaussian, mc);
        for (size_t dof = 0; dof < ndof; ++dof) {
            hist.add_entry(displ_gaussian[dof]);
        }
    }
    EXPECT_EQ(static_cast<size_t>(hist.entries()), ntot);
    EXPECT_LE(min, hist.min());
    EXPECT_LE(max, hist.max());
    EXPECT_DOUBLE_EQ(hist.bin(), bin);
    EXPECT_DOUBLE_EQ(hist.entries(), nsteps * ndof);
    EXPECT_DOUBLE_EQ(std::accumulate(hist.begin(), hist.end(), double(0)), hist.entries());
    const std::vector<double> vecdata = hist.get_vecdata();
    EXPECT_EQ(vecdata.size(), hist.size());
    EXPECT_DOUBLE_EQ(std::accumulate(vecdata.begin(), vecdata.end(), double(0)), hist.entries());
    const std::vector<double> vecdata_normalized = hist.get_vecdata_normalized();
    EXPECT_EQ(vecdata_normalized.size(), hist.size());
    EXPECT_DOUBLE_EQ(std::accumulate(vecdata_normalized.begin(), vecdata_normalized.end(), double(0)) * bin, 1);
    const std::vector<double> error = hist.get_vecdata_error();
    for (size_t ii = 0; ii < hist.size(); ++ii) {
        const double xi = hist.get_position(ii);
        const double true_i = exp(-pow((xi - sampler.expected_mean()), 2) / (2 * sampler.expected_variance(ss))) / sqrt(M_PI * 2 * sampler.expected_variance(ss));
        const double EPS = 1e-12;
        if (vecdata_normalized.at(ii) > EPS) {
            EXPECT_NEAR(vecdata_normalized.at(ii), true_i, 2 * error.at(ii));
        }
    }
}
