#include <cmath>
#include <random>
#include <gtest/gtest.h>

#include "mcpele/cdf_accumulator.h"

class TestCDFAcc : public ::testing::Test {
public:
    size_t nr_points;
    std::mt19937_64 generator;
    double low;
    double high;
    std::uniform_real_distribution<double> distribution;
    void SetUp()
    {
        nr_points = 100000;
        generator = std::mt19937_64(42);
        low = - 42;
        high = 42.42;
        distribution = std::uniform_real_distribution<double>(low, high);
    }
};

TEST_F(TestCDFAcc, Works)
{
    mcpele::CDFAccumulator cdf;
    for (size_t i = 0; i < nr_points; ++i) {
        cdf.add(distribution(generator));
    }
    const auto x = cdf.get_vecdata_x();
    const auto cdf_x = cdf.get_vecdata_cdf_x();
    EXPECT_EQ(x.size(), cdf_x.size());
    EXPECT_LE(x.size(), nr_points);
    EXPECT_LE(low, x.front());
    EXPECT_LE(x.back(), high);
    EXPECT_DOUBLE_EQ(cdf_x.front(), 1);
    EXPECT_DOUBLE_EQ(cdf_x.back(), double(1) / double(nr_points));
    for (size_t i = 0; i < x.size(); ++i) {
        EXPECT_NEAR(cdf_x.at(i), (high - x.at(i)) / (high - low), 1 / sqrt(nr_points));
    }
}
