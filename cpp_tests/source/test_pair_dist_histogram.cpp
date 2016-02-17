#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include <gtest/gtest.h>

#include "pele/lbfgs.h"

#include "mcpele/record_pair_dist_histogram.h"

using pele::Array;

#define EXPECT_NEAR_RELATIVE(A, B, T)  EXPECT_NEAR(fabs(A)/(fabs(A)+fabs(B)+1), fabs(B)/(fabs(A)+fabs(B)+1), T)

using std::shared_ptr;
using mcpele::ConfTest;
using mcpele::MC;

#define boxdim 3

struct TrivialPotential : public pele::BasePotential{
    size_t call_count;
    virtual double get_energy(Array<double> coords)
    {
        call_count++;
        return 0.;
    }
};

struct TrivialTakestep : public mcpele::TakeStep{
    size_t call_count;
    TrivialTakestep()
        : call_count(0)
    {}
    virtual void displace(Array<double> &coords, MC * mc)
    {
        call_count++;
    }
};

class TestPairDistHist: public ::testing::Test{
public:

    typedef pele::Array<double> arr_t;
    size_t nparticles;
    size_t ndof;
    arr_t origin;
    arr_t x;
    double stepsize;
    std::shared_ptr<TrivialPotential> potential;
    std::shared_ptr<TrivialTakestep> step;
    size_t max_iter;
    size_t nr_bins;
    arr_t boxvector;
    double boxlength;
    virtual void SetUp(){
        nparticles = 1e1;
        ndof = boxdim * nparticles;
        origin = Array<double>(ndof);
        std::fill(origin.data(), origin.data() + ndof, 0);
        x = Array<double>(ndof);
        std::fill(x.data(), x.data() + ndof, 0);
        stepsize = 1e-2;
        potential = std::make_shared<TrivialPotential>();
        step = std::make_shared<TrivialTakestep>();
        max_iter = 1e5;
        nr_bins = 100;
        boxvector = arr_t(boxdim);
        boxlength = 42;
        std::fill(boxvector.data(), boxvector.data() + boxdim, boxlength);
    }
};

TEST_F(TestPairDistHist, BasicFunctionality){
    std::shared_ptr<mcpele::MC> mc = std::make_shared<mcpele::MC>(potential, x, 1);
    mc->set_takestep(step);
    const size_t eqsteps(max_iter/1e1);
    std::shared_ptr<mcpele::RecordPairDistHistogram<boxdim> > record_gr = std::make_shared<mcpele::RecordPairDistHistogram<boxdim> >(boxvector, nr_bins, eqsteps, 1);
    std::shared_ptr<pele::GradientOptimizer> opt = std::make_shared<pele::LBFGS>(potential, x);
    std::shared_ptr<mcpele::RecordPairDistHistogram<boxdim> > record_gr_quench = std::make_shared<mcpele::RecordPairDistHistogram<boxdim> >(boxvector, nr_bins, eqsteps, 1, opt);
    mc->add_action(record_gr);
    mc->add_action(record_gr_quench);
    mc->run(max_iter);
    EXPECT_EQ(mc->get_iterations_count(), max_iter);
    EXPECT_TRUE(record_gr->get_eqsteps() == eqsteps);
    EXPECT_TRUE(record_gr_quench->get_eqsteps() == eqsteps);
    const double number_density = nparticles / pow(boxlength, boxdim);
    pele::Array<double> r = record_gr->get_hist_r();
    pele::Array<double> gr = record_gr->get_hist_gr(number_density, nparticles);
    EXPECT_EQ(r[0], (boxlength * 0.5 / nr_bins) * 0.5); // bin size OK
    // normalization OK
    if (boxdim == 2) {
        EXPECT_DOUBLE_EQ(gr[0], nparticles * (nparticles - 1) / (nparticles * number_density * 2 * M_PI * r[0] * (boxlength * 0.5 / nr_bins)));
    }
    if (boxdim == 3) {
        const double dr = (boxlength * 0.5 / nr_bins);
        EXPECT_DOUBLE_EQ(gr[0], nparticles * (nparticles - 1) / (nparticles * number_density * 4 * M_PI / 3 * (pow(r[0] + 0.5 * dr, 3) - pow(r[0] - 0.5 * dr, 3))));
    }
    for (size_t i = 1; i < gr.size(); ++i) {
        EXPECT_EQ(gr[i], 0);
    }
}

