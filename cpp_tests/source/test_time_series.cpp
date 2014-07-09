#include <iostream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <memory>
#include <gtest/gtest.h>

#include "pele/array.h"
#include "pele/harmonic.h"

#include "mcpele/actions.h"

struct TrivialTakestep : public mcpele::TakeStep{
    size_t call_count;
    TrivialTakestep()
        : call_count(0)
    {}
    virtual void takestep(pele::Array<double> &coords, double stepsize, mcpele::MC * mc=NULL)
    {
        call_count++;
    }
};

TEST(EnergyTimeseries, Basic){
    const size_t boxdim = 3;
    const size_t nparticles = 100;
    const size_t ndof = nparticles*boxdim;
    const size_t niter = 10000;
    const size_t record_every = 100;
    const double k = 400;
    const double stepsize = 1e-2;
    pele::Array<double> coords(ndof,1);
    pele::Array<double> origin(ndof,0);
    pele::Harmonic* potential = new pele::Harmonic(origin, k, boxdim);
    const auto enumerical = potential->get_energy(coords);
    double etrue(0);
    for (size_t i = 0; i < ndof; ++i) {
        const auto delta = coords[i] - origin[i];
        etrue += 0.5*k*delta*delta;
    }
    EXPECT_EQ(enumerical, etrue);
    std::shared_ptr<mcpele::MC> mc = std::make_shared<mcpele::MC>(potential, coords, 1, stepsize);
    mcpele::RecordEnergyTimeseries* ts = new mcpele::RecordEnergyTimeseries(niter, record_every);
    mc->add_action(std::shared_ptr<mcpele::RecordEnergyTimeseries>(ts));
    mc->set_takestep(std::make_shared<TrivialTakestep>());
    mc->run(niter);
    EXPECT_EQ(mc->get_iterations_count(), niter);
    pele::Array<double> series = ts->get_time_series();
    EXPECT_EQ(series.size(), niter/record_every);
    for (size_t i = 0; i < series.size(); ++i) {
        EXPECT_EQ(series[i], enumerical);
    }
    delete potential;
}

TEST(EVTimeseries, Works){
    const size_t boxdim = 3;
    const size_t nparticles = 100;
    const size_t ndof = nparticles*boxdim;
    const size_t niter = 10000;
    const size_t record_every = 100;
    const double k = 400;
    const double stepsize = 1e-2;
    pele::Array<double> coords(ndof,1);
    pele::Array<double> origin(ndof,0);
    pele::Harmonic* potential = new pele::Harmonic(origin, k, boxdim);
    const auto enumerical = potential->get_energy(coords);
    double etrue(0);
    for (size_t i = 0; i < ndof; ++i) {
        const auto delta = coords[i] - origin[i];
        etrue += 0.5*k*delta*delta;
    }
    EXPECT_EQ(enumerical, etrue);
    std::shared_ptr<mcpele::MC> mc = std::make_shared<mcpele::MC>(potential, coords, 1, stepsize);
    mcpele::RecordLowestEValueTimeseries* ts = new mcpele::RecordLowestEValueTimeseries(niter, record_every);
    mc->add_action(std::shared_ptr<mcpele::RecordLowestEValueTimeseries>(ts));
    mc->set_takestep(std::make_shared<TrivialTakestep>());
    mc->run(niter);
    EXPECT_EQ(mc->get_iterations_count(), niter);
    pele::Array<double> series = ts->get_time_series();
    EXPECT_EQ(series.size(), niter/record_every);
    for (size_t i = 0; i < series.size(); ++i) {
        EXPECT_EQ(series[i], enumerical);
    }
    delete potential;
}
