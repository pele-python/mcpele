#include <stdexcept>
#include <cmath>

#include "actions.h"

using pele::Array;
using std::sqrt;

namespace mcpele{

/*
 * Record energy histogram
*/

RecordEnergyHistogram::RecordEnergyHistogram(const double min, const double max, const double bin, const size_t eqsteps)
    : m_hist(min, max, bin),
      m_eqsteps(eqsteps),
      m_count(0)
    {}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) 
{
    m_count = mc->get_iterations_count();
    if (m_count > m_eqsteps){
        m_hist.add_entry(energy);
    }
}

/*
 * Record scalar time series, every record_every-th step.
 */

RecordScalarTimeseries::RecordScalarTimeseries(const size_t niter, const size_t record_every)
    : m_niter(niter),
      m_record_every(record_every)
{
    if (record_every==0) {
        throw std::runtime_error("RecordScalarTimeseries: record_every expected to be at least 1");
    }
    m_time_series.reserve(niter/record_every);
}

void RecordScalarTimeseries::action(Array<double> &coords, double energy, bool accepted, MC* mc)
{
    const size_t counter = mc->get_iterations_count();
    if (counter % m_record_every == 0) {
        m_record_scalar_value(this->get_recorded_scalar(coords, energy, accepted, mc));
    }
}

/*
 * Record energy time series, measuring every __record_every-th step.
 */

RecordEnergyTimeseries::RecordEnergyTimeseries(const size_t niter, const size_t record_every)
    : RecordScalarTimeseries(niter, record_every)
{}

double RecordEnergyTimeseries::get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)
{
    return energy;
}

/*
 * Record time series of lowest eigenvalue
 */

RecordLowestEValueTimeseries::RecordLowestEValueTimeseries(const size_t niter, const size_t record_every,
        std::shared_ptr<pele::BasePotential> landscape_potential, const size_t boxdimension,
        pele::Array<double> ranvec, const size_t lbfgsniter)
    : RecordScalarTimeseries(niter, record_every),
      m_lowest_ev(landscape_potential, boxdimension, ranvec, lbfgsniter)
{}

double RecordLowestEValueTimeseries::get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)
{
    return m_lowest_ev.compute_lowest_eigenvalue(coords);
}

/*
 * Record time series of root mean squared displacement (averaged over all particles)
 * Motivation: check if HS fluid is decorrelated between snapshots
 */

RecordDisplacementPerParticleTimeseries::RecordDisplacementPerParticleTimeseries(const size_t niter, const size_t record_every,
            pele::Array<double> initial_coords, const size_t boxdimension)
    : RecordScalarTimeseries(niter, record_every),
      m_rsm_displacement(initial_coords, boxdimension)
{}

double RecordDisplacementPerParticleTimeseries::get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)
{
    return m_rsm_displacement.compute_mean_particle_displacement(coords);
}

}//namespace mcpele
