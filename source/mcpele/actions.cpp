#include <stdexcept>
#include <cmath>

#include "actions.h"

using pele::Array;
using std::sqrt;

namespace mcpele{

/*Adjust Step
 *     factor is a multiplicative factor by which the stepsize is adjusted
 *     niter determines the number of steps for which the action should take effect (generally
 *     we want to adjust the step size only at the beginning of a simulation)
 *     navg is the number of steps over which the acceptance is averaged
 *     factor must be 0<f<1, if rejected make step shorter, if accepted make step longer
*/

AdjustStep::AdjustStep(double target, double factor, size_t niter, size_t navg)
    : _target(target), _factor(factor), _acceptedf(0),
    _niter(niter), _navg(navg), _count(0), _naccepted(0),
    _nrejected(0)
{
    if (factor > 1 || factor < 0) 
        throw std::runtime_error("AdjustStep: illegal input");
}

void AdjustStep::action(Array<double> &coords, double energy, bool accepted, MC* mc) 
{

    _count = mc->get_iterations_count();

    if (_count <= _niter) {
        if (accepted == true)
            ++_naccepted;
        else
            ++_nrejected;

        if(_count % _navg == 0)
        {
            _acceptedf = (double) _naccepted / (_naccepted + _nrejected);

            //std::cout<<"acceptance "<<_acceptedf<< "\n";
            //std::cout<<"stepsize before"<<mc->_stepsize<< "\n";
            if (_acceptedf < _target)
                mc->_stepsize *= _factor;
            else
                mc->_stepsize /= _factor;
            //std::cout<<"stepsize after"<<mc->_stepsize<< "\n";

            //now reset to zero memory of acceptance and rejection
            _naccepted = 0;
            _nrejected = 0;
        }

    }
}


/*
 * Record energy histogram
*/

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double
        bin, size_t eqsteps)
    : _hist(mcpele::Histogram(min, max, bin)), 
    _eqsteps(eqsteps), 
    _count(0)
{}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) 
{
    _count = mc->get_iterations_count();
    if (_count > _eqsteps){
        _hist.add_entry(energy);
    }
}

/*
 * Record scalar time series, every record_every-th step.
 */

RecordScalarTimeseries::RecordScalarTimeseries(const size_t niter, const size_t record_every)
    : _niter(niter),
      _record_every(record_every)
{
    if (record_every==0) {
        throw std::runtime_error("RecordScalarTimeseries: record_every expected to be at least 1");
    }
    _time_series.reserve(niter/record_every);
}

void RecordScalarTimeseries::action(Array<double> &coords, double energy, bool accepted, MC* mc)
{
    const size_t counter = mc->get_iterations_count();
    if (counter % _record_every == 0) {
        _record_scalar_value(this->get_recorded_scalar(coords, energy, accepted, mc));
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
      _lowest_ev(landscape_potential, boxdimension, ranvec, lbfgsniter)
{}

double RecordLowestEValueTimeseries::get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)
{
    return _lowest_ev.compute_lowest_eigenvalue(coords);
}

/*
 * Record time series of root mean squared displacement (averaged over all particles)
 * Motivation: check if HS fluid is decorrelated between snapshots
 */

RecordMeanRMSDisplacementTimeseries::RecordMeanRMSDisplacementTimeseries(const size_t niter, const size_t record_every,
            pele::Array<double> initial_coords, const size_t boxdimension)
    : RecordScalarTimeseries(niter, record_every),
      _rsm_displacement(initial_coords, boxdimension)
{}

double RecordMeanRMSDisplacementTimeseries::get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)
{
    return _rsm_displacement.compute_mean_rsm_displacement(coords);
}

}//namespace mcpele
