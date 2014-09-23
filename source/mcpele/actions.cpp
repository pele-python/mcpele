#include <stdexcept>
#include <cmath>
#include <algorithm>

#include "actions.h"
#include "moving_average.h"

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
    : m_target(target),
      m_factor(factor),
      m_acceptedf(0),
      m_niter(niter),
      m_navg(navg),
      m_count(0),
      m_naccepted(0),
      m_nrejected(0)
{
    if (factor > 1 || factor < 0) 
        throw std::runtime_error("AdjustStep: illegal input");
}

void AdjustStep::action(Array<double> &coords, double energy, bool accepted, MC* mc) 
{

    m_count = mc->get_iterations_count();

    if (m_count <= m_niter) {
        if (accepted == true)
            ++m_naccepted;
        else
            ++m_nrejected;

        if(m_count % m_navg == 0)
        {
            m_acceptedf = (double) m_naccepted / (m_naccepted + m_nrejected);

            //std::cout<<"acceptance "<<_acceptedf<< "\n";
            //std::cout<<"stepsize before"<<mc->_stepsize<< "\n";
            if (m_acceptedf < m_target)
                mc->_stepsize *= m_factor;
            else
                mc->_stepsize /= m_factor;
            //std::cout<<"stepsize after"<<mc->_stepsize<< "\n";

            //now reset to zero memory of acceptance and rejection
            m_naccepted = 0;
            m_nrejected = 0;
        }

    }
}


/*
 * Record energy histogram
*/

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps)
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
    if (record_every == 0) {
        throw std::runtime_error("RecordScalarTimeseries: record_every expected to be at least 1");
    }
    m_time_series.reserve(niter / record_every);
}

void RecordScalarTimeseries::action(Array<double> &coords, double energy, bool accepted, MC* mc)
{
    const size_t counter = mc->get_iterations_count();
    if (counter % m_record_every == 0) {
        m_record_scalar_value(this->get_recorded_scalar(coords, energy, accepted, mc));
    }
}

/**
 * Apply a basic check on the behaviour of the moving average of the recorded
 * time series to decide if it is sufficiently stable.
 * Parameters:
 * nr_steps_to_check is the number of steps, from the present into the past,
 * that are considered for the moving average check.
 * rel_std_threshold is the relative size of the fluctuations in the moving
 * average that is considered stable.
 * Decreasing that number should make the test more sensitive.
 * NOTE: we do not test the variance because the variance has some functional
 * dependence on the mean var = f(mu). In order to remove this one would require
 * a "variance stabilising transformation" which would require the implementation
 * of a different moments class where f(x) si stored giving an approximately constant
 * variance as a function of the mean. See Ascombe transform for instance.
 */
bool RecordScalarTimeseries::moving_average_is_stable(const size_t nr_steps_to_check, const double rel_std_threshold)
{
    std::pair<double,double> mean_ma = this->get_moving_average_mean(nr_steps_to_check);
    const double mu_mean_ma = mean_ma.first;
    const double std_mean_ma = mean_ma.second;
    //std::cout<<"std_mean_ma/mu_mean_ma "<<100*std_mean_ma/mu_mean_ma<<std::endl;

    return ( (std_mean_ma / mu_mean_ma) < rel_std_threshold);
    /*std::pair<double,double> var_ma = this->get_moving_average_variance(nr_steps_to_check);
    const double mu_var_ma = var_ma.first;
    const double std_var_ma = var_ma.second;
    std::cout<<"std_var_ma/mu_var_ma "<<100*std_var_ma/mu_var_ma<<std::endl;
    && (std_var_ma / mu_var_ma) < rel_std_threshold*/
}

std::pair<double,double> RecordScalarTimeseries::get_moving_average_mean(const size_t nr_steps_to_check)
{
    const size_t scale = 10;
    const size_t tmp = nr_steps_to_check / scale;
    const size_t window_size = tmp + (tmp % 2);
    MovingAverageAcc moving_average(m_time_series, nr_steps_to_check, window_size);
    Moments moving_average_mean;
    const size_t nr_steps_ma = moving_average.get_nr_steps_ma();
    for (size_t i = 0; i < nr_steps_ma; ++i, moving_average.shift_right()) {
        moving_average_mean.update(moving_average.get_mean());
    }
    return std::pair<double,double>(moving_average_mean.mean(),moving_average_mean.std());
}

std::pair<double,double> RecordScalarTimeseries::get_moving_average_variance(const size_t nr_steps_to_check)
{
    const size_t scale = 10;
    const size_t tmp = nr_steps_to_check / scale;
    const size_t window_size = tmp + (tmp % 2);
    MovingAverageAcc moving_average(m_time_series, nr_steps_to_check, window_size);
    Moments moving_average_var;
    const size_t nr_steps_ma = moving_average.get_nr_steps_ma();
    for (size_t i = 0; i < nr_steps_ma; ++i, moving_average.shift_right()) {
        moving_average_var.update(moving_average.get_variance());
    }
    return std::pair<double,double>(moving_average_var.mean(),moving_average_var.std());
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

RecordMeanRMSDisplacementTimeseries::RecordMeanRMSDisplacementTimeseries(const size_t niter, const size_t record_every,
            pele::Array<double> initial_coords, const size_t boxdimension)
    : RecordScalarTimeseries(niter, record_every),
      m_rsm_displacement(initial_coords, boxdimension)
{}

double RecordMeanRMSDisplacementTimeseries::get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)
{
    return m_rsm_displacement.compute_mean_rsm_displacement(coords);
}

}//namespace mcpele
