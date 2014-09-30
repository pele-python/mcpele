#include "record_scalar_timeseries.h"

using pele::Array;

namespace mcpele {

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

} // namespace mcpele
