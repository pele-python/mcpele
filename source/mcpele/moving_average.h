#ifndef _MCPELE_MOVING_AVERAGE_H
#define _MCPELE_MOVING_AVERAGE_H

#include <algorithm>
#include <vector>

#include "histogram.h"

namespace mcpele {

/**
 * Computes moving averages of time series.
 * Reference to time series, from small to large times, is in m_time_series.
 * Assume that one wants to compute moving averages for the rightmost (latest)
 * m_nr_steps_total elements of the time series.
 * The number of time series steps that go into one moving average is
 * m_nr_steps_ma.
 * The number of different moving averages that one can compute, given these
 * two parameters (m_nr_steps_total and m_window_size) is m_nr_steps_ma.
 * Moving average is initialised to the leftmost moving average.
 * After that shift_right() moves the moving average window one time series
 * step further to the right (to the future).
 * In case shift_right() reaches the last step in the time series, the moving
 * average window returns to the initial position.
 */
class MovingAverageAcc {
private:
    const std::vector<double>& m_time_series;
    const size_t m_nr_steps_total;
    const size_t m_window_size;
    const size_t m_nr_steps_ma;
    std::vector<double>::const_iterator m_begin;
    std::vector<double>::const_iterator m_end;
    mcpele::Moments m_moments;
public:
    MovingAverageAcc(const std::vector<double>& time_series, const size_t nr_steps_total, const size_t nr_steps_ma)
        : m_time_series(time_series),
          m_nr_steps_total(nr_steps_total),
          m_window_size(nr_steps_ma), //window size
          m_nr_steps_ma(nr_steps_total - m_window_size + 1), //number of steps to move window from left to right end
          m_begin(m_time_series.end() - nr_steps_total),
          m_end(m_begin + m_window_size),
          m_moments()
    {
        if (nr_steps_ma % 2 != 0) {
            throw std::runtime_error("MovingAverageAcc: illegal input: nr_steps_ma");
        }
        if (time_series.size() < nr_steps_total) {
            throw std::runtime_error("MovingAverageAcc: illegal input: time series too short");
        }
        //initialise moments
        for(auto it = m_begin; it != m_end; ++it) {
            m_moments(*it);
        }
    }
    double get_mean() const { return m_moments.mean(); }
    double get_variance() const { return m_moments.variance(); }
    size_t get_nr_steps_ma() const { return m_nr_steps_ma; }
    void shift_right()
    {
        ++m_begin;
        ++m_end;
        if (m_end == m_time_series.end()) {
            reset();
        }
        else {
            m_moments.replace(*(m_begin - 1), *(m_end - 1));
        }
    }
    void reset()
    {
        m_begin = m_time_series.end() - m_nr_steps_total;
        m_end = m_begin + m_window_size;
        m_moments = mcpele::Moments();
        //initialise moments
        for(auto it = m_begin; it != m_end; ++it) {
            m_moments(*it);
        }
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_MOVING_AVERAGE_H
