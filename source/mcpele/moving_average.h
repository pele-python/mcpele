#ifndef _MCPELE_MOVING_AVERAGE_H
#define _MCPELE_MOVING_AVERAGE_H

#include <vector>
#include <algorithm>

#include "histogram.h"

namespace mcpele {

/**
 * Computes moving averages of time series.
 * Refernece to time series, from small to large times, is in m_time_series.
 * Assume that one wants to compute moving averages for the rightmost (latest)
 * m_nr_steps_total elements of the time series.
 * The number of time series steps that go into one moving average is
 * m_nr_steps_ma.
 * The number of different moving averages that one can compute, given these
 * two parameters (m_nr_steps_total and m_nr_steps_ma) is m_nr_different_ma.
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
    const size_t m_nr_steps_ma;
    const size_t m_nr_different_ma;
    std::vector<double>::const_iterator m_begin;
    std::vector<double>::const_iterator m_end;
    double m_this_ma;
public:
    MovingAverageAcc(const std::vector<double>& time_series, const size_t nr_steps_total, const size_t nr_steps_ma)
        : m_time_series(time_series),
          m_nr_steps_total(nr_steps_total),
          m_nr_steps_ma(nr_steps_ma),
          m_nr_different_ma(nr_steps_total - nr_steps_ma + 1),
          m_begin(m_time_series.end() - nr_steps_total),
          m_end(m_begin + m_nr_steps_ma),
          m_this_ma(std::accumulate(m_begin, m_end, 0) / static_cast<double>(m_nr_steps_ma))
    {
        if (nr_steps_ma % 2 != 0) {
            throw std::runtime_error("MovingAverageAcc: illegal input: nr_steps_ma");
        }
        if (time_series.size() < nr_steps_total) {
            throw std::runtime_error("MovingAverageAcc: illegal input: time series too short");
        }
    }
    double get() const
    {
        return m_this_ma;
    }
    size_t get_nr_different_ma() const
    {
        return m_nr_different_ma;
    }
    void shift_right()
    {
        ++m_begin;
        ++m_end;
        if (m_end == m_time_series.end()) {
            reset();
        }
        else {
            m_this_ma -= *(m_begin - 1) / static_cast<double>(m_nr_steps_ma);
            m_this_ma += *(m_end - 1) / static_cast<double>(m_nr_steps_ma);
        }
    }
    void reset()
    {
        m_begin = m_time_series.end() - m_nr_steps_total;
        m_end = m_begin + m_nr_steps_ma;
        m_this_ma = std::accumulate(m_begin, m_end, 0) / static_cast<double>(m_nr_steps_ma);
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_MOVING_AVERAGE_H
