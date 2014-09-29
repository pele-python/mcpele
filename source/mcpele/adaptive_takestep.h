#ifndef _MCPELE_ADAPTIVE_TAKESTEP_H__
#define _MCPELE_ADAPTIVE_TAKESTEP_H__

#include "mcpele/mc.h"

namespace mcpele {

class AdaptiveTakeStep : public TakeStep {
protected:
    std::shared_ptr<TakeStep> m_ts;
    size_t m_interval;
    size_t m_total_steps;
    size_t m_accepted_steps;
public:
    virtual ~AdaptiveTakeStep() {}
    AdaptiveTakeStep(std::shared_ptr<TakeStep> ts, const size_t interval=100)
        : m_ts(ts),
          m_interval(interval),
          m_total_steps(0),
          m_accepted_steps(0)
    {}
    void displace(pele::Array<double> &coords, MC * mc)
    {
        m_ts->displace(coords, mc);
    }
    void report(const MC* mc)
    {
        const double success = mc->get_success();
        ++m_total_steps;
        if (success) {
            ++m_accepted_steps;
        }
        if (mc->get_iterations_count() % m_interval == 0) {
            const double acceptance_fraction = static_cast<double>(m_accepted_steps) / static_cast<double>(m_total_steps);
            m_accepted_steps = 0;
            m_total_steps = 0;
            if (acceptance_fraction < m_ts->get_min_acceptance_ratio()) {
                m_ts->increase_acceptance();
            }
            else if (acceptance_fraction > m_ts->get_max_acceptance_ratio()) {
                m_ts->decrease_acceptance();
            }
        }
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_ADAPTIVE_TAKESTEP_H__
