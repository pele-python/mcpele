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
    AdaptiveTakeStep(std::shared_ptr<TakeStep> ts, const size_t interval=100);
    void displace(pele::Array<double> &coords, MC * mc) { m_ts->displace(coords, mc); }
    void report(pele::Array<double>& old_coords, const double old_energy,
            pele::Array<double>& new_coords, const double new_energy,
            const bool success, MC* mc);
};

} // namespace mcpele

#endif // #ifndef _MCPELE_ADAPTIVE_TAKESTEP_H__
