#ifndef _MCPELE_TAKE_STEP_PATTERN_H__
#define _MCPELE_TAKE_STEP_PATTERN_H__

#include "pattern_iterator.h"

namespace mcpele {

class TakeStepPattern : public TakeStep {
private:
    PatternIterator<std::shared_ptr<TakeStep> > m_steps;
public:
    virtual ~TakeStepPattern() {}
    void add_step(std::shared_ptr<TakeStep> step_input,
            const size_t every_input=1) { m_steps.add(step_input, every_input); }
    void displace(pele::Array<double>& coords, MC* mc);
    void report(pele::Array<double>& old_coords, const double old_energy,
            pele::Array<double>& new_coords, const double new_energy,
            const bool success, MC* mc) { m_steps.get_step_idx()->report(
                    old_coords, old_energy, new_coords, new_energy, success, mc); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_TAKE_STEP_PATTERN_H__
