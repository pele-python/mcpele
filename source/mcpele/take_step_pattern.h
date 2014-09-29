#ifndef _MCPELE_TAKE_STEP_PATTERN_H__
#define _MCPELE_TAKE_STEP_PATTERN_H__

#include "mc.h"
#include "pattern_iterator.h"

namespace mcpele {

class TakeStepPattern : public TakeStep {
private:
    std::vector<std::shared_ptr<TakeStep> > m_step_definitions;
    PatternIterator<size_t> m_step_pattern;
public:
    virtual ~TakeStepPattern() {}
    void add_step(std::shared_ptr<TakeStep> step_input, const size_t every_input);
    void displace(pele::Array<double>& coords, MC* mc);
    void report(pele::Array<double>& old_coords, const double old_energy,
            pele::Array<double>& new_coords, const double new_energy,
            const bool success, MC* mc) { get_step().report(old_coords,
                    old_energy, new_coords, new_energy, success, mc); }
    std::shared_ptr<TakeStep> get_step() const { return m_step_definitions.at(m_step_pattern.get_step_idx()); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_TAKE_STEP_PATTERN_H__
