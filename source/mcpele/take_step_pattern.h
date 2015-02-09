#ifndef _MCPELE_TAKE_STEP_PATTERN_H__
#define _MCPELE_TAKE_STEP_PATTERN_H__

#include "pattern_manager.h"

namespace mcpele {

/**
 * Create a pattern of TakeSteps to be exectuted by MC.
 *
 * Example
 * -------
 *
 *      To have MC looping over a pattern consisting of 99 times step_A and
 *      1 time step_B, do e.g.:
 *          auto pot = std::make_shared<SomePotential>(SomeParameters);
 *          auto mc = std::make_shared<mcpele::MC>(pot, SomeCoordinates, SomeTemperature);
 *          auto step_pattern = std::make_shared<mcpele::TakeStepPattern>();
 *          auto step_A = std::make_shared<SomeTakeStep>(SomeParameters);
 *          auto step_B = std::make_shared<SomeOtherTakeStep>(SomeOtherParameters);
 *          step_pattern->add_step(step_A, 99);
 *          step_pattern->add_step(step_B, 1);
 *          mc->set_takestep(step_pattern);
 *          mc->run(1e6);
 */
class TakeStepPattern : public TakeStep {
private:
    PatternManager<size_t> m_steps;
    std::vector<std::shared_ptr<TakeStep> > m_step_storage;
public:
    virtual ~TakeStepPattern() {}
    void add_step(std::shared_ptr<TakeStep> step_input,
            const size_t repetitions_input=1)
    {
        m_steps.add(m_step_storage.size(), repetitions_input);
        m_step_storage.push_back(step_input);
    }
    void displace(pele::Array<double>& coords, MC* mc);
    void report(pele::Array<double>& old_coords, const double old_energy,
                pele::Array<double>& new_coords, const double new_energy,
                const bool success, MC* mc)
    {
        m_step_storage.at(m_steps.get_step_ptr())->report(old_coords,
                        old_energy, new_coords, new_energy, success, mc);
    }
    std::vector<size_t> get_pattern() const { return m_steps.get_pattern(); }
    std::vector<size_t> get_pattern_direct() { return m_steps.get_pattern_direct(); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_TAKE_STEP_PATTERN_H__
