#include "take_step_pattern.h"

namespace mcpele {

void TakeStepPattern::add_step(std::shared_ptr<TakeStep> step_input, const size_t every_input)
{
    m_step_definitions.push_back(step_input);
    m_step_pattern.add(every_input);
}

void TakeStepPattern::displace(pele::Array<double>& coords, MC* mc)
{
    get_step().displace(coords, mc);
    ++m_step_pattern;
}

} // namespace mcpele
