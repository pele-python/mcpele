#include "take_step_pattern.h"

namespace mcpele {

void TakeStepPattern::displace(pele::Array<double>& coords, MC* mc)
{
    m_steps.get_step_index()->displace(coords, mc);
    ++m_steps;
}

} // namespace mcpele