#include "mcpele/take_step_pattern.h"

namespace mcpele {

void TakeStepPattern::displace(pele::Array<double>& coords, MC* mc)
{
    ++m_steps;
    m_step_storage.at(m_steps.get_step_ptr())->displace(coords, mc);
}

} // namespace mcpele
