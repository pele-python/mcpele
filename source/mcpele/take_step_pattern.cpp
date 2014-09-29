#include "take_step_pattern.h"

namespace mcpele {

void TakeStepPattern::displace(pele::Array<double>& coords, MC* mc)
{
    get_step().displace(coords, mc);
    ++m_step_pattern;
}

} // namespace mcpele
