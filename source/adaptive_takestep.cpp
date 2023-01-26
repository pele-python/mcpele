#include "mcpele/adaptive_takestep.h"

namespace mcpele {

AdaptiveTakeStep::AdaptiveTakeStep(std::shared_ptr<TakeStep> ts,
                                   const size_t interval, const double factor,
                                   const double min_acceptance_ratio,
                                   const double max_acceptance_ratio)
    : m_ts(ts),
      m_interval(interval),
      m_total_steps(0),
      m_accepted_steps(0),
      m_factor(factor),
      m_min_acceptance_ratio(min_acceptance_ratio),
      m_max_acceptance_ratio(max_acceptance_ratio) {
  if (factor <= 0 || factor >= 1) {
    throw std::runtime_error(
        "AdaptiveTakeStep::AdaptiveTakeStep: should be between 0 and 1");
  }
}

void AdaptiveTakeStep::report(pele::Array<double> &, const double,
                              pele::Array<double> &, const double,
                              const bool success, MC *mc) {
  ++m_total_steps;
  if (success) {
    ++m_accepted_steps;
  }
  if (mc->get_iterations_count() % m_interval == 0) {
    const double acceptance_fraction = static_cast<double>(m_accepted_steps) /
                                       static_cast<double>(m_total_steps);
    m_accepted_steps = 0;
    m_total_steps = 0;
    if (acceptance_fraction < get_min_acceptance_ratio()) {
      m_ts->increase_acceptance(m_factor);
    } else if (acceptance_fraction > get_max_acceptance_ratio()) {
      m_ts->decrease_acceptance(m_factor);
    }
  }
}

}  // namespace mcpele
