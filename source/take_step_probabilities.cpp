#include "mcpele/take_step_probabilities.h"

namespace mcpele {

TakeStepProbabilities::TakeStepProbabilities(const size_t seed)
    : m_generator(seed)
{}

void TakeStepProbabilities::add_step(std::shared_ptr<TakeStep> step_input, const double weight_input)
{
    m_steps.push_back(step_input);
    m_weights.push_back(weight_input);
    m_steps.swap(m_steps);
    m_weights.swap(m_weights);
    m_distribution = std::discrete_distribution<size_t>(m_weights.begin(), m_weights.end());
}

void TakeStepProbabilities::displace(pele::Array<double>& coords, MC* mc)
{
    if (m_steps.size() == 0) {
        throw std::runtime_error("TakeStepProbabilities::displace: no step specified");
    }
    m_current_index = m_distribution(m_generator);
    m_steps.at(m_current_index)->displace(coords, mc);
}

void TakeStepProbabilities::report(pele::Array<double>& old_coords, const double old_energy, pele::Array<double>& new_coords, const double new_energy, const bool success, MC* mc)
{
    m_steps.at(m_current_index)->report(old_coords, old_energy, new_coords, new_energy, success, mc);
}

} // namespace mcpele
