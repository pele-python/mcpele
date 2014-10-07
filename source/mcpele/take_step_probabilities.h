#ifndef _MCPELE_TAKE_STEP_PROBABILITIES_H__
#define _MCPELE_TAKE_STEP_PROBABILITIES_H__

#include <vector>
#include <random>

#include "mc.h"

namespace mcpele {

/**
 * Create a step pattern similar to TakeStepPattern.
 * However, the steps are specified together with their relative weights and
 * exectured accoringly.
 *
 * Reference
 * ---------
 * http://www.cplusplus.com/reference/random/discrete_distribution/
 */
class TakeStepProbabilities : public TakeStep {
private:
    std::vector<std::shared_ptr<TakeStep> > m_steps;
    std::vector<double> m_weights;
    std::discrete_distribution<size_t> m_distribution;
    std::mt19937_64 m_generator;
    size_t m_current_index;
public:
    virtual ~TakeStepProbabilities() {}
    TakeStepProbabilities(const size_t seed);
    void add_step(std::shared_ptr<TakeStep> step_input, const double weight_input=1);
    void displace(pele::Array<double>& coords, MC* mc);
    void report(pele::Array<double>& old_coords, const double old_energy,
            pele::Array<double>& new_coords, const double new_energy,
            const bool success, MC* mc);
    std::vector<double> get_weights() const { return m_weights; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_TAKE_STEP_PROBABILITIES_H__
