#ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__
#define _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__

#include "random_coords_displacement.h"

namespace mcpele {

class RandomCoordsDisplacementAdaptive : public RandomCoordsDisplacement {
protected:
    const double m_factor;
    const double m_min_acc_frac;
    const double m_max_acc_frac;
public:
    virtual ~RandomCoordsDisplacementAdaptive() {}
    RandomCoordsDisplacementAdaptive(const size_t seed,
            const double stepsize=1, const double factor=0.9,
            const double min_acc_frac=0.2, const double max_acc_frac=0.5)
        : RandomCoordsDisplacement(seed, stepsize),
          m_factor(factor),
          m_min_acc_frac(min_acc_frac),
          m_max_acc_frac(max_acc_frac)
    {}
    void increase_acceptance() { m_stepsize *= m_factor; }
    void decrease_acceptance() { m_stepsize /= m_factor; }
    double get_min_acceptance_ratio() const { return m_min_acc_frac; }
    double get_max_acceptance_ratio() const { return m_max_acc_frac; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__
