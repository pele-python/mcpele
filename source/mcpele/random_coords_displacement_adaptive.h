#ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__
#define _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__

#include "random_coords_displacement.h"

namespace mcpele {

class RandomCoordsDisplacementAdaptive : public RandomCoordsDisplacement {
public:
    virtual ~RandomCoordsDisplacementAdaptive() {}
    RandomCoordsDisplacementAdaptive(const size_t seed,
            const double stepsize=1)
        : RandomCoordsDisplacement(seed, stepsize)
    {}
    void increase_acceptance(const double factor) { m_stepsize *= factor; }
    void decrease_acceptance(const double factor) { m_stepsize /= factor; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__
