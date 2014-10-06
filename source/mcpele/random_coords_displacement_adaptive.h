#ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__
#define _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__

#include "random_coords_displacement.h"

namespace mcpele {

class RandomCoordsDisplacementAllAdaptive : public RandomCoordsDisplacementAll {
public:
    virtual ~RandomCoordsDisplacementAllAdaptive() {}
    RandomCoordsDisplacementAllAdaptive(const size_t seed,
            const double stepsize=1)
        : RandomCoordsDisplacementAll(seed, stepsize)
    {}
    void increase_acceptance(const double factor) { m_stepsize *= factor; }
    void decrease_acceptance(const double factor) { m_stepsize /= factor; }
};

class RandomCoordsDisplacementSingleAdaptive : public RandomCoordsDisplacementSingle {
public:
    virtual ~RandomCoordsDisplacementSingleAdaptive() {}
    RandomCoordsDisplacementSingleAdaptive(const size_t seed,
            const size_t nparticles, const size_t ndim, const double stepsize=1)
        : RandomCoordsDisplacementSingle(seed, nparticles, ndim, stepsize)
    {}
    void increase_acceptance(const double factor) { m_stepsize *= factor; }
    void decrease_acceptance(const double factor) { m_stepsize /= factor; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RANDOM_COORDS_DISPLACEMENT_ADAPTIVE_H__
