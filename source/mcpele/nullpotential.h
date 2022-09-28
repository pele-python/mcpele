#ifndef _MCPELE_NULLPOTENTIAL_H
#define _MCPELE_NULLPOTENTIAL_H

#include <algorithm>
#include <functional>

#include "pele/base_potential.hpp"

namespace mcpele {

class NullPotential : public pele::BasePotential {
public:
    NullPotential() {};
    virtual ~NullPotential(){}
    virtual double inline get_energy(pele::Array<double> const & x){return 0.;};
    virtual double inline get_energy_gradient(pele::Array<double> const & x, pele::Array<double> & grad){return 0.;};
};

} // namespace mcpele

#endif // #ifndef _MCPELE_NULLPOTENTIAL_H
