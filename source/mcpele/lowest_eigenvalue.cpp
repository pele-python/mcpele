#include <cmath>

#include "pele/lbfgs.h"
#include "pele/lowest_eig_potential.h"

#include "lowest_eigenvalue.h"

namespace mcpele{


LowestEigenvalue::LowestEigenvalue(pele::BasePotential* landscape_potential, const size_t boxdimension,
        pele::Array<double> ranvec, const double lbfgstol, const size_t lbfgsM, const size_t lbfgsniter,
        const double lbfgsmaxstep, const double H0)
    : _landscape_potential(landscape_potential),
      _boxdimension(boxdimension),
      _ranvec(ranvec.copy()),
      _lbfgstol(lbfgstol),
      _lbfgsM(lbfgsM),
      _lbfgsniter(lbfgsniter),
      _lbfgsmaxstep(lbfgsmaxstep),
      _H0(H0)
{
    if (isinf(double(1) / norm(_ranvec))) {
        throw std::runtime_error("LowestEigenvalue: 1/norm(_ranvec) is isinf");
    }
    _ranvec /= norm(_ranvec);
}

double LowestEigenvalue::get_lowest_eigenvalue(pele::Array<double> coords)
{
    pele::LowestEigPotential lowesteigpot(_landscape_potential, coords, _boxdimension);
    pele::LBFGS lbfgs(&lowesteigpot, _ranvec.copy(), _lbfgstol, _lbfgsM);
    lbfgs.set_maxstep(_lbfgsmaxstep);
    lbfgs.set_H0(_H0);
    lbfgs.set_use_relative_f(1);
    lbfgs.run(_lbfgsniter);
    _H0 = lbfgs.get_H0();
    const double lowesteig = lbfgs.get_f();
    return lowesteig;
}


}//namespace mcpele
