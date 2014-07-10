#include <cmath>

#include "pele/lbfgs.h"
#include "pele/lowest_eig_potential.h"

#include "lowest_eigenvalue.h"

namespace mcpele{


FindLowestEigenvalue::FindLowestEigenvalue(std::shared_ptr<pele::BasePotential> landscape_potential, const size_t boxdimension,
        pele::Array<double> ranvec, const size_t lbfgsniter)
    : _landscape_potential(landscape_potential),
      _boxdimension(boxdimension),
      _ranvec(ranvec.copy()),
      _lbfgsniter(lbfgsniter),
      _H0(1)
{
    if (isinf(double(1) / norm(_ranvec))) {
        throw std::runtime_error("FindLowestEigenvalue: 1/norm(_ranvec) is isinf");
    }
    _ranvec /= norm(_ranvec);
}

double FindLowestEigenvalue::get_lowest_eigenvalue(pele::Array<double> coords)
{
    std::shared_ptr<pele::LowestEigPotential> lowesteigpot = std::make_shared<pele::LowestEigPotential>(_landscape_potential, coords, _boxdimension);
    pele::LBFGS lbfgs(lowesteigpot, _ranvec.copy());
    lbfgs.set_H0(_H0);
    lbfgs.set_use_relative_f(1);
    lbfgs.run(_lbfgsniter);
    _H0 = lbfgs.get_H0();
    const double lowesteig = lbfgs.get_f();
    return lowesteig;
}


}//namespace mcpele
