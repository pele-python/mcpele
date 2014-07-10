#ifndef _MCPELE_LOWEST_EIGENVALUE_H
#define _MCPELE_LOWEST_EIGENVALUE_H

#include "pele/base_potential.h"
#include "pele/lbfgs.h"
#include "pele/lowest_eig_potential.h"

namespace mcpele{


class FindLowestEigenvalue{
private:
    std::shared_ptr<pele::LowestEigPotential> _lowesteigpot;
    size_t _lbfgsniter;
    double _H0;
    pele::Array<double> _ranvec;
    pele::LBFGS _lbfgs;
public:
    FindLowestEigenvalue(std::shared_ptr<pele::BasePotential> landscape_potential, const size_t boxdimension,
            pele::Array<double> ranvec, const size_t lbfgsniter);
    double get_lowest_eigenvalue(pele::Array<double> coords);
    double operator()(pele::Array<double> coords)
    {
        return get_lowest_eigenvalue(coords);
    }
};


}//namespace mcpele

#endif//#ifndef _MCPELE_LOWEST_EIGENVALUE_H
