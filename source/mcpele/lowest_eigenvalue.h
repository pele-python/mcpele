#ifndef _MCPELE_LOWEST_EIGENVALUE_H
#define _MCPELE_LOWEST_EIGENVALUE_H

#include "pele/base_potential.h"

namespace mcpele{


class FindLowestEigenvalue{
private:
    pele::BasePotential* _landscape_potential;
    const size_t _boxdimension;
    pele::Array<double> _ranvec;
    const double _lbfgstol;
    const size_t _lbfgsM;
    const size_t _lbfgsniter;
    const double _lbfgsmaxstep;
    double _H0;
public:
    FindLowestEigenvalue(pele::BasePotential* landscape_potential, const size_t boxdimension,
            pele::Array<double> ranvec, const double lbfgstol, const size_t lbfgsM,
            const size_t lbfgsniter, const double lbfgsmaxstep, const double H0);
    double get_lowest_eigenvalue(pele::Array<double> coords);
    double operator()(pele::Array<double> coords)
    {
        return get_lowest_eigenvalue(coords);
    }
};


}//namespace mcpele

#endif//#ifndef _MCPELE_LOWEST_EIGENVALUE_H
