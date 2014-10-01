#include <cmath>

#include "mcpele/lowest_eigenvalue.h"

namespace mcpele{

FindLowestEigenvalue::FindLowestEigenvalue(std::shared_ptr<pele::BasePotential> landscape_potential, const size_t boxdimension,
        const pele::Array<double> ranvec, const size_t lbfgsniter)
    : m_lowesteigpot(std::make_shared<pele::LowestEigPotential>(landscape_potential, ranvec.copy(), boxdimension)),
      m_ranvec((ranvec.copy() /= norm(ranvec))),
      m_lbfgs(m_lowesteigpot, m_ranvec.copy())
{
    if (isinf(double(1) / norm(ranvec))) {
        throw std::runtime_error("FindLowestEigenvalue: 1/norm(ranvec) is isinf");
    }
    m_lbfgs.set_max_iter(lbfgsniter);
}

double FindLowestEigenvalue::compute_lowest_eigenvalue(pele::Array<double> coords)
{
    m_lowesteigpot->reset_coords(coords);
    m_lbfgs.reset(m_ranvec);
    m_lbfgs.set_use_relative_f(1);
    m_lbfgs.run();
    const double lowesteig = m_lbfgs.get_f();
    return lowesteig;
}

}//namespace mcpele
