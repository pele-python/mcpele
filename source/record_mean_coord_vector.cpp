#include <stdexcept>
#include "mcpele/record_mean_coord_vector.h"

namespace mcpele{
/*
 * m_mcv mean coordinate vector
 * m_mcv2 element wise square mean coordinate vector, useful to compute the variance in each coordinate
 * */
RecordMeanCoordVector::RecordMeanCoordVector(const size_t ndof, const size_t eqsteps)
    : m_mcv(ndof,0), //it is important that m_mcv is initialised to vec{0} for the first update step to work correctly
      m_mcv2(ndof,0),
      m_ndof(ndof),
      m_count(0),
      m_eqsteps(eqsteps)
    {}

double RecordMeanCoordVector::m_update_average(double avg, double x)
{
    double n = m_count;
    return (avg * n + x) / (n+1);
}

void RecordMeanCoordVector::m_update_mean_coord_vector(pele::Array<double> &new_coords){
    for(size_t i=0; i<m_ndof; ++i){
        double newx = new_coords[i];
        m_mcv[i] = this->m_update_average(m_mcv[i], newx);
        m_mcv2[i] = this->m_update_average(m_mcv2[i], newx*newx);
    }
    ++m_count;
}

void RecordMeanCoordVector::action(pele::Array<double> &coords, double energy, bool accepted, MC* mc)
{

    if (coords.size() != m_ndof) {
        throw std::runtime_error("RecordMeanCoordVector::action: ndof and new coords have different size");
    }
    size_t count = mc->get_iterations_count();
    if (count > m_eqsteps){
        if (accepted){
            this->m_update_mean_coord_vector(coords);
        }
        else{
            this->m_update_mean_coord_vector(m_mcv);
        }
    }
}

} //namespace mcpele
