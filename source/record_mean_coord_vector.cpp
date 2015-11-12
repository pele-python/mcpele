#include "mcpele/record_mean_coord_vector.h"

namespace mcpele{

RecordMeanCoordVector::RecordMeanCoordVector(pele::Array<double> initial_coordinates, const size_t eqsteps)
    : m_mcv(initial_coordinates.size(),0),
      m_ndof(initial_coordinates.size()),
      m_count(0),
      m_eqsteps(eqsteps)
    {}

double RecordMeanCoordVector::_update_average(double avg, double x)
{
    double n = m_count;
    return (avg * n + x) / (n+1);
}

void RecordMeanCoordVector::_update_mean_coord_vector(pele::Array<double> new_coords){
    for(size_t i=0; i<m_ndof; ++i){
        m_mcv[i] = this->_updade_average(m_mcv[i], new_coords[i]);
    }
    ++m_count;
}

void RecordMeanCoordVector::action(Array<double> &coords, double energy, bool accepted, MC* mc)
{
    size_t count = mc->get_iterations_count();
    if (count > m_eqsteps){
        if (accepted){
            this->_update_mean_coord_vector(coords);
        }
        else{
            this->_update_mean_coord_vector(m_mcv);
        }
    }
}

} //namespace mcpele
