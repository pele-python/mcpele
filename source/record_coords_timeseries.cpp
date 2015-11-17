#include <stdexcept>
#include "mcpele/record_coords_timeseries.h"

namespace mcpele{
/*
 * m_mcv mean coordinate vector
 * m_mcv2 element wise square mean coordinate vector, useful to compute the variance in each coordinate
 * */
RecordCoordsTimeseries::RecordCoordsTimeseries(const size_t ndof, const size_t record_every, const size_t eqsteps)
    : RecordVectorTimeseries(record_every, eqsteps),
      m_mcv(ndof,0), //it is important that m_mcv is initialised to vec{0} for the first update step to work correctly
      m_mcv2(ndof,0),
      m_ndof(ndof),
      m_count(0)
    {}

double RecordCoordsTimeseries::m_update_average(double avg, double x)
{
    double n = m_count;
    return (avg * n + x) / (n+1.);
}

void RecordCoordsTimeseries::m_update_mean_coord_vector(pele::Array<double> &new_coords){
    for(size_t i=0; i<m_ndof; ++i){
        double newx = new_coords[i];
        m_mcv[i] = this->m_update_average(m_mcv[i], newx);
        m_mcv2[i] = this->m_update_average(m_mcv2[i], newx*newx);
    }
    ++m_count;
}

void RecordCoordsTimeseries::action(pele::Array<double> &coords, double energy, bool accepted, MC* mc)
{

    if (coords.size() != m_ndof) {
        throw std::runtime_error("RecordCoordsTimeseries::action: ndof and new coords have different size");
    }
    size_t counter = mc->get_iterations_count();

    if (counter > m_eqsteps){
        this->m_update_mean_coord_vector(coords);
        if (counter % m_record_every == 0) {
            m_record_vector_value(this->get_recorded_vector(coords, energy, accepted, mc));
        }
    }
}

} //namespace mcpele
