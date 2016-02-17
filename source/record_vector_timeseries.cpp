#include "mcpele/record_vector_timeseries.h"

using pele::Array;

namespace mcpele {

RecordVectorTimeseries::RecordVectorTimeseries(const size_t record_every, const size_t eqsteps)
    : m_record_every(record_every),
      m_eqsteps(eqsteps)
{
    if (record_every == 0) {
        throw std::runtime_error("RecordVectorTimeseries: record_every expected to be at least 1");
    }
}

void RecordVectorTimeseries::action(Array<double> &coords, double energy, bool accepted, MC* mc)
{
    const size_t counter = mc->get_iterations_count();
    if (counter % m_record_every == 0 && counter > m_eqsteps) {
        m_record_vector_value(this->get_recorded_vector(coords, energy, accepted, mc));
    }
}

} // namespace mcpele
