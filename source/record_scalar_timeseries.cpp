#include "mcpele/record_scalar_timeseries.h"

#include "mcpele/moving_average.h"

using pele::Array;

namespace mcpele {

RecordScalarTimeseries::RecordScalarTimeseries(const size_t niter,
                                               const size_t record_every)
    : m_record_every(record_every) {
  if (record_every == 0) {
    throw std::runtime_error(
        "RecordScalarTimeseries: record_every expected to be at least 1");
  }
  m_time_series.reserve(niter / record_every);
}

void RecordScalarTimeseries::action(Array<double> &coords, double energy,
                                    bool accepted, MC *mc) {
  const size_t counter = mc->get_iterations_count();
  if (counter % m_record_every == 0) {
    m_record_scalar_value(
        this->get_recorded_scalar(coords, energy, accepted, mc));
  }
}

}  // namespace mcpele
