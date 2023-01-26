#include "mcpele/record_energy_histogram.h"

using pele::Array;

namespace mcpele {

RecordEnergyHistogram::RecordEnergyHistogram(const double min, const double max,
                                             const double bin,
                                             const size_t eqsteps)
    : m_hist(min, max, bin), m_eqsteps(eqsteps), m_count(0) {}

void RecordEnergyHistogram::action(Array<double> &coords, double energy,
                                   bool accepted, MC *mc) {
  m_count = mc->get_iterations_count();
  if (m_count > m_eqsteps) {
    m_hist.add_entry(energy);
  }
}

pele::Array<double> RecordEnergyHistogram::get_histogram() const {
  std::vector<double> vecdata(m_hist.get_vecdata());
  pele::Array<double> histogram(vecdata);
  return histogram.copy();
}

}  // namespace mcpele
