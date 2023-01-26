#ifndef _MCPELE_RECORD_LOWEST_EVALUE_TIMESERIES_H__
#define _MCPELE_RECORD_LOWEST_EVALUE_TIMESERIES_H__

#include "lowest_eigenvalue.h"
#include "record_scalar_timeseries.h"

namespace mcpele {

/**
 * Record time series of lowest eigenvalue
 */

class RecordLowestEValueTimeseries : public RecordScalarTimeseries {
 private:
  FindLowestEigenvalue m_lowest_ev;

 public:
  RecordLowestEValueTimeseries(
      const size_t niter, const size_t record_every,
      std::shared_ptr<pele::BasePotential> landscape_potential,
      const size_t boxdimension, pele::Array<double> ranvec,
      const size_t lbfgsniter = 30)
      : RecordScalarTimeseries(niter, record_every),
        m_lowest_ev(landscape_potential, boxdimension, ranvec, lbfgsniter) {}
  virtual ~RecordLowestEValueTimeseries() {}
  virtual double get_recorded_scalar(pele::Array<double> &coords,
                                     const double energy, const bool accepted,
                                     MC *mc) {
    return m_lowest_ev.compute_lowest_eigenvalue(coords);
  }
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_RECORD_LOWEST_EVALUE_TIMESERIES_H__
