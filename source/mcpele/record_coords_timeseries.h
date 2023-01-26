#ifndef _MCPELE_RECORD_COORDS_TIMESERIES_H__
#define _MCPELE_RECORD_COORDS_TIMESERIES_H__

#include "record_vector_timeseries.h"

namespace mcpele {

class RecordCoordsTimeseries : public RecordVectorTimeseries {
 private:
  pele::Array<double> m_mcv, m_mcv2;
  const size_t m_ndof;
  size_t m_count;
  double m_update_average(double avg, double x);
  void m_update_mean_coord_vector(pele::Array<double> &new_coords);

 public:
  RecordCoordsTimeseries(const size_t ndof, const size_t record_every,
                         const size_t eqsteps);
  virtual ~RecordCoordsTimeseries() {}
  virtual pele::Array<double> get_recorded_vector(pele::Array<double> &coords,
                                                  const double energy,
                                                  const bool accepted, MC *mc) {
    return coords;
  }
  virtual void action(pele::Array<double> &coords, double energy, bool accepted,
                      MC *mc);
  pele::Array<double> get_mean_coordinate_vector() { return m_mcv.copy(); }
  pele::Array<double> get_mean2_coordinate_vector() { return m_mcv2.copy(); }
  pele::Array<double> get_variance_coordinate_vector() {
    pele::Array<double> var = m_mcv2.copy();
    for (size_t i = 0; i < m_ndof; ++i) {
      var[i] -= m_mcv[i] * m_mcv[i];
    }
    return var.copy();
  }
  size_t get_count() { return m_count; }
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_RECORD_COORDS_TIMESERIES_H__
