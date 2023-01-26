#ifndef _MCPELE__RECORD_ENERGY_HISTOGRAM_H__
#define _MCPELE__RECORD_ENERGY_HISTOGRAM_H__

#include "histogram.h"
#include "mc.h"

namespace mcpele {

/*
 * Record energy histogram
 */

class RecordEnergyHistogram : public Action {
 protected:
  mcpele::Histogram m_hist;

 private:
  const size_t m_eqsteps;
  size_t m_count;

 public:
  RecordEnergyHistogram(const double min, const double max, const double bin,
                        const size_t eqsteps);
  virtual ~RecordEnergyHistogram(){};
  virtual void action(pele::Array<double> &coords, double energy, bool accepted,
                      MC *mc);
  pele::Array<double> get_histogram() const;
  void print_terminal() const { m_hist.print_terminal(); }
  double get_max() const { return m_hist.max(); }
  double get_min() const { return m_hist.min(); }
  size_t get_eqsteps() const { return m_eqsteps; }
  double get_mean() const { return m_hist.get_mean(); }
  double get_variance() const { return m_hist.get_variance(); }
  int get_count() const { return m_hist.get_count(); }
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE__RECORD_ENERGY_HISTOGRAM_H__
