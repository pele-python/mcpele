#ifndef _MCPELE_RECORD_SCALAR_TIMESERIES_H__
#define _MCPELE_RECORD_SCALAR_TIMESERIES_H__

#include "mc.h"

namespace mcpele {

/**
 * Record scalar time series, every record_every-th step.
 */
class RecordScalarTimeseries : public Action {
private:
    const size_t m_record_every;
    std::vector<double> m_time_series;
    void m_record_scalar_value(const double input)
    {
        m_time_series.push_back(input);
    }
public:
    RecordScalarTimeseries(const size_t, const size_t);
    virtual ~RecordScalarTimeseries(){}
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)=0;
    pele::Array<double> get_time_series()
    {
        m_time_series.shrink_to_fit();
        return pele::Array<double>(m_time_series).copy();
    }
    void clear() { m_time_series.clear(); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RECORD_SCALAR_TIMESERIES_H__
