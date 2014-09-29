#ifndef _MCPELE_RECORD_ENERGY_TIMESERIES_H__
#define _MCPELE_RECORD_ENERGY_TIMESERIES_H__

#include "record_scalar_timeseries.h"

namespace mcpele {

class RecordEnergyTimeseries : public RecordScalarTimeseries{
public:
    RecordEnergyTimeseries(const size_t niter, const size_t record_every)
        : RecordScalarTimeseries(niter, record_every)
    {}
    virtual ~RecordEnergyTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords,
            const double energy, const bool accepted, MC* mc) { return energy; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RECORD_ENERGY_TIMESERIES_H__
