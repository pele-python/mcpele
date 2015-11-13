#ifndef _MCPELE_RECORD_COORDS_TIMESERIES_H__
#define _MCPELE_RECORD_COORDS_TIMESERIES_H__

#include "record_vector_timeseries.h"

namespace mcpele {

class RecordCoordsTimeseries : public RecordVectorTimeseries {
public:
    RecordCoordsTimeseries(const size_t record_every)
        : RecordVectorTimeseries(record_every)
    {}
    virtual ~RecordCoordsTimeseries() {}
    virtual pele::Array<double> get_recorded_vector(pele::Array<double> &coords,
          const double energy, const bool accepted, MC* mc) { return coords; }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RECORD_COORDS_TIMESERIES_H__
