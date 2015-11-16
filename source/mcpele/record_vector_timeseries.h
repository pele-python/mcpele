#ifndef _MCPELE_RECORD_VECTOR_TIMESERIES_H__
#define _MCPELE_RECORD_VECTOR_TIMESERIES_H__

#include "mc.h"
#include <deque>
#include <cstdlib>

namespace mcpele {

/**
 * Record vector time series, every record_every-th step.
 */
class RecordVectorTimeseries : public Action {
protected:
    const size_t m_record_every, m_eqsteps;
    std::deque<pele::Array<double>> m_time_series;
    void m_record_vector_value(pele::Array<double> input)
    {
        try{
            m_time_series.push_back(input.copy());
        }
        catch(std::bad_alloc &ba){
            std::cerr<< "mcpele::RecordVectorTimeseries: bad_alloc caught: " << ba.what() << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }
public:
    RecordVectorTimeseries(const size_t record_every, const size_t eqsteps);
    virtual ~RecordVectorTimeseries(){}
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
    virtual pele::Array<double> get_recorded_vector(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc)=0;
    std::deque<pele::Array<double>> get_time_series()
    {
        m_time_series.shrink_to_fit();
        return m_time_series;
    }
    void clear() { m_time_series.clear(); }
    size_t get_record_every(){return m_record_every;}
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RECORD_VECTOR_TIMESERIES_H__
