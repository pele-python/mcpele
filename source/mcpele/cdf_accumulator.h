#ifndef _MCPELE_CDF_ACCUMULATOR_H
#define _MCPELE_CDF_ACCUMULATOR_H

#include <map>
#include <vector>

namespace mcpele{

class CDFAccumulator {
private:
    std::map<double, size_t> m_data_x;
    size_t m_total_number;
public:
    virtual ~CDFAccumulator() {}
    CDFAccumulator();
    void add(const double);
    std::vector<double> get_vecdata_x() const;
    std::vector<double> get_vecdata_cdf_x();
};

} //namespace mcpele

#endif //#ifndef _MCPELE_CDF_ACCUMULATOR_H
