#include <utility>
#include <cassert>
#include <stdexcept>
#include <limits>
#include <algorithm>

#include "cdf_accumulator.h"

namespace mcpele {

CDFAccumulator::CDFAccumulator()
    : m_total_number(0)
{}

void CDFAccumulator::add(const double inp)
{
    const bool already_in = m_data_x.count(inp);
    if (already_in) {
        auto it = m_data_x.find(inp);
        assert(it != m_data_x.end());
        assert(it->first == inp);
        if (it->second == std::numeric_limits<size_t>::max()) {
            throw std::range_error("CDFAccumulator::add: add would create integer overflow");
        }
        ++(it->second);
    }
    else {
        m_data_x.insert(std::make_pair(inp, size_t(1)));
    }
    if (m_total_number == std::numeric_limits<size_t>::max()) {
        throw std::range_error("CDFAccumulator::add: add would create integer overflow");
    }
    ++m_total_number;
}

std::vector<double> CDFAccumulator::get_vecdata_x() const
{
    std::vector<double> result;
    result.reserve(m_data_x.size());
    for (auto it = m_data_x.begin(); it != m_data_x.end(); ++it) {
        result.push_back(it->first);
    }
    result.swap(result);
    return result;
}

std::vector<double> CDFAccumulator::get_vecdata_cdf_x()
{
    std::vector<double> result;
    result.reserve(m_data_x.size());
    size_t remaining_x = m_total_number;
    while (!m_data_x.empty()) {
        auto data_it = m_data_x.begin();
        const double this_number = data_it->second;
        const double this_fraction = static_cast<double>(remaining_x) / static_cast<double>(m_total_number);
        remaining_x -= this_number;
        result.push_back(this_fraction);
        m_data_x.erase(data_it);
    }
    result.swap(result);
    return result;
}

} //namespace mcpele
