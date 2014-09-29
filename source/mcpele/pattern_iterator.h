#ifndef _MCPELE_PATTERN_ITERATOR_H__
#define _MCPELE_PATTERN_ITERATOR_H__

#include <vector>
#include <utility>
#include <stdexcept>

#include "mc.h"

namespace mcpele {

template <class T>
class PatternIterator {
public:
    typedef typename std::vector<std::pair<size_t, T> > vec_t;
private:
    vec_t m_repetitions_indices;
    typename vec_t::const_iterator m_current;
    size_t m_current_count;
    bool m_initialized;
public:
    virtual ~PatternIterator() {}
    PatternIterator()
        : m_initialized(false)
    {}
    const PatternIterator& operator ++ ()
    {
        if (m_current_count == 0) {
            ++m_current;
            if (m_current == m_repetitions_indices.end()) {
                m_current = m_repetitions_indices.begin();
            }
            m_current_count = m_current->first;
        }
        --m_current_count;
        return *this;
    }
    void add(const T index_input, const size_t every_input=1)
    {
        if (every_input < 1) {
            throw std::range_error("PatternIterator::add: illegal input");
        }
        m_repetitions_indices.push_back(std::make_pair(every_input, index_input));
        m_current = m_repetitions_indices.begin();
        if (!m_initialized) {
            m_initialized = true;
        }
    }
    T get_step_idx() const
    {
        if (!m_initialized) {
            throw std::runtime_error("PatternIterator::get_step_idx: illegal access");
        }
        return m_current->second;
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PATTERN_ITERATOR_H__
