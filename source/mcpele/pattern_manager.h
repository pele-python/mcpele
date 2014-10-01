#ifndef _MCPELE_PATTERN_MANAGER_H__
#define _MCPELE_PATTERN_MANAGER_H__

#include <vector>
#include <utility>
#include <stdexcept>

#include "mc.h"

namespace mcpele {

template <class T>
class PatternManager {
public:
    typedef typename std::vector<std::pair<size_t, T> > vec_t;
private:
    vec_t m_repetitions_indices;
    typename vec_t::const_iterator m_current;
    size_t m_current_count;
    bool m_initialized;
public:
    virtual ~PatternManager() {}
    PatternManager()
        : m_initialized(false)
    {}
    const PatternManager& operator ++ ()
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
    void add(const T index_input, const size_t repetitions_input=1)
    {
        if (repetitions_input < 1) {
            throw std::range_error("PatternManager::add: illegal input");
        }
        m_repetitions_indices.push_back(std::make_pair(repetitions_input, index_input));
        m_current = m_repetitions_indices.begin();
        if (!m_initialized) {
            m_initialized = true;
        }
    }
    T get_step_index() const
    {
        if (!m_initialized) {
            throw std::runtime_error("PatternManager::get_step_index: illegal access");
        }
        return m_current->second;
    }
    /**
     * Return visualization of the step pattern.
     * Steps are repeseted by integer labels, starting from 0, in the order of
     * addition to the pattern.
     */
    std::vector<size_t> get_pattern() const
    {
        std::vector<size_t> result;
        for (typename vec_t::const_iterator i = m_repetitions_indices.begin(); i != m_repetitions_indices.end(); ++i) {
            const std::vector<size_t> tmp(i->first, static_cast<size_t>(i - m_repetitions_indices.begin()));
            result.insert(result.end(), tmp.begin(), tmp.end());
        }
        result.swap(result);
        return result;
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PATTERN_MANAGER_H__