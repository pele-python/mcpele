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
    vec_t m_step_repetitions;
    size_t m_current_step_index;
    size_t m_current_step_count;
    bool m_initialized;
public:
    virtual ~PatternManager() {}
    PatternManager()
        : m_initialized(false)
    {}
    const PatternManager& operator ++ ()
    {
        --m_current_step_count;
        if (m_current_step_count == 0) {
            ++m_current_step_index;
            if (m_current_step_index == m_step_repetitions.size()) {
                m_current_step_index = 0;
            }
            m_current_step_count = m_step_repetitions.at(m_current_step_index).first;
            assert(m_current_step_count);
        }
        return *this;
    }
    void add(const T index_input, const size_t repetitions_input=1)
    {
        if (repetitions_input < 1) {
            throw std::range_error("PatternManager::add: illegal input");
        }
        m_step_repetitions.push_back(std::make_pair(repetitions_input, index_input));
        m_current_step_index = 0;
        m_current_step_count = m_step_repetitions.at(0).first;
        if (!m_initialized) {
            m_initialized = true;
        }
    }
    T get_step_ptr() const
    {
        if (!m_initialized) {
            throw std::runtime_error("PatternManager::get_step_index: illegal access");
        }
        return m_step_repetitions.at(m_current_step_index).second;
    }
    /**
     * Return visualization of the step pattern.
     * Steps are represented by integer labels, starting from 0, in the order of
     * addition to the pattern.
     */
    std::vector<size_t> get_pattern() const
    {
        std::vector<size_t> result;
        for (typename vec_t::const_iterator i = m_step_repetitions.begin(); i != m_step_repetitions.end(); ++i) {
            const std::vector<size_t> tmp(i->first, static_cast<size_t>(i - m_step_repetitions.begin()));
            result.insert(result.end(), tmp.begin(), tmp.end());
        }
        result.swap(result);
        assert(result.front() == 0);
        assert(result.back() == m_step_repetitions.size() - 1);
        return result;
    }
    std::vector<size_t> get_pattern_direct()
    {
        m_current_step_index = 0;
        m_current_step_count = m_step_repetitions.at(0).first;
        assert(m_current_step_index == 0);
        size_t pattern_length = 0;
        for (typename vec_t::const_iterator i = m_step_repetitions.begin(); i != m_step_repetitions.end(); ++i) {
            pattern_length += i->first;
        }
        assert(m_current_step_index == 0);
        std::vector<size_t> result;
        for (size_t i = 0; i < pattern_length; ++i, ++*this) {
            result.push_back(m_current_step_index);
        }
        assert(result.size() == pattern_length);
        assert(result.front() == 0);
        assert(result.back() == m_step_repetitions.size() - 1);
        return result;
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_PATTERN_MANAGER_H__
