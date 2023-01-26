#ifndef _MCPELE_HISTOGRAM_H
#define _MCPELE_HISTOGRAM_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <list>

#include "pele/array.hpp"

namespace mcpele {

/*Dynamic histogram class that expand if energies outside of the initial bounds
 * are found. Being generous on the initial bounds saves a lot of time in
 * reallocation of memory, at the cost of memory preallocation. Notes:
 * ->floor always casts towards minus infinity
 * ->a list is used instead of a vector because more efficient at pushing
 * forward
 * ->begin and end return list iterators point to the beginning and the end of
 * the histogram respectively.
 * -> the most basic test that histogram must satisfy is that there must be as
 * many beads as the number of iterations (commented out at the end of the
 * script)
 * */

class Moments {
 public:
  typedef double data_t;
  typedef size_t index_t;

 private:
  data_t m_mean;
  data_t m_mean2;
  index_t m_count;

 public:
  Moments() : m_mean(0), m_mean2(0), m_count(0) {}
  void update(const data_t input) {
    m_mean = (m_mean * m_count + input) / (m_count + 1);
    m_mean2 = (m_mean2 * m_count + (input * input)) / (m_count + 1);
    if (m_count == std::numeric_limits<index_t>::max()) {
      throw std::runtime_error("Moments: update: integer overflow");
    }
    ++m_count;
  }
  /**
   * replace a data point with another one
   */
  void replace(const data_t old_data, const data_t new_data) {
    m_mean += (new_data - old_data) / m_count;
    m_mean2 += (new_data * new_data - old_data * old_data) / m_count;
  }
  void operator()(const data_t input) { update(input); }
  index_t count() const { return m_count; }
  data_t mean() const { return m_mean; }
  data_t variance() const { return (m_mean2 - m_mean * m_mean); }
  data_t std() const { return sqrt(variance()); }
};

class Histogram {
 private:
  double m_max;
  double m_min;
  double m_bin;
  double m_eps;
  int m_N;
  std::vector<double> m_hist;
  int m_niter;
  Moments m_moments;

 public:
  Histogram(const double min, const double max, const double bin);
  ~Histogram() {}
  void add_entry(double entry);
  double max() const { return m_max; }
  double min() const { return m_min; }
  double bin() const { return m_bin; }
  size_t size() const { return m_N; }
  int get_count() const { return m_niter; }
  double get_mean() const { return m_moments.mean(); }
  double get_variance() const { return m_moments.variance(); }
  std::vector<double>::iterator begin() { return m_hist.begin(); }
  std::vector<double>::iterator end() { return m_hist.end(); }
  double get_position(const size_t bin_index) const {
    return m_min + (0.5 + bin_index) * m_bin;
  }
  std::vector<double> get_vectics() const;
  std::vector<double> get_vecdata() const { return m_hist; }
  double get_entry(const size_t bin_index) const {
    return m_hist.at(bin_index);
  }
  std::vector<double> get_vecdata_error() const;
  std::vector<double> get_vecdata_normalized() const;
  void print_terminal() const;
  void resize(const double E, const int i);
};

}  // namespace mcpele

#endif  // #ifndef _MCPELE_HISTOGRAM_H
