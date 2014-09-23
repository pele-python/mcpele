#ifndef _MCPELE_HISTOGRAM_H
#define _MCPELE_HISTOGRAM_H

#include <cmath>
#include <algorithm>
#include <list>
#include <iostream>
#include <limits>

#include "pele/array.h"

namespace mcpele{

/*Dynamic histogram class that expand if energies outside of the initial bounds are found.
 * Being generous on the initial bounds saves a lot of time in reallocation of memory, at
 * the cost of memory preallocation.
 * Notes:
 * ->floor always casts towards minus infinity
 * ->a list is used instead of a vector because more efficient at pushing forward
 * ->begin and end return list iterators point to the beginning and the end of the
 *   histogram respectively.
 * -> the most basic test that histogram must satisfy is that there must be as many
 *      beads as the number of iterations (commented out at the end of the script)
 * */

class Moments {
public:
    typedef double data_t;
    typedef size_t index_t;
private:
    data_t _mean;
    data_t _mean2;
    index_t _count;
public:
    Moments():_mean(0), _mean2(0), _count(0) {}
    void update(const data_t input)
    {
        _mean = (_mean * _count + input) / (_count + 1);
        _mean2 = (_mean2 * _count + (input * input)) / (_count + 1);
        if (_count == std::numeric_limits<index_t>::max()) {
            throw std::runtime_error("Moments: update: integer overflow");
        }
        ++_count;
    }

    /*replace a data point with another one*/
    void replace(data_t old_data, data_t new_data)
    {
        _mean += (new_data-old_data)/_count;
        _mean2 += (new_data*new_data - old_data*old_data)/_count;
    }

    void operator() (const data_t input){ update(input); }
    index_t count() const { return _count; }
    data_t mean() const { return _mean; }
    data_t variance() const{ return (_mean2 - _mean * _mean); }
    data_t std() const { return sqrt(variance()); }
};

class Histogram{
private:
    double _max, _min, _bin, _eps;
    int _N;
    std::vector<double> _hist;
    int _niter;
    Moments moments;
public:
    Histogram(double min, double max, double bin);
    ~Histogram() {}
    void add_entry(double entry);
    double max() const {return _max;}
    double min() const {return _min;}
    double bin() const {return _bin;}
    size_t size() const {return _N;}
    int entries() const {return _niter;}
    double get_mean() const {return moments.mean();}
    double get_variance() const {return moments.variance();}
    std::vector<double>::iterator begin(){return _hist.begin();}
    std::vector<double>::iterator end(){return _hist.end();}
    double get_position(const size_t bin_index) const {return _min + (0.5 + bin_index) * _bin;}
    std::vector<double> get_vecdata() const {return _hist;}
    double get_entry(const size_t bin_index) const {return _hist.at(bin_index);}
    std::vector<double> get_vecdata_error() const;
    void print_terminal(size_t ntot) const 
    {
        for(size_t i = 0; i < _hist.size(); ++i) {
            std::cout << i << "-" << (i + 1) << ": ";
            std::cout << std::string(_hist[i] * 10000 / ntot, '*') <<  "\n";
        }
    };
    void resize(double E, int i);
};

}//namespace mcpele

#endif
