#include "histogram.h"

namespace mcpele{

Histogram::Histogram(double min, double max, double bin)
    : _max(floor((max/bin)+1)*bin),
    _min(floor((min/bin))*bin),
    _bin(bin),
    _eps(std::numeric_limits<double>::epsilon()),
    _N((_max - _min) / bin),
    _hist(_N,0),_niter(0)
{
#ifdef DEBUG
    std::cout<<"histogram is of size "<<_N<< "\n";
#endif
}

void Histogram::add_entry(double E)
{
    moments(E);
    int i;
    E = E + _eps; //this is a dirty hack, not entirely sure of its generality and possible consequences, tests seem to be fine
    i = floor((E-_min)/_bin);
    if (i < _N && i >= 0) {
        _hist[i] += 1;
        ++_niter;
    }
    else
        this->resize(E,i);

    /*THIS IS A TEST*/
    /*int renorm = 0;
     * for(vector<size_t>::iterator it = _hist.begin();it != _hist.end();++it)
      {
          renorm += *it;
      }

    if (renorm != _niter)
    {
        std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n renorm "<<renorm<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<< "\n";
        assert(renorm == _niter);
    }*/
}

void Histogram::resize(double E, int i)
{
    int newlen;
    if (i >= _N) {
        newlen = (i + 1) - _N;
        _hist.insert(_hist.end(), (newlen-1), 0);
        _hist.push_back(1);
        ++_niter;
        _max = floor((E/_bin)+1)*_bin; //round to nearest increment
        _N = round((_max - _min) / _bin); //was round
        if ( (int) _hist.size() != _N) {
            std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n size "<<_hist.size()<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<< "\n";
            assert( (int) _hist.size() == _N);
            exit (EXIT_FAILURE);
        }
        std::cout<<"resized above at niter "<<_niter<< "\n";
    } else if (i < 0) {
        newlen = -1*i;
        _hist.insert(_hist.begin(), (newlen-1), 0);
        _hist.insert(_hist.begin(),1);
        ++_niter;
        _min = floor((E/_bin))*_bin; //round to nearest increment
        _N = round((_max-_min)/_bin); //was round
        if ( (int) _hist.size() != _N) {
            std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n size "<<_hist.size()<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<< "\n";
            assert( (int) _hist.size() == _N);
            exit (EXIT_FAILURE);
        }
        std::cout<<"resized below at niter "<<_niter<< "\n";
    } else {
        std::cerr<<"histogram encountered unexpected condition"<< "\n";
        std::cout<<" E "<<E<<"\n niter "<<_niter<<"\n min "<<_min<<"\n max "<<_max<<"\n i "<<i<<"\n N "<<_N<< "\n";
    }
}

/*
 * Note: This gives the error bar on a bin of width _bin, under the assumption that the sum of all bin areas is 1.
 * */
std::vector<double> Histogram::get_vecdata_error() const
{
    std::vector<double> result(_hist.size(), 0);
    for (size_t i = 0; i < result.size(); ++i) {
        const double this_fraction = static_cast<double>(_hist.at(i)) / static_cast<double>(entries());
        result.at(i) = sqrt(this_fraction * (1 - this_fraction) / _bin) / sqrt(entries());
    }
    return result;
}


}//namespace mcpele
