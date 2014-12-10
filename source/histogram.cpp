#include "mcpele/histogram.h"

namespace mcpele{

Histogram::Histogram(const double min, const double max, const double bin)
    : m_max(floor((max / bin) + 1) * bin),
      m_min(floor((min / bin)) * bin),
      m_bin(bin),
      m_eps(std::numeric_limits<double>::epsilon()),
      m_N((m_max - m_min) / bin),
      m_hist(m_N, 0),
      m_niter(0)
{
#ifdef DEBUG
    std::cout << "histogram is of size " << m_N << "\n";
#endif
}

void Histogram::add_entry(double E)
{
    m_moments(E);
    int i;
    E = E + m_eps; //this is a dirty hack, not entirely sure of its generality and possible consequences, tests seem to be fine
    i = floor((E - m_min) / m_bin);
    if (i < m_N && i >= 0) {
        m_hist[i] += 1;
        ++m_niter;
    }
    else {
        this->resize(E, i);
    }

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

void Histogram::resize(const double E, const int i)
{
    int newlen;
    if (i >= m_N) {
        newlen = (i + 1) - m_N;
        m_hist.insert(m_hist.end(), (newlen - 1), 0);
        m_hist.push_back(1);
        ++m_niter;
        m_max = floor((E / m_bin) + 1) * m_bin; //round to nearest increment
        m_N = round((m_max - m_min) / m_bin); //was round
        if (static_cast<int>(m_hist.size()) != m_N) {
            std::cout<< " E " << E << "\n niter " << m_niter<< "\n size " << m_hist.size() << "\n min " << m_min << "\n max " << m_max << "\n i " << i << "\n N " << m_N << "\n";
            assert(static_cast<int>(m_hist.size()) == m_N);
            exit (EXIT_FAILURE);
        }
        std::cout<< "resized above at niter " << m_niter << "\n";
    }
    else if (i < 0) {
        newlen = -1 * i;
        m_hist.insert(m_hist.begin(), (newlen - 1), 0);
        m_hist.insert(m_hist.begin(),1);
        ++m_niter;
        m_min = floor((E / m_bin)) * m_bin; //round to nearest increment
        m_N = round((m_max - m_min) / m_bin); //was round
        if (static_cast<int>(m_hist.size()) != m_N) {
            std::cout<<" E "<< E << "\n niter " << m_niter << "\n size " << m_hist.size() << "\n min " << m_min << "\n max " << m_max << "\n i " << i << "\n N " << m_N << "\n";
            assert(static_cast<int>(m_hist.size()) == m_N);
            exit (EXIT_FAILURE);
        }
        std::cout<< "resized below at niter " << m_niter << "\n";
    }
    else {
        std::cerr << "histogram encountered unexpected condition\n";
        std::cout << " E " << E << "\n niter " << m_niter << "\n min " << m_min << "\n max " << m_max << "\n i " << i << "\n N " << m_N << "\n";
    }
}

/*
 * Note: This gives the error bar on a bin of width _bin, under the assumption that the sum of all bin areas is 1.
 * */
std::vector<double> Histogram::get_vecdata_error() const
{
    std::vector<double> result(m_hist.size(), 0);
    for (size_t i = 0; i < result.size(); ++i) {
        const double this_fraction = static_cast<double>(m_hist.at(i)) / static_cast<double>(entries());
        result.at(i) = sqrt(this_fraction * (1 - this_fraction) / m_bin) / sqrt(entries());
    }
    return result;
}

std::vector<double> Histogram::get_vecdata_normalized() const
{
    std::vector<double> result(m_hist.size(), 0);
    for (size_t i = 0; i < m_hist.size(); ++i) {
        const double this_fraction = static_cast<double>(m_hist.at(i)) / static_cast<double>(m_niter);
        result.at(i) = this_fraction / m_bin;
    }
    return result;
}

void Histogram::print_terminal() const
{
    for(size_t i = 0; i < m_hist.size(); ++i) {
        std::cout << i << "-" << (i + 1) << ": ";
        std::cout << std::string(m_hist[i] * 10000 / m_niter, '*') <<  "\n";
    }
}


}//namespace mcpele
