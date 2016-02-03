#ifndef _MCPELE_RECORD_PAIR_DIST_HISTOGRAM_H__
#define _MCPELE_RECORD_PAIR_DIST_HISTOGRAM_H__

#include "pele/optimizer.h"

#include "mc.h"
#include "pair_dist_histogram.h"

namespace mcpele {

/**
 * Record pair-distance distribution (radial distribution function)
 * Templated on boxdimension, should work fine with pele::periodic_distance
 * Input parameters:
 * --- boxvector: defines the (periodic) simlation box
 * --- nr_bins: number of bins for g(r) histogram
 * --- eqsteps: number of equilibration steps to be excluded from g(r) computation
 * --- record_every: after more than eqsteps steps have been done, record every record_everyth step
 * Everytime the action is called, it accumulates the present configuration into the same g(r) histogram.
 * The action function calls add_configuration which accumulates the current configuration into the g(r) histogram.
 * The g(r) histogram can be read out at any point after that.
 * To read out the data, two functions are used:
 * --- get_hist_r() gives the r value array for the g(r) histogram
 * --- get_hist_gr() gives the corresponding g(r) value array, normalized using the input number of particles and number density
 *     (Admittedly number density could have been reconstructed independently of that input.)
 */

template<size_t BOXDIM>
class RecordPairDistHistogram : public Action {
private:
    mcpele::PairDistHistogram<BOXDIM> m_hist_gr;
    const size_t m_eqsteps;
    const size_t m_record_every;
public:
    RecordPairDistHistogram(pele::Array<double> boxvector, const size_t nr_bins, const size_t eqsteps, const size_t record_every)
        : m_hist_gr(boxvector, nr_bins),
          m_eqsteps(eqsteps),
          m_record_every(record_every)
    {}
    virtual ~RecordPairDistHistogram() {}
    virtual void action(pele::Array<double>& coords, double energy, bool accepted, MC* mc)
    {
        const size_t count = mc->get_iterations_count();
        if (count > m_eqsteps) {
            if (count % m_record_every == 0) {
                process_add_configuration(coords);
            }
        }
    }
    virtual void process_add_configuration(pele::Array<double>& coords)
    {
        m_hist_gr.add_configuration(coords);
    }
    size_t get_eqsteps() const
    {
        return m_eqsteps;
    }
    pele::Array<double> get_hist_r() const
    {
        std::vector<double> vecdata(m_hist_gr.get_vecdata_r());
        return pele::Array<double>(vecdata).copy();
    }
    pele::Array<double> get_hist_gr(const double number_density, const size_t nr_particles) const
    {
        std::vector<double> vecdata(m_hist_gr.get_vecdata_gr(number_density, nr_particles));
        return pele::Array<double>(vecdata).copy();
    }
};

template<size_t BOXDIM>
class RecordPairDistHistogramQuench : public RecordPairDistHistogram<BOXDIM> {
private:
    std::shared_ptr<pele::GradientOptimizer> m_opt;
public:
    RecordPairDistHistogramQuench(pele::Array<double> boxvector, const size_t nr_bins, const size_t eqsteps, const size_t record_every, std::shared_ptr<pele::GradientOptimizer> opt)
        : RecordPairDistHistogram<BOXDIM>(boxvector, nr_bins, eqsteps, record_every),
          m_opt(opt)
    {}
    void process_add_configuration(pele::Array<double>& coords)
    {
        auto c = coords.copy();
        m_opt->reset(c);
        m_opt->run();
        pele::Array<double> tmp = m_opt->get_x();
        RecordPairDistHistogram<BOXDIM>::process_add_configuration(tmp);
    }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RECORD_PAIR_DIST_HISTOGRAM_H__
