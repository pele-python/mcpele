#ifndef _MCPELE_ACTIONS_H
#define _MCPELE_ACTIONS_H

#include <algorithm>
#include <list>

#include "pele/array.h"
#include "pele/distance.h"

#include "mc.h"
#include "histogram.h"
#include "lowest_eigenvalue.h"
#include "rsm_displacement.h"
#include "pair_dist_histogram.h"

namespace mcpele{

/*Adjust Step
 *     factor is a multiplicative factor by which the stepsize is adjusted
 *     niter determines the number of steps for which the action should take effect (generally
 *     we want to adjust the step size only at the beginning of a simulation)
 *     navg is the number of steps over which the acceptance is averaged
 *     factor must be 0<f<1, if rejected make step shorter, if accepted make step longer
*/

class AdjustStep : public Action {
protected:
    const double m_target;
    const double m_factor;
    double m_acceptedf;
    const size_t m_niter;
    const size_t m_navg;
    size_t m_count;
    size_t m_naccepted;
    size_t m_nrejected;
public:
    AdjustStep(double target, double factor, size_t niter, size_t navg);
    virtual ~AdjustStep() {}
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
};


/*
 * Record energy histogram
*/

class RecordEnergyHistogram : public Action {
protected:
    mcpele::Histogram m_hist;
private:
    const size_t m_eqsteps;
    size_t m_count;
public:
    RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps);
    virtual ~RecordEnergyHistogram(){};

    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);

    pele::Array<double> get_histogram() const {
        std::vector<double> vecdata(m_hist.get_vecdata());
        pele::Array<double> histogram(vecdata);
        return histogram.copy();
    }

    void print_terminal(size_t ntot) const {
                m_hist.print_terminal(ntot);};

    double get_max() const {
        double max_;
        max_ = m_hist.max();
        return max_;
    };

    double get_min() const {
        double min_;
        min_ = m_hist.min();
        return min_;
    };

    size_t get_eqsteps() const {return m_eqsteps;}
    double get_mean() const {return m_hist.get_mean();}
    double get_variance() const {return m_hist.get_variance();}
    int get_entries() const {return m_hist.entries();}
};

/*
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
                m_hist_gr.add_configuration(coords);
            }
        }
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

/*
 * Record scalar time series, every record_every-th step.
 */

class RecordScalarTimeseries : public Action{
private:
    const size_t m_niter;
    const size_t m_record_every;
    std::vector<double> m_time_series;
    void m_record_scalar_value(const double input)
    {
        m_time_series.push_back(input);
    }
public:
    RecordScalarTimeseries(const size_t, const size_t);
    virtual ~RecordScalarTimeseries(){}
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc) = 0;
    bool moving_average_is_stable(const size_t nr_steps_to_check = 1000, const double rel_std_threshold = 0.1);
    pele::Array<double> get_time_series()
    {
        m_time_series.shrink_to_fit();
        return pele::Array<double>(m_time_series).copy();
    }
    void clear()
    {
        m_time_series.clear();
    }
};

/*
 * Record energy time series
 */

class RecordEnergyTimeseries : public RecordScalarTimeseries{
public:
    RecordEnergyTimeseries(const size_t niter, const size_t record_every);
    virtual ~RecordEnergyTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc);
};

/*
 * Record time series of lowest eigenvalue
 */

class RecordLowestEValueTimeseries : public RecordScalarTimeseries{
private:
    FindLowestEigenvalue m_lowest_ev;
public:
    RecordLowestEValueTimeseries(const size_t niter, const size_t record_every,
            std::shared_ptr<pele::BasePotential> landscape_potential, const size_t boxdimension,
            pele::Array<double> ranvec, const size_t lbfgsniter = 30);
    virtual ~RecordLowestEValueTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy,
            const bool accepted, MC* mc);
};

/*
 * Record time series of mean displacement per particle.
 * Motivation: check if HS fluid is decorrelated between snapshots
 * This is useful if particles are not placed back into periodic box.
 */

class RecordDisplacementPerParticleTimeseries : public RecordScalarTimeseries{
private:
    GetDisplacementPerParticle m_rsm_displacement;
public:
    RecordDisplacementPerParticleTimeseries(const size_t niter, const size_t record_every,
            pele::Array<double> initial_coords, const size_t boxdimension);
    virtual ~RecordDisplacementPerParticleTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy,
                const bool accepted, MC* mc);
};

}//namespace mcpele
#endif//#ifndef _MCPELE_ACTIONS_H
