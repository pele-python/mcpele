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
    double _target, _factor, _acceptedf;
    size_t _niter, _navg, _count, _naccepted, _nrejected;
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
    mcpele::Histogram _hist;
private:
    size_t _eqsteps, _count;
public:
    RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps);
    virtual ~RecordEnergyHistogram(){};

    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);

    pele::Array<double> get_histogram() const {
        std::vector<double> vecdata(_hist.get_vecdata());
        pele::Array<double> histogram(vecdata);
        return histogram.copy();
    }

    void print_terminal(size_t ntot) const {
                _hist.print_terminal(ntot);};

    double get_max() const {
        double max_;
        max_ = _hist.max();
        return max_;
    };

    double get_min() const {
        double min_;
        min_ = _hist.min();
        return min_;
    };

    size_t get_eqsteps() const {return _eqsteps;}
    double get_mean() const {return _hist.get_mean();}
    double get_variance() const {return _hist.get_variance();}
    int get_entries() const {return _hist.entries();}
};

/*
 * Record scalar time series, every record_every-th step.
 */

class RecordScalarTimeseries : public Action{
private:
    const size_t _niter, _record_every;
    std::vector<double> _time_series;
    void _record_scalar_value(const double input)
    {
        _time_series.push_back(input);
    }
public:
    RecordScalarTimeseries(const size_t, const size_t);
    virtual ~RecordScalarTimeseries(){}
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy, const bool accepted, MC* mc) = 0;
    pele::Array<double> get_time_series()
    {
        _time_series.shrink_to_fit();
        return pele::Array<double>(_time_series).copy();
    }
    void clear()
    {
        _time_series.clear();
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
    FindLowestEigenvalue _lowest_ev;
public:
    RecordLowestEValueTimeseries(const size_t niter, const size_t record_every,
            std::shared_ptr<pele::BasePotential> landscape_potential, const size_t boxdimension,
            pele::Array<double> ranvec, const size_t lbfgsniter = 30);
    virtual ~RecordLowestEValueTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy,
            const bool accepted, MC* mc);
};

/*
 * Record time series of root mean squared displacement (averaged over all particles)
 * Motivation: check if HS fluid is decorrelated between snapshots
 * Probably this is most useful if particles are not placed back into periodic box.
 */

class RecordMeanRMSDisplacementTimeseries : public RecordScalarTimeseries{
private:
    GetMeanRMSDisplacement _rsm_displacement;
public:
    RecordMeanRMSDisplacementTimeseries(const size_t niter, const size_t record_every,
            pele::Array<double> initial_coords, const size_t boxdimension);
    virtual ~RecordMeanRMSDisplacementTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords, const double energy,
                const bool accepted, MC* mc);
};

}//namespace mcpele
#endif//#ifndef _MCPELE_ACTIONS_H
