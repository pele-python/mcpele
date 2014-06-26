#ifndef _MCPELE_ACTIONS_H
#define _MCPELE_ACTIONS_H

#include <math.h>
#include <algorithm>
#include <list>
#include "pele/array.h"
#include "pele/distance.h"
#include "mc.h"
#include "histogram.h"

using std::runtime_error;
using pele::Array;
using std::sqrt;

namespace mcpele{

/*Adjust Step
 * 	factor is a multiplicative factor by which the stepsize is adjusted
 * 	niter determines the number of steps for which the action should take effect (generally
 * 	we want to adjust the step size only at the beginning of a simulation)
 * 	navg is the number of steps over which the acceptance is averaged
 * 	factor must be 0<f<1, if rejected make step shorter, if accepted make step longer
*/

class AdjustStep : public Action {
protected:
	double _target, _factor, _acceptedf;
	size_t _niter, _navg, _count, _naccepted, _nrejected;
public:
	AdjustStep(double target, double factor, size_t niter, size_t navg);
	virtual ~AdjustStep() {}
	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
};

AdjustStep::AdjustStep(double target, double factor, size_t niter, size_t navg):
			_target(target),_factor(factor),_acceptedf(0),
			_niter(niter),_navg(navg),_count(0),_naccepted(0),
			_nrejected(0){}


void AdjustStep::action(Array<double> &coords, double energy, bool accepted, MC* mc) {

	_count = mc->get_iterations_count();

	if (_count <= _niter)
		{
			if (accepted == true)
				++_naccepted;
			else
				++_nrejected;

			if(_count % _navg == 0)
			{
				_acceptedf = (double) _naccepted / (_naccepted + _nrejected);

				//std::cout<<"acceptance "<<_acceptedf<<std::endl;
				//std::cout<<"stepsize before"<<mc->_stepsize<<std::endl;
				if (_acceptedf < _target)
					mc->_stepsize *= _factor;
				else
					mc->_stepsize /= _factor;
				//std::cout<<"stepsize after"<<mc->_stepsize<<std::endl;

				//now reset to zero memory of acceptance and rejection
				_naccepted = 0;
				_nrejected = 0;
			}

		}
}

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

	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);

	Array<double> get_histogram() const {
		std::vector<double> vecdata(_hist.get_vecdata());
		Array<double> histogram(vecdata);
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

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps):
			_hist(mcpele::Histogram(min, max, bin)),
			_eqsteps(eqsteps),_count(0){}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
	_count = mc->get_iterations_count();
	if (_count > _eqsteps)
	{
	    _hist.add_entry(energy);
	}
}

/*
 * Record energy time series, measuring every __record_every-th step.
 */

class RecordEnergyTimeseries : public Action{
    private:
        void _record_energy_value(const double energy);
        const size_t _niter, _record_every;
        std::vector<double> _time_series;
    public:
        RecordEnergyTimeseries(const size_t niter, const size_t record_every);
        virtual ~RecordEnergyTimeseries(){}
        virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
        pele::Array<double> get_time_series();
        void clear(){_time_series.clear();}
};

RecordEnergyTimeseries::RecordEnergyTimeseries(const size_t niter, const size_t record_every)
    :_niter(niter),_record_every(record_every)
    {
        _time_series.reserve(niter);
        if (record_every==0) throw std::runtime_error("RecordEnergyTimeseries: __record_every expected to be at least 1");
    }

void RecordEnergyTimeseries::action(Array<double> &coords, double energy, bool accepted, MC* mc){
    size_t counter = mc->get_iterations_count();
    if (counter % _record_every == 0)
        _record_energy_value(energy);
}

void RecordEnergyTimeseries::_record_energy_value(const double energy){
    _time_series.push_back(energy);
}

pele::Array<double> RecordEnergyTimeseries::get_time_series(){
    _time_series.shrink_to_fit();
    return pele::Array<double>(_time_series).copy();
}

}
#endif
