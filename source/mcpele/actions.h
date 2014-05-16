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
	double _bin, _mean, _mean2; //first and second moment of the distribution
	size_t _eqsteps, _count;
public:
	RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps);
	virtual ~RecordEnergyHistogram(){};

	virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);

	virtual Array<double> get_histogram(){
		std::vector<double> vecdata =_hist.get_vecdata();
		Array<double> histogram(vecdata);
		Array<double> histogram2(histogram.copy());
		return histogram2;
	}

	virtual void print_terminal(size_t ntot){
				_hist.print_terminal(ntot);};

	virtual double get_max(){
		double max_;
		max_ = _hist.max();
		return max_;
	};

	virtual double get_min(){
			double min_;
			min_ = _hist.min();
			return min_;
		};

	virtual double get_mean(){return _mean;};
	virtual double get_variance(){return (_mean2 - _mean*_mean);};
};

RecordEnergyHistogram::RecordEnergyHistogram(double min, double max, double bin, size_t eqsteps):
			_hist(mcpele::Histogram(min, max, bin)),_bin(bin),_mean(0.),_mean2(0.),
			_eqsteps(eqsteps),_count(0){}

void RecordEnergyHistogram::action(Array<double> &coords, double energy, bool accepted, MC* mc) {
	_count = mc->get_iterations_count();
	if (_count > _eqsteps){
		_hist.add_entry(energy);
	    double count = _count - _eqsteps + 1;
		_mean = (_mean*(count-1)+energy)/count;
	    _mean2 = (_mean2*(count-1)+(energy*energy))/count;
	}
}

/*
 * Record energy time series, measuring every __record_every-th step.
 */

class RecordEnergyTimeseries : public Action{
    private:
        void _record_energy_value(const double energy);
        const size_t _record_every;
        size_t _counter;
        std::vector<double> _time_series;
    public:
        RecordEnergyTimeseries(const size_t record_every);
        virtual ~RecordEnergyTimeseries(){}
        virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc);
        pele::Array<double> get_time_series();
};

RecordEnergyTimeseries::RecordEnergyTimeseries(const size_t record_every)
    :_record_every(record_every),_counter(0)
    {
        if (record_every==0) throw std::runtime_error("RecordEnergyTimeseries: __record_every expected to be at least 1");
    }

void RecordEnergyTimeseries::action(Array<double> &coords, double energy, bool accepted, MC* mc){
    ++_counter;
    if (_counter % _record_every == 0)
        _record_energy_value(energy);
}

void RecordEnergyTimeseries::_record_energy_value(const double energy){
    _time_series.push_back(energy);
}

pele::Array<double> RecordEnergyTimeseries::get_time_series(){
    _time_series.shrink_to_fit();
    return pele::Array<double>(_time_series);
}

}
#endif
