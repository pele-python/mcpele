#ifndef _MCPELE_MC_H
#define _MCPELE_MC_H

#include <math.h>
#include <algorithm>
#include <list>
#include <memory>
#include "pele/array.h"
#include "pele/base_potential.h"
#include <stdexcept>

using std::list;
using std::runtime_error;
using std::shared_ptr;
using pele::Array;
using std::sqrt;

namespace mcpele{

class MC;

/*
 * Action
 */

class Action {
public:
    Action(){std::cout<< "Action()" << std::endl;}
    virtual ~Action(){std::cout << "~Action()" << std::endl;}
    virtual void action(Array<double> &coords, double energy, bool accepted, MC* mc) =0;
};

/*
 * Accept Test
 */

class AcceptTest{
public:
    AcceptTest(){std::cout << "AcceptTest()" << std::endl;}
    //virtual ~AcceptTest(){}
    virtual ~AcceptTest(){std::cout << "~AcceptTest()" << std::endl;}
    virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc) =0;
};

/*
 * Accept Test
 */

class ConfTest{
public:
    ConfTest(){std::cout << "ConfTest()" << std::endl;}
    //virtual ~ConfTest(){}
    virtual ~ConfTest(){std::cout << "~ConfTest()" << std::endl;}
    virtual bool test(Array<double> &trial_coords, MC * mc) =0;
};

/*
 * Take Step
 */

class TakeStep{
public:
    TakeStep(){std::cout << "TakeStep()" << std::endl;}
    //virtual ~TakeStep(){}
    virtual ~TakeStep(){std::cout << "TakeStep()" << std::endl;}
    virtual void takestep(Array<double> &coords, double stepsize, MC * mc) =0;
};

/*
 * Monte Carlo
 * _coords and _trialcoords are arrays that store coordinates and trial coordinates respectively
 * _potential is on object of Pele::BasePotential type that defines the interaction potential
 * _E_reject_count is the count of rejections due to an energy test (e.g. Metropolis)
 * _conf_reject_count is the count of rejections due to a configuration test (e.g. spherical container)
 * _niter is the count of steps whithin a MCMC run, it is reset to zero at the end of the run
 * _nitercount is the cumulative number of MCMC steps taken by the class
 * _neval is the number of energy evaluations
 * _stepsize the is the stepsize to pass to takestep
 * _temperature is the temperature at which the simulation is performed
 * _energy is the current energy of the system
 * _success records whether the step has been accepted or rejected
 */

class MC {
public:
    typedef std::list<Action*> actions_t;
    typedef std::list<AcceptTest*> accept_t;
    typedef std::list<ConfTest*> conf_t;
protected:
    pele::BasePotential * _potential;
    Array<double> _coords, _trial_coords;
    actions_t _actions;
    accept_t _accept_tests;
    conf_t _conf_tests, _late_conf_tests;
    TakeStep* _takestep;
    size_t _nitercount, _accept_count, _E_reject_count, _conf_reject_count;
    bool _success;
    /*nitercount is the cumulative count, it does not get reset at the end of run*/
public:
    /*need to keep these public to make them accessible to tests and actions, be careful though!*/
    size_t _niter, _neval;
    double _stepsize, _temperature, _energy, _trial_energy;

    MC(pele::BasePotential * potential, Array<double>& coords, double temperature, double stepsize);

    virtual ~MC(){}

    void one_iteration();
    void run(size_t max_iter);
    void set_temperature(double T){_temperature = T;}
    void set_stepsize(double stepsize){_stepsize = stepsize;}
    void add_action(Action* action){_actions.push_back(action);}
    void add_accept_test(AcceptTest* accept_test){_accept_tests.push_back(accept_test);}
    void add_conf_test(ConfTest* conf_test){_conf_tests.push_back(conf_test);}
    void add_late_conf_test(ConfTest* conf_test){_late_conf_tests.push_back(conf_test);}
    void set_takestep( TakeStep* takestep){_takestep = takestep;}
    void set_coordinates(Array<double>& coords, double energy){
	    _coords = coords.copy();
	    _energy = energy;
    }
    double get_energy() const {return _energy;}
    //this function is necessary if for example some potential parameter has been varied
    void reset_energy(){
	if(_niter > 0){throw std::runtime_error("MC::reset_energy after first iteration is forbidden");}
	_energy = _potential->get_energy(_coords);
	++_neval;
    }
    double get_trial_energy() const {return _trial_energy;}
    Array<double> get_coords() const {return _coords.copy();}
    Array<double> get_trial_coords() const {return _trial_coords.copy();}
    double get_norm_coords() const {return norm(_coords);}
    double get_accepted_fraction() const {return ((double) _accept_count)/_nitercount;};
    double get_conf_rejection_fraction() const {return ((double)_conf_reject_count)/_nitercount;};
    double get_E_rejection_fraction() const {return ((double)_E_reject_count)/_nitercount;};
    size_t get_iterations_count() const {return _nitercount;};
    size_t get_neval() const {return _neval;};
    double get_stepsize() const {return _stepsize;};
    pele::BasePotential * get_potential_ptr(){return _potential;}
};

MC::MC(pele::BasePotential * potential, Array<double>& coords, double temperature, double stepsize):
            _potential(potential), _coords(coords.copy()),_trial_coords(_coords.copy()),
			_nitercount(0), _accept_count(0), _E_reject_count(0), _conf_reject_count(0),
			_success(true), _niter(0), _neval(0), _stepsize(stepsize), _temperature(temperature)

		{
			_energy = _potential->get_energy(_coords);
			_trial_energy = _energy;
			++_neval;
			/*std::cout<<"mcrunner Energy is "<<_energy<<std::endl;
			std::cout<<"mcrunner potential ptr is "<<_potential<<std::endl;*/
		}

void MC::one_iteration()
{
    _success = true;
    ++_niter;
    ++_nitercount;

    _trial_coords.assign(_coords);

    _takestep->takestep(_trial_coords, _stepsize, this);

    //std::cout<<"_conf_test size "<<_conf_tests.size()<<std::endl; //debug

    //for (auto& test : _conf_tests ){
    for (conf_t::iterator test = _conf_tests.begin(); test != _conf_tests.end(); ++test){
	//_success = test->test(_trial_coords, this);
	_success = (*test)->test(_trial_coords, this);
	    if (_success == false){
		    ++_conf_reject_count;
		    break;
	    }
    }

    if (_success == true)
    {
	    _trial_energy = _potential->get_energy(_trial_coords);
	    ++_neval;

	    for (accept_t::iterator test = _accept_tests.begin(); test != _accept_tests.end(); ++test){
		    _success = (*test)->test(_trial_coords, _trial_energy, _coords, _energy, _temperature, this);
		    if (_success == false){
			    ++_E_reject_count;
			    break;
		    }
	    }

	    if (_success == true){

		for (conf_t::iterator test = _late_conf_tests.begin(); test != _late_conf_tests.end(); ++test){
		    _success = (*test)->test(_trial_coords, this);
		    if (_success == false){
			++_conf_reject_count;
			break;
		    }
		}

		if (_success == true){
		    _coords.assign(_trial_coords);
		    _energy = _trial_energy;
		    ++_accept_count;
		}
	    }
    }
    for (actions_t::iterator action = _actions.begin(); action != _actions.end(); ++action){
	(*action)->action(_coords, _energy, _success, this);
    }
}

void MC::run(size_t max_iter)
{
    while(_niter < max_iter)
	    this->one_iteration();
    _niter = 0;
}

}

#endif
