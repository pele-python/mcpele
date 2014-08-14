#ifndef _MCPELE_MC_H
#define _MCPELE_MC_H

#include <cmath>
#include <algorithm>
#include <list>
#include <memory>
#include <stdexcept>

#include "pele/array.h"
#include "pele/base_potential.h"

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
    //Action(){std::cout<< "Action()" <<  "\n";}
    //virtual ~Action(){std::cout << "~Action()" <<  "\n";}
    virtual ~Action(){}
    virtual void action(Array<double> &coords, double energy, bool accepted,
            MC* mc) =0;
};

/*
 * Accept Test
 */

class AcceptTest{
public:
    //AcceptTest(){std::cout << "AcceptTest()" <<  "\n";}
    //virtual ~AcceptTest(){std::cout << "~AcceptTest()" <<  "\n";}
    virtual ~AcceptTest(){}
    virtual bool test(Array<double> &trial_coords, double trial_energy,
            Array<double> & old_coords, double old_energy, double temperature,
            MC * mc) =0;
};

/*
 * Accept Test
 */

class ConfTest{
public:
    //ConfTest(){std::cout << "ConfTest()" <<  "\n";}
    //virtual ~ConfTest(){std::cout << "~ConfTest()" <<  "\n";}
    virtual ~ConfTest(){}
    virtual bool conf_test(Array<double> &trial_coords, MC * mc) =0;
};

/*
 * Take Step
 */

class TakeStep{
public:
    //TakeStep(){std::cout << "TakeStep()" <<  "\n";}
    //virtual ~TakeStep(){std::cout << "TakeStep()" <<  "\n";}
    virtual ~TakeStep(){}
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
    typedef std::vector<shared_ptr<Action>> actions_t;
    typedef std::vector<shared_ptr<AcceptTest>> accept_t;
    typedef std::vector<shared_ptr<ConfTest>> conf_t;
protected:
    std::shared_ptr<pele::BasePotential> _potential;
    Array<double> _coords, _trial_coords;
    actions_t _actions;
    accept_t _accept_tests;
    conf_t _conf_tests, _late_conf_tests;
    shared_ptr<TakeStep> _takestep;
    size_t _nitercount, _accept_count, _E_reject_count, _conf_reject_count;
    bool _success;
    /*nitercount is the cumulative count, it does not get reset at the end of run*/
    bool _print_progress;
public:
    /*need to keep these public to make them accessible to tests and actions, be careful though!*/
    size_t _niter, _neval;
    double _stepsize, _temperature, _energy, _trial_energy;

    MC(std::shared_ptr<pele::BasePotential> potential, Array<double>& coords, double temperature, double stepsize);

    virtual ~MC() {}

    void one_iteration();
    void run(size_t max_iter);
    void set_temperature(double T) { _temperature = T; }
    void set_stepsize(double stepsize){ _stepsize = stepsize; }
    void add_action(shared_ptr<Action> action) { _actions.push_back(action); }
    void add_accept_test(shared_ptr<AcceptTest> accept_test)
    {
        _accept_tests.push_back(accept_test);
    }
    void add_conf_test(shared_ptr<ConfTest> conf_test)
    {
        _conf_tests.push_back(conf_test);
    }
    void add_late_conf_test(shared_ptr<ConfTest> conf_test)
    {
        _late_conf_tests.push_back(conf_test);
    }
    void set_takestep(shared_ptr<TakeStep> takestep){_takestep = takestep;}
    void set_coordinates(Array<double>& coords, double energy)
    {
        _coords = coords.copy();
        _energy = energy;
    }
    double get_energy() const { return _energy; }
    //this function is necessary if for example some potential parameter has been varied
    void reset_energy()
    {
        if(_niter > 0){
            throw std::runtime_error("MC::reset_energy after first iteration is forbidden");
        }
        _energy = compute_energy(_coords);
    }
    double get_trial_energy() const { return _trial_energy; }
    Array<double> get_coords() const { return _coords.copy(); }
    Array<double> get_trial_coords() const { return _trial_coords.copy(); }
    double get_norm_coords() const { return norm(_coords); }
    size_t get_naccept() const { return _accept_count; };
    size_t get_nreject() const { return _nitercount - _accept_count; };
    double get_accepted_fraction() const { return ((double) _accept_count)/_nitercount; };
    double get_conf_rejection_fraction() const
    {
        return ((double)_conf_reject_count)/_nitercount;
    }
    double get_E_rejection_fraction() const
    {
        return ((double)_E_reject_count)/_nitercount;
    }
    size_t get_iterations_count() const { return _nitercount; }
    size_t get_neval() const {return _neval;};
    double get_stepsize() const {return _stepsize;};
    std::shared_ptr<pele::BasePotential> get_potential_ptr() { return _potential; }
    bool take_step_specified() const { return (_takestep!=NULL); }
    void check_input();
    void set_print_progress(const bool input) { _print_progress=input; }
    void set_print_progress() { set_print_progress(true); }

protected:
    inline double compute_energy(Array<double> x)
    {
        ++_neval;
        return _potential->get_energy(x);
    }
    bool do_conf_tests(Array<double> x);
    bool do_accept_tests(Array<double> xtrial, double etrial, Array<double> xold, double eold);
    bool do_late_conf_tests(Array<double> x);
    void do_actions(Array<double> x, double energy, bool success);
};

}//namespace mcpele

#endif//#ifndef _MCPELE_MC_H
