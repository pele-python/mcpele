#ifndef _MCPELE_MC_H
#define _MCPELE_MC_H

#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>

#include "pele/array.h"
#include "pele/base_potential.h"

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
    virtual void action(pele::Array<double> &coords, double energy, bool accepted,
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
    virtual bool test(pele::Array<double> &trial_coords, double trial_energy,
            pele::Array<double> & old_coords, double old_energy, double temperature,
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
    virtual bool conf_test(pele::Array<double> &trial_coords, MC * mc) =0;
};

/*
 * Take Step
 */

class TakeStep {
public:
    virtual ~TakeStep() {}
    virtual void displace(pele::Array<double>& coords, MC* mc) = 0;
    virtual void report(const bool success) {}
    virtual void increase_acceptance() {}
    virtual void decrease_acceptance() {}
};

class AdaptiveTakeStep : public TakeStep {
protected:
    std::shared_ptr<TakeStep> m_ts;
    size_t m_interval;
    size_t m_total_steps;
    size_t m_accepted_steps;
public:
    virtual ~AdaptiveTakeStep() {}
    AdaptiveTakeStep(std::shared_ptr<TakeStep> ts, const size_t interval=100)
        : m_ts(ts),
          m_interval(interval),
          m_total_steps(0),
          m_accepted_steps(0)
    {}
    void displace(pele::Array<double> &coords, MC * mc)
    {
        m_ts->displace(coords, mc);
    }
    void report(const bool success)
    {
        ++m_total_steps;
        if (success) {
            ++m_accepted_steps;
        }
        if (mc->get_iterations_count() % m_interval == 0) {
            const double acceptance_fraction = static_cast<double>(m_accepted_steps) / static_cast<double>(m_total_steps);
            m_accepted_steps = 0;
            m_total_steps = 0;
            if (acceptance_fraction < ts->get_min_acc_fraction()) {
                ts->increase_acceptance();
            }
            else if (acceptance_fraction > ts->get_max_acc_fraction()) {
                ts->decrease_acceptance();
            }
        }
    }
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
    typedef std::vector<std::shared_ptr<Action> > actions_t;
    typedef std::vector<std::shared_ptr<AcceptTest> > accept_t;
    typedef std::vector<std::shared_ptr<ConfTest> > conf_t;
protected:
    std::shared_ptr<pele::BasePotential> _potential;
    pele::Array<double> _coords, _trial_coords;
    actions_t _actions;
    accept_t _accept_tests;
    conf_t _conf_tests;
    conf_t _late_conf_tests;
    std::shared_ptr<TakeStep> _take_step;
    size_t _nitercount, _accept_count, _E_reject_count, _conf_reject_count;
    bool _success;
    /*nitercount is the cumulative count, it does not get reset at the end of run*/
    bool _print_progress;
public:
    /*need to keep these public to make them accessible to tests and actions, be careful though!*/
    size_t _niter, _neval;
    double _stepsize, _temperature, _energy, _trial_energy;

    MC(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>& coords, double temperature, double stepsize);

    virtual ~MC() {}

    void one_iteration();
    void run(size_t max_iter);
    void set_temperature(double T) { _temperature = T; }
    double get_temperature(){return _temperature;}
    void set_stepsize(double stepsize){ _stepsize = stepsize; }
    void add_action(std::shared_ptr<Action> action) { _actions.push_back(action); }
    void add_accept_test(std::shared_ptr<AcceptTest> accept_test)
    {
        _accept_tests.push_back(accept_test);
    }
    void add_conf_test(std::shared_ptr<ConfTest> conf_test)
    {
        _conf_tests.push_back(conf_test);
    }
    void add_late_conf_test(std::shared_ptr<ConfTest> conf_test)
    {
        _late_conf_tests.push_back(conf_test);
    }
    void set_takestep(std::shared_ptr<TakeStep> takestep)
    {
        _take_step = takestep;
    }
    void set_coordinates(pele::Array<double>& coords, double energy)
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
    pele::Array<double> get_coords() const { return _coords.copy(); }
    pele::Array<double> get_trial_coords() const { return _trial_coords.copy(); }
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
    bool take_step_specified() const { return (_steps.size() > 0); }
    void check_input();
    void set_print_progress(const bool input) { _print_progress=input; }
    void set_print_progress() { set_print_progress(true); }

protected:
    inline double compute_energy(pele::Array<double> x)
    {
        ++_neval;
        return _potential->get_energy(x);
    }
    bool do_conf_tests(pele::Array<double> x);
    bool do_accept_tests(pele::Array<double> xtrial, double etrial, pele::Array<double> xold, double eold);
    bool do_late_conf_tests(pele::Array<double> x);
    void do_actions(pele::Array<double> x, double energy, bool success);
    void take_steps();
};

}//namespace mcpele

#endif//#ifndef _MCPELE_MC_H
