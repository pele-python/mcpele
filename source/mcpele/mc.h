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
    virtual void report(const MC* mc) {}
    virtual void increase_acceptance() {}
    virtual void decrease_acceptance() {}
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
    std::shared_ptr<pele::BasePotential> m_potential;
    pele::Array<double> m_coords;
    pele::Array<double> m_trial_coords;
    actions_t m_actions;
    accept_t m_accept_tests;
    conf_t m_conf_tests;
    conf_t m_late_conf_tests;
    std::shared_ptr<TakeStep> m_take_step;
    size_t m_nitercount;
    size_t m_accept_count;
    size_t m_E_reject_count;
    size_t m_conf_reject_count;
    bool m_success;
    /*nitercount is the cumulative count, it does not get reset at the end of run*/
    bool m_print_progress;
public:
    /*need to keep these public to make them accessible to tests and actions, be careful though!*/
    size_t m_niter;
    size_t m_neval;
    double m_stepsize;
    double m_temperature;
    double m_energy;
    double m_trial_energy;
    MC(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>& coords, const double temperature);
    virtual ~MC() {}
    void one_iteration();
    void run(const size_t max_iter);
    void set_temperature(const double T) { m_temperature = T; }
    double get_temperature() const { return m_temperature; }
    void set_stepsize(const double stepsize) { m_stepsize = stepsize; }
    void add_action(std::shared_ptr<Action> action) { m_actions.push_back(action); }
    void add_accept_test(std::shared_ptr<AcceptTest> accept_test) { m_accept_tests.push_back(accept_test); }
    void add_conf_test(std::shared_ptr<ConfTest> conf_test) { m_conf_tests.push_back(conf_test); }
    void add_late_conf_test(std::shared_ptr<ConfTest> conf_test) { m_late_conf_tests.push_back(conf_test); }
    void set_takestep(std::shared_ptr<TakeStep> takestep) { m_take_step = takestep; }
    void set_coordinates(pele::Array<double>& coords, double energy);
    double get_energy() const { return m_energy; }
    void reset_energy();
    double get_trial_energy() const { return m_trial_energy; }
    pele::Array<double> get_coords() const { return m_coords.copy(); }
    pele::Array<double> get_trial_coords() const { return m_trial_coords.copy(); }
    double get_norm_coords() const { return norm(m_coords); }
    size_t get_naccept() const { return m_accept_count; }
    size_t get_nreject() const { return m_nitercount - m_accept_count; }
    double get_accepted_fraction() const { return static_cast<double>(m_accept_count) /
            static_cast<double>(m_nitercount); }
    double get_conf_rejection_fraction() const { return static_cast<double>(m_conf_reject_count) /
            static_cast<double>(m_nitercount); }
    double get_E_rejection_fraction() const { return static_cast<double>(m_E_reject_count) /
            static_cast<double>(m_nitercount); }
    size_t get_iterations_count() const { return m_nitercount; }
    size_t get_neval() const { return m_neval; }
    double get_stepsize() const { return m_stepsize; }
    std::shared_ptr<pele::BasePotential> get_potential_ptr() { return m_potential; }
    bool take_step_specified() const { return (m_take_step != NULL); }
    void check_input();
    void set_print_progress(const bool input) { m_print_progress = input; }
    void set_print_progress() { set_print_progress(true); }
    bool get_success() const { return m_success; }
protected:
    inline double compute_energy(pele::Array<double> x)
    {
        ++m_neval;
        return m_potential->get_energy(x);
    }
    bool do_conf_tests(pele::Array<double> x);
    bool do_accept_tests(pele::Array<double> xtrial, double etrial, pele::Array<double> xold, double eold);
    bool do_late_conf_tests(pele::Array<double> x);
    void do_actions(pele::Array<double> x, double energy, bool success);
    void take_steps();
};

}//namespace mcpele

#endif//#ifndef _MCPELE_MC_H
