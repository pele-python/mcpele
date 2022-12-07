#ifndef _MCPELE_MC_H
#define _MCPELE_MC_H

#include <cmath>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <typeinfo>
#include <cxxabi.h>

#include "pele/array.hpp"
#include "pele/base_potential.hpp"
#include "success_container.h"
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



inline std::string demangle(const char* mangled)
{
      int status;
      std::unique_ptr<char[], void (*)(void*)> result(
        abi::__cxa_demangle(mangled, 0, 0, &status), std::free);
      return result.get() ? std::string(result.get()) : "error occurred";
}


/*
 * Take Step
 */

class TakeStep {
public:
    virtual ~TakeStep() {}
    virtual void displace(pele::Array<double>& coords, MC* mc) = 0;
    virtual void report(pele::Array<double>&, const double,
            pele::Array<double>&, const double, const bool, MC*) {}
    virtual void increase_acceptance(const double) {}
    virtual void decrease_acceptance(const double) {}
    virtual const std::vector<long> get_changed_atoms() const { return std::vector<long>(); }
    virtual const std::vector<double> get_changed_coords_old() const { return std::vector<double>(); }
    virtual void set_current_step_name(MC * mc);

};

/**
 * Monte Carlo
 * _coords and _trialcoords are arrays that store coordinates and trial coordinates respectively
 * _potential is on object of Pele::BasePotential type that defines the interaction potential
 * _E_reject_count is the count of rejections due to an energy test (e.g. Metropolis)
 * _conf_reject_count is the count of rejections due to a configuration test (e.g. spherical container)
 * _niter is the count of steps whithin a MCMC run, it is reset to zero at the end of the run
 * _nitercount is the cumulative number of MCMC steps taken by the class
 * _neval is the number of energy evaluations
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
    SuccessAccumulator m_success_accumulator;
    size_t m_nitercount;
    size_t m_accept_count;
    size_t m_E_reject_count;
    size_t m_conf_reject_count;
    bool m_success, m_last_success;
    bool m_print_progress;
public:
    /*need to keep these public to make them accessible to tests and actions, be careful though!*/
    /*nitercount is the cumulative count, it does not get reset at the end of run*/
    size_t m_niter;
    size_t m_neval;
    double m_temperature;
    double m_energy;
    double m_trial_energy;
private:
    size_t m_report_steps;
    bool m_enable_input_warnings;
public:
    MC(std::shared_ptr<pele::BasePotential> potential, pele::Array<double>& coords, const double temperature);
    virtual ~MC() {}
    void one_iteration();

    void run(const size_t max_iter);
    void set_temperature(const double T) { m_temperature = T; }
    double get_temperature() const { return m_temperature; }
    void set_report_steps(const size_t report_steps) { m_report_steps = report_steps; }
    size_t get_report_steps() const { return m_report_steps; }
    void add_action(std::shared_ptr<Action> action) { m_actions.push_back(action); }
    void add_accept_test(std::shared_ptr<AcceptTest> accept_test) { m_accept_tests.push_back(accept_test); }
    void add_conf_test(std::shared_ptr<ConfTest> conf_test) { m_conf_tests.push_back(conf_test); }
    void add_late_conf_test(std::shared_ptr<ConfTest> conf_test) { m_late_conf_tests.push_back(conf_test); }
    void set_takestep(std::shared_ptr<TakeStep> takestep) { m_take_step = takestep; }
    std::shared_ptr<TakeStep> get_takestep() const { return m_take_step; }
    void set_coordinates(pele::Array<double>& coords, double energy);
    double get_energy() const { return m_energy; }
    void reset_energy();
    double get_trial_energy() const { return m_trial_energy; }
    pele::Array<double> get_coords() const { return m_coords.copy(); }
    pele::Array<double> get_trial_coords() const { return m_trial_coords.copy(); }
    double get_norm_coords() const { return norm(m_coords); }
    size_t get_naccept() const { return m_accept_count; }
    size_t get_nreject() const { return m_nitercount - m_accept_count; }

    // get acceptance and rejection fractions
    double get_accepted_fraction() const { return static_cast<double>(m_accept_count) /
            static_cast<double>(m_nitercount); }
    double get_conf_rejection_fraction() const { return static_cast<double>(m_conf_reject_count) /
            static_cast<double>(m_nitercount); }
    double get_E_rejection_fraction() const { return static_cast<double>(m_E_reject_count) /
            static_cast<double>(m_nitercount); }

    double get_norm_trial_coords() const { return norm(m_trial_coords); }
    

    size_t get_iterations_count() const { return m_nitercount; }
    size_t get_neval() const { return m_neval; }
    std::shared_ptr<pele::BasePotential> get_potential_ptr() { return m_potential; }
    bool take_step_specified() const { return m_take_step != NULL; }
    bool report_steps_specified() const { return get_report_steps() > 0; }
    void check_input();

    void set_print_progress(const bool input) { m_print_progress = input; }
    void set_print_progress() { set_print_progress(true); }

    // Assuming that these are running one iteration
    bool get_success() const { return m_success_accumulator.get_current_success(); }
    bool get_last_success() const { return m_success_accumulator.get_last_success(); }

    // this helps map the step name to the success accumulator
    void add_step_name_to_success_accumulator(const std::string& name) { m_success_accumulator.add_step_taken(name); }



    pele::Array<size_t> get_counters() const
    {
        pele::Array<size_t> counters(5);
        counters[0] = m_nitercount;
        counters[1] = m_accept_count;
        counters[2] = m_E_reject_count;
        counters[3] = m_conf_reject_count;
        counters[4] = m_neval;
        return counters;
    }
    void set_counters(pele::Array<size_t> const & counters)
    {
        m_nitercount = counters[0];
        m_accept_count = counters[1];
        m_E_reject_count = counters[2];
        m_conf_reject_count = counters[3];
        m_neval = counters[4];
    }
    const std::vector<long> get_changed_atoms() const { return m_take_step->get_changed_atoms(); }
    const std::vector<double> get_changed_coords_old() const { return m_take_step->get_changed_coords_old(); }
    /**
     * this will trigger premature exit from the MC run loop
     */
    void abort() { m_niter = std::numeric_limits<size_t>::max(); }
    void enable_input_warnings() { m_enable_input_warnings = true; }
    void disable_input_warnings() { m_enable_input_warnings = false; }
protected:
    inline double compute_energy(pele::Array<double> & x)
    {
        ++m_neval;
        return m_potential->get_energy(x);
    }
    bool do_conf_tests(pele::Array<double> & x);
    bool do_accept_tests(pele::Array<double> & xtrial, double etrial, pele::Array<double> & xold, double eold);
    bool do_late_conf_tests(pele::Array<double> & x);
    void do_actions(pele::Array<double> & x, double energy, bool success);
    void take_steps();
};

}//namespace mcpele

#endif//#ifndef _MCPELE_MC_H
