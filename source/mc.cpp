#include "mcpele/mc.h"
#include "mcpele/progress.h"

using pele::Array;

namespace mcpele {

MC::MC(std::shared_ptr<pele::BasePotential> potential, Array<double>& coords, const double temperature)
    : m_potential(potential),
      m_coords(coords.copy()),
      m_trial_coords(m_coords.copy()),
      m_take_step(NULL),
      m_nitercount(0),
      m_accept_count(0),
      m_E_reject_count(0),
      m_conf_reject_count(0),
      m_success(true),
      m_print_progress(false),
      m_niter(0),
      m_neval(0),
      m_temperature(temperature),
      m_report_steps(0),
      m_enable_input_warnings(true)
{
    m_energy = compute_energy(m_coords);
    m_trial_energy = m_energy;
    /*std::cout<<"mcrunner Energy is "<<_energy<< "\n";
    std::cout<<"mcrunner potential ptr is "<<_potential<< "\n";*/
}

/**
 * perform the configuration tests.  Stop as soon as one of them fails
 */
bool MC::do_conf_tests(Array<double> & x)
{
    bool result;
    for (auto & test : m_conf_tests) {
        result = test->conf_test(x, this);
        if (not result) {
            ++m_conf_reject_count;
            return false;
        }
    }
    return true;
}

/**
 * perform the acceptance tests.  Stop as soon as one of them fails
 */
bool MC::do_accept_tests(Array<double> & xtrial, double etrial, Array<double> & xold, double eold)
{
    bool result;
    for (auto & test : m_accept_tests) {
        result = test->test(xtrial, etrial, xold, eold, m_temperature, this);
        if (not result) {
            ++m_E_reject_count;
            return false;
        }
    }
    return true;
}

/**
 * perform the configuration tests.  Stop as soon as one of them fails
 */
bool MC::do_late_conf_tests(Array<double> & x)
{
    bool result;
    for (auto & test : m_late_conf_tests) {
        result = test->conf_test(x, this);
        if (not result) {
            ++m_conf_reject_count;
            return false;
        }
    }
    return true;
}

void MC::do_actions(Array<double> & x, double energy, bool success)
{
    for (auto & action : m_actions) {
        action->action(x, energy, success, this);
    }
}

void MC::take_steps()
{
    m_take_step->displace(m_trial_coords, this);
}


void MC::one_iteration()
{
    m_success = true;
    ++m_niter;
    ++m_nitercount;

    m_trial_coords.assign(m_coords);

    // take a step with the trial coords
    //_takestep->takestep(_trial_coords, _stepsize, this);
    take_steps();

    // perform the initial configuration tests
    m_success = do_conf_tests(m_trial_coords);

    // if the trial configuration is OK, compute the energy, and run the acceptance tests
    if (m_success) {
        // compute the energy
        m_trial_energy = compute_energy(m_trial_coords);

        // perform the acceptance tests.  Stop as soon as one of them fails
        m_success = do_accept_tests(m_trial_coords, m_trial_energy, m_coords, m_energy);
    }

    // Do some final checks to ensure the configuration is OK.
    // These come last because they might be computationally demanding.
    if (m_success) {
        m_success = do_late_conf_tests(m_trial_coords);
    }

    // adapt stepsize etc.
    if (get_iterations_count() <= m_report_steps) {
        m_take_step->report(m_coords, m_energy, m_trial_coords, m_trial_energy, m_success, this);
    }

    // if the step is accepted, copy the coordinates and energy
    if (m_success) {
        m_coords.assign(m_trial_coords);
        m_energy = m_trial_energy;
        ++m_accept_count;
    }

    // perform the actions on the new configuration
    do_actions(m_coords, m_energy, m_success);
}

void MC::check_input()
{
    if (!take_step_specified()) {
        throw std::runtime_error("MC::check_input: takestep not set");
    }
    if (m_enable_input_warnings) {
        if (m_conf_tests.size()==0 && m_late_conf_tests.size()==0) {
        std::cout << "warning: no conf tests set" <<"\n";
        }
        if (m_actions.size()==0) {
            std::cout << "warning: no actions set" <<  "\n";
        }
        if (m_accept_tests.size()==0) {
            std::cout << "warning: no accept tests set" <<  "\n";
        }
    }
}

void MC::set_coordinates(pele::Array<double>& coords, double energy)
{
    m_coords = coords.copy();
    m_energy = energy;
}

//this function is necessary if for example some potential parameter has been varied
void MC::reset_energy()
{
    if(m_niter > 0) {
        throw std::runtime_error("MC::reset_energy after first iteration is forbidden");
    }
    m_energy = compute_energy(m_coords);
}

void MC::run(size_t max_iter)
{
    check_input();
    progress stat(max_iter);
    while(m_niter < max_iter) {
        this->one_iteration();
        if (m_print_progress) {
            stat.next(m_niter);
        }
    }
    m_niter = 0;
}

} // namespace mcpele
