#include "mc.h"

namespace mcpele{

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

}//namespace mcpele
