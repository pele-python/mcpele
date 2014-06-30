#include "mc.h"

namespace mcpele{

MC::MC(pele::BasePotential * potential, Array<double>& coords, double temperature, double stepsize):
            _potential(potential), _coords(coords.copy()),_trial_coords(_coords.copy()), _takestep(NULL),
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

    //for (auto& test : _conf_tests ){
    for (conf_t::iterator test1 = _conf_tests.begin(); test1 != _conf_tests.end(); ++test1){
    //for (size_t test1 = 0; test1 != _conf_tests.size(); ++test1){
	//_success = test->test(_trial_coords, this);
	_success = (*test1)->test(_trial_coords, this);
	//_success = _conf_tests.at(test1)->test(_trial_coords, this);
	    if (_success == false){
		    ++_conf_reject_count;
		    break;
	    }
    }

    if (_success == true)
    {
	    _trial_energy = _potential->get_energy(_trial_coords);
	    ++_neval;

	    for (accept_t::iterator test2 = _accept_tests.begin(); test2 != _accept_tests.end(); ++test2){
	    //for (size_t test2 = 0; test2 < _accept_tests.size(); ++test2){
		    _success = (*test2)->test(_trial_coords, _trial_energy, _coords, _energy, _temperature, this);
		    //_success = _accept_tests.at(test2)->test(_trial_coords, _trial_energy, _coords, _energy, _temperature, this);
		    if (_success == false){
			    ++_E_reject_count;
			    break;
		    }
	    }

	    if (_success == true){

		for (conf_t::iterator test3 = _late_conf_tests.begin(); test3 != _late_conf_tests.end(); ++test3){
		//for (size_t test3 = 0; test3 < _late_conf_tests.size(); ++test3){
		    _success = (*test3)->test(_trial_coords, this);
		    //_success = _late_conf_tests.at(test3)->test(_trial_coords, this);
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
    for (actions_t::iterator action1 = _actions.begin(); action1 != _actions.end(); ++action1){
    //for (size_t action1 = 0; action1 < _actions.size(); ++action1){
	(*action1)->action(_coords, _energy, _success, this);
	//_actions.at(action1)->action(_coords, _energy, _success, this);
    }
}

void MC::check_input(){
    if (!take_step_specified())
	throw std::runtime_error("MC::check_input: takestep not set");
}

void MC::run(size_t max_iter)
{
    std::cout << "_conf_tests.size(): " << _conf_tests.size() << std::endl; //debug
    std::cout << "_late_conf_tests.size(): " << _late_conf_tests.size() << std::endl; //debug
    std::cout << "_actions.size(): " << _actions.size() << std::endl; //debug
    std::cout << "_accept_tests.size(): " << _accept_tests.size() << std::endl; //debug
    //throw std::runtime_error("TEST!");
    check_input();
    while(_niter < max_iter)
	    //std::cout << "done: " << double(_niter)/double(max_iter) << std::endl;
	    this->one_iteration();
    _niter = 0;
}

}//namespace mcpele
