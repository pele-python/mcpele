#ifndef _MCPELE_ACCEPT_TEST_H__
#define _MCPELE_ACCEPT_TEST_H__

#include <limits>
#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "pele/array.h"
#include "mc.h"


using std::runtime_error;
using pele::Array;

namespace mcpele{

/*Metropolis acceptance criterion
 * */

class MetropolisTest:public AcceptTest{
protected:
	size_t _seed;
	std::mt19937_64 _generator;
	std::uniform_real_distribution<double> _distribution;

public:
	MetropolisTest(size_t rseed);
	virtual ~MetropolisTest() {}
	virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc);
	size_t get_seed() const {return _seed;}
};

/*ENERGY WINDOW TEST
 * this test checks that the energy of the system stays within a certain energy window
 * */

class EnergyWindowTest:public AcceptTest{
protected:
	double _min_energy, _max_energy;
public:
	EnergyWindowTest(double min_energy, double max_energy);
	virtual ~EnergyWindowTest() {}
	virtual bool test(Array<double> &trial_coords, double trial_energy, Array<double> & old_coords, double old_energy, double temperature, MC * mc);
};

}
#endif
