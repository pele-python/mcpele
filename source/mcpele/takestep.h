#ifndef _MCPELE_TAKESTEP_H__
#define _MCPELE_TAKESTEP_H__

#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "pele/array.h"
#include "mc.h"


using std::runtime_error;
using pele::Array;
using std::sqrt;

namespace mcpele{

/*Random coords displacement, generates a random displacement for a N dimensional system
 * sampling from a N-dimensional sphere
 * the stepsize is defined per coordinates, that's why the maximum stepsize is sqrt(N)*stepsize
 * */

class RandomCoordsDisplacement:public TakeStep{
protected:
	size_t _seed;
	static std::mt19937_64 _generator;
	std::uniform_real_distribution<double> _distribution;
public:
	RandomCoordsDisplacement(size_t rseed);
	RandomCoordsDisplacement();
	virtual ~RandomCoordsDisplacement(){}
	//virtual ~RandomCoordsDisplacement() { std::cout << "destruct RandomCoordsDisplacement" << std::endl; }
	virtual void takestep(Array<double>& coords, double stepsize, MC * mc=NULL);
	size_t get_seed() const {return _seed;}
	static void set_generator_seed(const size_t inp){_generator.seed(inp);}
};

/*
 * Uniform Gaussian step
 * this step samples first from the standard normal N(0,1) and outputs a random variate sampled from N(0,stepsize)
 */

class GaussianCoordsDisplacement:public TakeStep{
protected:
    size_t _seed;
    double _mean, _stdev;
    static std::mt19937_64 _generator;
    std::normal_distribution<double> _distribution;
public:
    GaussianCoordsDisplacement(size_t rseed);
    GaussianCoordsDisplacement();
    virtual ~GaussianCoordsDisplacement(){}
    virtual void takestep(Array<double>& coords, double stepsize, MC * mc=NULL);
    size_t get_seed() const {return _seed;}
    static void set_generator_seed(const size_t inp){_generator.seed(inp);}
};

}
#endif
