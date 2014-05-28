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
	std::mt19937_64 _generator;
	std::uniform_real_distribution<double> _distribution;
public:
	RandomCoordsDisplacement(size_t rseed);
	virtual ~RandomCoordsDisplacement() {}
	virtual void takestep(Array<double>& coords, double stepsize, MC * mc);
	size_t get_seed(){return _seed;}
};

RandomCoordsDisplacement::RandomCoordsDisplacement(size_t rseed):
		_seed(rseed), _generator(_seed), _distribution(0.0,1.0)
		{
        #ifdef DEBUG
			std::cout<<"seed TakeStep:"<<_seed<<std::endl;
        #endif
		}

void RandomCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc){
	double rand;
	//assert(coords.size() == _N);
	for(size_t i=0; i<coords.size();++i){
	    rand = _distribution(_generator);
		coords[i] += (0.5-rand)*stepsize;
	}
}

/*
 * Uniform Gaussian step
 * this step samples first from the standard normal N(0,1) and outputs a random variate sampled from N(0,stepsize)
 */

class GaussianCoordsDisplacement:public TakeStep{
protected:
    size_t _seed;
    double _mean, _stdev;
    std::mt19937_64 _generator;
    std::normal_distribution<double> _distribution;
public:
    GaussianCoordsDisplacement(size_t rseed);
    virtual ~GaussianCoordsDisplacement() {}
    virtual void takestep(Array<double>& coords, double stepsize, MC * mc);
    size_t get_seed(){return _seed;}
};

GaussianCoordsDisplacement::GaussianCoordsDisplacement(size_t rseed):
        _seed(rseed), _mean(0.0), _stdev(1.0),
        _generator(_seed), _distribution(_mean,_stdev)
        {
        #ifdef DEBUG
            std::cout<<"seed TakeStep:"<<_seed<<std::endl;
        #endif
        }

void GaussianCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc){
    //assert(coords.size() == _N);
    for(size_t i=0; i<coords.size();++i){
        double randz = _distribution(_generator); //this is sample from N(0,1)
        coords[i] += randz*stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
}

}
#endif
