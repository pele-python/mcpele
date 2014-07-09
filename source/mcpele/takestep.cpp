#include "takestep.h"

namespace mcpele{

//static members of RandomCoordsDisplacement
std::mt19937_64 RandomCoordsDisplacement::_generator;

//other members
RandomCoordsDisplacement::RandomCoordsDisplacement(size_t rseed)
    : _seed(rseed), _distribution(0.0, 1.0)
{
    set_generator_seed(_seed);
#ifdef DEBUG
    std::cout<<"seed TakeStep:"<<_seed<< "\n";
#endif
    //std::cout << "construct RandomCoordsDisplacement" <<  "\n";
}

RandomCoordsDisplacement::RandomCoordsDisplacement()
    : _seed(42), _distribution(0.0, 1.0)
{}

void RandomCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc)
{
    double rand;
    //assert(coords.size() == _N);
    for (size_t i=0; i<coords.size();++i) {
        rand = _distribution(_generator);
        coords[i] += (0.5-rand)*stepsize;
    }
}

//static members of GaussianCoordsDisplacement
std::mt19937_64 GaussianCoordsDisplacement::_generator;

//other members
GaussianCoordsDisplacement::GaussianCoordsDisplacement(size_t rseed)
    : _seed(rseed), _mean(0.0), _stdev(1.0), _distribution(_mean, _stdev)
{
    set_generator_seed(_seed);
#ifdef DEBUG
    std::cout<<"seed TakeStep:"<<_seed<< "\n";
#endif
}

GaussianCoordsDisplacement::GaussianCoordsDisplacement()
    : _seed(42), _mean(0.0), _stdev(1.0), _distribution(_mean, _stdev)
{}

void GaussianCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc)
{
    //assert(coords.size() == _N);
    for(size_t i=0; i<coords.size(); ++i){
        double randz = _distribution(_generator); //this is sample from N(0,1)
        coords[i] += randz * stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
}

}//namespace mcpele
