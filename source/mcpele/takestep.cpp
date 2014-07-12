#include "takestep.h"

namespace mcpele{

RandomCoordsDisplacement::RandomCoordsDisplacement(size_t rseed)
    : _seed(rseed),
      _generator(rseed),
      _distribution(0.0, 1.0)
{
    #ifdef DEBUG
        std::cout<<"seed TakeStep:"<<_seed<< "\n";
    #endif
}

void RandomCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc)
{
    double rand;
    //assert(coords.size() == _N);
    for (size_t i=0; i<coords.size();++i) {
        rand = _distribution(_generator);
        coords[i] += (0.5-rand)*stepsize;
    }
}

GaussianCoordsDisplacement::GaussianCoordsDisplacement(size_t rseed)
    : _seed(rseed),
      _mean(0.0),
      _stdev(1.0),
      _generator(rseed),
      _distribution(_mean, _stdev)
{
    #ifdef DEBUG
        std::cout<<"seed TakeStep:"<<_seed<< "\n";
    #endif
}

void GaussianCoordsDisplacement::takestep(Array<double>& coords, double stepsize, MC * mc)
{
    //assert(coords.size() == _N);
    for(size_t i=0; i<coords.size(); ++i){
        double randz = _distribution(_generator); //this is sample from N(0,1)
        coords[i] += randz * stepsize; //here the stepsize plays the same role as the stdev. This is sampled from N(0,stepsize)
    }
}

}//namespace mcpele
