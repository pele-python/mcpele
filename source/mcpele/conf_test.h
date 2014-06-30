#ifndef _MCPELE_CONF_TEST_H__
#define _MCPELE_CONF_TEST_H__

#include <iostream>
#include <math.h>
#include <algorithm>
#include <random>
#include <chrono>
#include "pele/array.h"
#include "pele/optimizer.h"
#include "pele/distance.h"
#include "mc.h"

using std::runtime_error;
using pele::Array;

namespace mcpele{

class CheckSphericalContainer:public ConfTest{
protected:
	double _radius2;
	size_t _ndim;
public:

	CheckSphericalContainer(double radius, size_t ndim):
			_radius2(radius*radius),_ndim(ndim)
	    {
		//std::cout << "construct CheckSphericalContainer" << std::endl;
	    }

	virtual bool test(Array<double> &trial_coords, MC * mc)
	{
	  size_t N = trial_coords.size();

	  for (size_t i=0; i<N; i+=_ndim)
	  {
	      double r2 = 0;
		  for (size_t j=i; j<i+_ndim; ++j)
		  {
			r2 += trial_coords[j] * trial_coords[j];
		  }
		  if (r2 > _radius2)
		  {
		    //printf("fail spherical container %d %f %f %f %f\n", i, sqrt(r2), x[i], x[i+1], x[i+2]);
		    //an atom is outside the spherical container
		    return false;
		  }
	  }

	  //printf("check spherical OK ");
	  return true;
	}

	virtual ~CheckSphericalContainer(){}
};

}//namespace mcpele


#endif
