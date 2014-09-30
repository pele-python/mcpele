#ifndef _MCPELE_RSM_DISPLACEMENT_H
#define _MCPELE_RSM_DISPLACEMENT_H

#include "pele/array.h"

namespace mcpele{


class GetDisplacementPerParticle{
private:
    pele::Array<double> m_initial_coordinates;
    const size_t m_boxdimension;
    const size_t m_nr_particles;
public:
    virtual ~GetDisplacementPerParticle(){}
    GetDisplacementPerParticle(pele::Array<double>, const size_t);
    double compute_mean_particle_displacement(pele::Array<double>);
    double get_particle_displ(const size_t, pele::Array<double>);
};


} //namespace mcpele

#endif //#ifndef _MCPELE_RSM_DISPLACEMENT_H
