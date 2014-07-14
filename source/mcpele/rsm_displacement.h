#ifndef _MCPELE_RSM_DISPLACEMENT_H
#define _MCPELE_RSM_DISPLACEMENT_H

#include "pele/array.h"

namespace mcpele{


class GetMeanRMSDisplacement{
private:
    pele::Array<double> m_initial_coordinates;
    const size_t m_boxdimension;
    const size_t m_nr_particles;
public:
    virtual ~GetMeanRMSDisplacement(){}
    GetMeanRMSDisplacement(pele::Array<double>, const size_t);
    double compute_mean_rsm_displacement(pele::Array<double>);
    double get_particle_rsm_displ(const size_t, pele::Array<double>);
};


} //namespace mcpele

#endif //#ifndef _MCPELE_RSM_DISPLACEMENT_H
