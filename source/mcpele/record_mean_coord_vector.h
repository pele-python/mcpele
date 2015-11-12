#ifndef _MCPELE_RECORD_MEAN_COORD_VECTOR_H
#define _MCPELE_RECORD_MEAN_COORD_VECTOR_H

#include "mc.h"
#include "pele/array.h"

namespace mcpele{


class RecordMeanCoordVector : public Action{
private:
    pele::Array<double> m_mcv;
    const size_t m_ndof, m_count, m_eqsteps;
    double _update_average(double avg, double x);
    void _update_mean_coord_vector(pele::Array<double> new_coords);
public:
    RecordMeanCoordVector(pele::Array<double> initial_coords, const size_t eqsteps);
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
    virtual ~RecordMeanCoordVector(){}
    pele::Array<double> get_mean_coordinate_vector(){return m_mcv.copy();}
};


} //namespace mcpele

#endif //#ifndef _MCPELE_RECORD_MEAN_COORD_VECTOR_H
