#ifndef _MCPELE_RECORD_MEAN_COORD_VECTOR_H__
#define _MCPELE_RECORD_MEAN_COORD_VECTOR_H__

#include "mc.h"
#include "pele/array.h"

namespace mcpele{

class RecordMeanCoordVector : public Action{
private:
    pele::Array<double> m_mcv, m_mcv2;
    const size_t m_ndof, m_eqsteps;
    size_t m_count;
    double m_update_average(double avg, double x);
    void m_update_mean_coord_vector(pele::Array<double> &new_coords);
public:
    RecordMeanCoordVector(const size_t ndof, const size_t eqsteps);
    virtual ~RecordMeanCoordVector(){}
    virtual void action(pele::Array<double> &coords, double energy, bool accepted, MC* mc);
    pele::Array<double> get_mean_coordinate_vector(){return m_mcv.copy();}
    pele::Array<double> get_mean2_coordinate_vector(){return m_mcv2.copy();}
    pele::Array<double> get_variance_coordinate_vector(){
        pele::Array<double> var = m_mcv2.copy();
        for(size_t i=0; i<m_ndof; ++i){
            var[i] -= m_mcv[i]*m_mcv[i];
        }
        return var.copy();
    }
    size_t get_count(){return m_count;}
};

} //namespace mcpele

#endif //#ifndef _MCPELE_RECORD_MEAN_COORD_VECTOR_H__
