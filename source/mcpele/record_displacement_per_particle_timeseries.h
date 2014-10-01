#ifndef _MCPELE_RECORD_DISPLACEMENT_PER_PARTICLE_TIMESERIES_H__
#define _MCPELE_RECORD_DISPLACEMENT_PER_PARTICLE_TIMESERIES_H__

#include "record_scalar_timeseries.h"
#include "rsm_displacement.h"

namespace mcpele {

/**
 * Record time series of mean displacement per particle.
 * Motivation: check if HS fluid is decorrelated between snapshots
 * This is useful if particles are not placed back into periodic box.
 */
class RecordDisplacementPerParticleTimeseries : public RecordScalarTimeseries{
private:
    GetDisplacementPerParticle m_rsm_displacement;
public:
    RecordDisplacementPerParticleTimeseries(const size_t niter,
            const size_t record_every, pele::Array<double> initial_coords,
            const size_t boxdimension)
        : RecordScalarTimeseries(niter, record_every),
          m_rsm_displacement(initial_coords, boxdimension)
    {}
    virtual ~RecordDisplacementPerParticleTimeseries(){}
    virtual double get_recorded_scalar(pele::Array<double> &coords,
            const double energy, const bool accepted, MC* mc)
            { return m_rsm_displacement.compute_mean_particle_displacement(coords); }
};

} // namespace mcpele

#endif // #ifndef _MCPELE_RECORD_DISPLACEMENT_PER_PARTICLE_TIMESERIES_H__
