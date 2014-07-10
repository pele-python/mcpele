#include <cmath>
#include <stdexcept>

#include "histogram.h"
#include "rsm_displacement.h"

namespace mcpele{


GetMeanRMSDisplacement::GetMeanRMSDisplacement(pele::Array<double> initial_coordinates_,
        const size_t boxdimension_)
    : initial_coordinates(initial_coordinates_.copy()),
      boxdimension(boxdimension_),
      nr_particles(initial_coordinates_.size() / boxdimension_)
{}

double GetMeanRMSDisplacement::compute_mean_rsm_displacement(pele::Array<double> new_coords)
{
    if (new_coords.size() != initial_coordinates.size()) {
        throw std::runtime_error("GetMeanRMSDisplacement::compute_mean_rsm_displacement: illegal new coords");
    }
    Moments mom;
    for (size_t particle = 0; particle < nr_particles; ++particle) {
        mom.update(get_particle_rsm_displ(particle, new_coords));
    }
    return mom.mean();
}

double GetMeanRMSDisplacement::get_particle_rsm_displ(const size_t particle_idx, pele::Array<double> new_coords)
{
    const size_t particle_start = particle_idx * boxdimension;
    double sumsq(0);
    for (size_t i = particle_start; i < particle_start + boxdimension; ++i) {
        const double deltai = initial_coordinates[i] - new_coords[i];
        sumsq += deltai*deltai;
    }
    return sqrt(sumsq);
}


} //namespace mcpele
