#include "mcpele/check_spherical_container.h"

using std::runtime_error;
using pele::Array;

namespace mcpele{

CheckSphericalContainer::CheckSphericalContainer(const double radius, const size_t ndim)
    : m_radius2(radius * radius),
      m_ndim(ndim)
{}

bool CheckSphericalContainer::conf_test(Array<double> &trial_coords, MC * mc)
{
  const size_t N = trial_coords.size();
  for (size_t i = 0; i < N; i += m_ndim) {
      double r2 = 0;
      for (size_t j = i; j < i + m_ndim; ++j) {
          r2 += trial_coords[j] * trial_coords[j];
      }
      if (r2 > m_radius2) {
          //printf("fail spherical container %d %f %f %f %f\n", i, sqrt(r2), x[i], x[i+1], x[i+2]);
          //an atom is outside the spherical container
          return false;
      }
  }
  //printf("check spherical OK ");
  return true;
}

} // namespace mcpele
