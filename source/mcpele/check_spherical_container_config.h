#ifndef _MCPELE_CHECK_SPHERICAL_CONTAINER_CONFIG_H__
#define _MCPELE_CHECK_SPHERICAL_CONTAINER_CONFIG_H__

#include "mc.h"

namespace mcpele {

class CheckSphericalContainerConfig : public ConfTest {
protected:
    double m_radius2;
public:
    CheckSphericalContainerConfig(const double radius) : m_radius2(radius * radius) {}
    bool conf_test(pele::Array<double> &trial_coords, MC * mc) { return pele::dot(trial_coords, trial_coords) <= m_radius2; }
    virtual ~CheckSphericalContainerConfig() {}
};

} // namespace mcpele

#endif // #ifndef _MCPELE_CHECK_SPHERICAL_CONTAINER_CONFIG_H__
