#ifndef _MCPELE_CONF_TEST_OR_H__
#define _MCPELE_CONF_TEST_OR_H__

#include <vector>
#include <random>

#include "mc.h"

namespace mcpele {

/**
 * Create union of two configurational tests,
 * it is sufficient for one of them to be true in order to pass the overall test.
 *
 */
class ConfTestOR : public ConfTest {
private:
    std::vector<std::shared_ptr<ConfTest> > m_tests;
public:
    virtual ~ConfTestOR(){}
    ConfTestOR();
    void add_test(std::shared_ptr<ConfTest> test_input);
    bool conf_test(pele::Array<double> &trial_coords, MC * mc);
};

} // namespace mcpele

#endif // #ifndef _MCPELE_CONF_TEST_OR_H__
