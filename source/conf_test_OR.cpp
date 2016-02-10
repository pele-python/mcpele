#include "mcpele/conf_test_OR.h"

namespace mcpele {

ConfTestOR::ConfTestOR(){}

void ConfTestOR::add_test(std::shared_ptr<ConfTest> test_input)
{
    m_tests.push_back(test_input);
    m_tests.swap(m_tests);
}

bool ConfTestOR::conf_test(pele::Array<double> &trial_coords, MC * mc)
{
    if (m_tests.size() == 0) {
        throw std::runtime_error("ConfTestOR::conf_test: no conf test specified");
    }
    for(auto & test : m_tests){
        bool result = test->conf_test(trial_coords, mc);
        if (result){
            return true;
        }
    }
    return false;
}

} // namespace mcpele
