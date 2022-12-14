#include "mcpele/mc.h"




namespace mcpele {

void TakeStep::set_current_step_name(MC *mc, std::string prefix) {
        if (!mc->m_record_acceptance_rate) {
            return;
        }
        std::string name = demangle(typeid(*this).name());
        if (prefix != "") {
            name = prefix + "_" + name;
        }
        mc->add_step_name_to_success_accumulator(name);
}
}