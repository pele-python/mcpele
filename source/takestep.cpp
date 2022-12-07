#include "mcpele/mc.h"




namespace mcpele {

void TakeStep::set_current_step_name(MC *mc) {
        std::string name = demangle(typeid(*this).name());
        mc->add_step_name_to_success_accumulator(name);
}
}