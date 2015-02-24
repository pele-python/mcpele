from libcpp cimport bool as cbool
from _pele_mc cimport cppAcceptTest,_Cdef_AcceptTest, shared_ptr

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef extern from "mcpele/metropolis_test.h" namespace "mcpele":
    cdef cppclass cppMetropolisTest "mcpele::MetropolisTest":
        cppMetropolisTest(size_t) except +
        size_t get_seed() except +
        void set_generator_seed(size_t) except +
