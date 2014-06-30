from libcpp cimport bool as cbool
from _pele_mc cimport cppAcceptTest,_Cdef_AcceptTest

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef extern from "mcpele/accept_test.h" namespace "mcpele":
    cdef cppclass cppMetropolisTest "mcpele::MetropolisTest":
        cppMetropolisTest(size_t) except +