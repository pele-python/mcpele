# distutils: language = c++
# distutils: sources = accept_test.cpp

cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport cppAcceptTest,_Cdef_AcceptTest

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef extern from "mcpele/accept_test.h" namespace "mcpele":
    cdef cppclass cppMetropolisTest "mcpele::MetropolisTest":
        cppMetropolisTest(size_t) except +

cdef class _Cdef_Metropolis(_Cdef_AcceptTest):
    """This class is the python interface for the c++ pele::MetropolisTest acceptance test class implementation
    """
    cdef cppMetropolisTest* newptr
    def __cinit__(self, rseed):
        self.thisptr = <cppAcceptTest*>new cppMetropolisTest(rseed)
        self.newptr = <cppMetropolisTest*> self.thisptr

    def __dealloc__(self):
        del self.thisptr
        
class MetropolisTest(_Cdef_Metropolis):
    """This class is the python interface for the c++ Metropolis implementation.
    """