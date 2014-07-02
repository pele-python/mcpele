# distutils: language = c++
# distutils: sources = accept_test.cpp

import sys

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

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