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
        self.thisptr = shared_ptr[cppAcceptTest](<cppAcceptTest*> new cppMetropolisTest(rseed))
        self.newptr = <cppMetropolisTest*> self.thisptr.get()
    
    def get_seed(self):
        cdef res = self.newptr.get_seed()
        return res
    
    def set_generator_seed(self, input):
        cdef inp = input
        self.newptr.set_generator_seed(inp)
        
class MetropolisTest(_Cdef_Metropolis):
    """This class is the python interface for the c++ Metropolis implementation.
    """