# distutils: language = c++
# distutils: sources = accept_test.cpp

import sys

#===============================================================================
# Metropolis acceptance criterion
#===============================================================================

cdef class _Cdef_Metropolis(_Cdef_AcceptTest):
    cdef cppMetropolisTest* newptr
    def __cinit__(self, rseed):
        self.thisptr = shared_ptr[cppAcceptTest](<cppAcceptTest*> new cppMetropolisTest(rseed))
        self.newptr = <cppMetropolisTest*> self.thisptr.get()
    
    def get_seed(self):
        """return random number generator seed
        
        Returns
        -------
        int 
            random number generator seed
        """
        cdef res = self.newptr.get_seed()
        return res
    
    def set_generator_seed(self, input):
        """sets the random number generator seed
        
        Parameters
        ----------
        input : pos int
            random number generator seed
        """
        cdef inp = input
        self.newptr.set_generator_seed(inp)
        
class MetropolisTest(_Cdef_Metropolis):
    """Metropolis acceptance criterion
    
    This class is the Python interface for the c++ pele::MetropolisTest 
    acceptance test class implementation. The Metropolis acceptance criterion
    accepts each move with probability
    
    .. math:: P( x_{old} \Rightarrow x_{new}) = min \{ 1, \exp [- \\beta (E_{new} - E_{old})] \}
    
    where :math:`\\beta` is the reciprocal of the temperature.
    """