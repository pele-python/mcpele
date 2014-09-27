# distutils: language = c++
# distutils: sources = takestep.cpp

import sys

#===============================================================================
# RandomCoordsDisplacement
#===============================================================================

cdef class _Cdef_RandomCoordsDisplacement(_Cdef_TakeStep):
    """This class is the python interface for the c++ pele::RandomCoordsDisplacement take step class implementation
    """
    cdef cppRandomCoordsDisplacement* newptr
    def __cinit__(self, rseed, stepsize):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppRandomCoordsDisplacement(rseed, stepsize))
        self.newptr = <cppRandomCoordsDisplacement*> self.thisptr.get()
    
    def get_seed(self):
        cdef res = self.newptr.get_seed()
        return res
    
    def set_generator_seed(self, input):
        cdef inp = input
        self.newptr.set_generator_seed(inp)
        
class RandomCoordsDisplacement(_Cdef_RandomCoordsDisplacement):
    """This class is the python interface for the c++ RandomCoordsDisplacement implementation.
    """

#===============================================================================
# GaussianCoordsDisplacement
#===============================================================================

cdef class _Cdef_GaussianCoordsDisplacement(_Cdef_TakeStep):
    """This class is the python interface for the c++ pele::RandomCoordsDisplacement take step class implementation
    """
    cdef cppGaussianCoordsDisplacement* newptr
    def __cinit__(self, rseed, stepsize):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppGaussianCoordsDisplacement(rseed, stepsize))
        self.newptr = <cppGaussianCoordsDisplacement*> self.thisptr.get()
    
    def get_seed(self):
        cdef res = self.newptr.get_seed()
        return res
    
    def set_generator_seed(self, input):
        cdef inp = input
        self.newptr.set_generator_seed(inp)
    
class GaussianCoordsDisplacement(_Cdef_GaussianCoordsDisplacement):
    """This class is the python interface for the c++ RandomCoordsDisplacement implementation.
    """
    
#
# ParticlePairSwap
#

cdef class _Cdef_ParticlePairSwap(_Cdef_TakeStep):
    """
    Python interface for c++ ParticlePairSwap
    
    Parameters
    ----------
    seed : pos integer
        Seed for random number generator.
    nr_particles : pos integer
        Number of particles.
    swap_every : pos integer
        Spacing for swapping attempts: particle pair swap is attemted every
        swap_every move.
    """
    cdef cppParticlePairSwap* newptr
    def __cinit__(self, seed, nr_particles, swap_every):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppParticlePairSwap(seed, nr_particles, swap_every))
        self.newptr = <cppParticlePairSwap*> self.thisptr.get()

    def get_seed(self):
        cdef res = self.newptr.get_seed()
        return res
    
    def set_generator_seed(self, input):
        cdef inp = input
        self.newptr.set_generator_seed(inp)

class ParticlePairSwap(_Cdef_ParticlePairSwap):
    """
    Python interface for c++ ParticlePairSwap
    """
