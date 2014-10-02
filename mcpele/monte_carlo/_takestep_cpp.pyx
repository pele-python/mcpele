# distutils: language = c++
# distutils: sources = takestep.cpp

import sys

#===============================================================================
# RandomCoordsDisplacement
#===============================================================================

cdef class _Cdef_RandomCoordsDisplacement(_Cdef_TakeStep):
    """This class is the python interface for the c++ pele::RandomCoordsDisplacement take step class implementation
    """
    cdef cppRandomCoordsDisplacementAdaptive* newptr
    def __cinit__(self, rseed, stepsize, report_interval=100, factor=0.9, min_acc_ratio=0.2, max_acc_ratio=0.5):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppAdaptiveTakeStep(shared_ptr[cppTakeStep](<cppTakeStep*> new cppRandomCoordsDisplacementAdaptive(rseed, stepsize)), report_interval, factor, min_acc_ratio, max_acc_ratio))
        self.newptr = <cppRandomCoordsDisplacementAdaptive*> self.thisptr.get()
    
    def get_seed(self):
        cdef res = self.newptr.get_seed()
        return res
    
    def set_generator_seed(self, input):
        cdef inp = input
        self.newptr.set_generator_seed(inp)
        
    def get_stepsize(self):
        return self.newptr.get_stepsize()
        
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
        
    def get_count(self):
        return self.newptr.get_count()
    
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
    def __cinit__(self, seed, nr_particles):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppParticlePairSwap(seed, nr_particles))
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
    
#
# TakeStepPattern
#

cdef class _Cdef_TakeStepPattern(_Cdef_TakeStep):
    """
    Python interface for c++ TakeStepPattern
    """
    cdef cppTakeStepPattern* newptr
    def __cinit__(self):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppTakeStepPattern())
        self.newptr = <cppTakeStepPattern*> self.thisptr.get()
    
    def add_step(self, _Cdef_TakeStep step, nr_repetitions):
        self.newptr.add_step(step.thisptr, nr_repetitions)

class TakeStepPattern(_Cdef_TakeStepPattern):
    """
    Python interface for c++ TakeStepPattern
    """
    
#
# TakeStepProbabilities
#

cdef class _Cdef_TakeStepProbabilities(_Cdef_TakeStep):
    """
    Python interface for c++ TakeStepProbabilities
    """
    cdef cppTakeStepProbabilities* newptr
    def __cinit__(self, seed):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppTakeStepProbabilities(seed))
        self.newptr = <cppTakeStepProbabilities*> self.thisptr.get()
    def add_step(self, _Cdef_TakeStep step, weight):
        self.newptr.add_step(step.thisptr, weight)
        
class TakeStepProbabilities(_Cdef_TakeStepProbabilities):
    """
    Python interface for c++ TakeStepProbabilities
    """
