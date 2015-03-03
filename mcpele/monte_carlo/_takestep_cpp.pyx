# distutils: language = c++
# distutils: sources = takestep.cpp

import sys
from pele.potentials import _pele
from pele.potentials._pele cimport array_wrap_np

#===============================================================================
# RandomCoordsDisplacement
#===============================================================================

cdef class _Cdef_RandomCoordsDisplacement(_Cdef_TakeStep):
    cdef cppRandomCoordsDisplacement* newptr
    def __cinit__(self, rseed, stepsize, report_interval=100, factor=0.9,
                  min_acc_ratio=0.2, max_acc_ratio=0.5, single=False,
                  nparticles=0, bdim=0):
        if not single:
            self.newptr = <cppRandomCoordsDisplacement*> new cppRandomCoordsDisplacementAll(rseed, stepsize)
            self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> 
                   new cppAdaptiveTakeStep(shared_ptr[cppTakeStep](<cppTakeStep*> self.newptr), 
                                           report_interval, factor, min_acc_ratio, max_acc_ratio))
        else:
            assert(bdim > 0 and nparticles > 0)
            self.newptr = <cppRandomCoordsDisplacement*> new cppRandomCoordsDisplacementSingle(rseed, nparticles, bdim, stepsize)
            self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> 
                   new cppAdaptiveTakeStep(shared_ptr[cppTakeStep](<cppTakeStep*> self.newptr), 
                                           report_interval, factor, min_acc_ratio, max_acc_ratio))
        
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
        
    def get_count(self):
        """get the total count of the number of steps taken
        
        Returns
        -------
        int
            total count of steps taken
        """
        return self.newptr.get_count()
        
    def get_stepsize(self):
        """get the step size
        
        Returns
        -------
        double
            stepsize
        """
        return self.newptr.get_stepsize()
        
class RandomCoordsDisplacement(_Cdef_RandomCoordsDisplacement):
    """Take a uniform random step in a ``bdim`` dimensional hypercube
    
    this class is the Python interface for the c++ RandomCoordsDisplacement implementation.
    Takes a step by sampling uniformly a ``bdim`` dimensional box
    
    Parameters
    ----------
    rseed : pos int
        seed for the random number generator (std:library 64 bits Merseene Twister)
    stepsize : double
        size of step in each dimension
    report_interval : int
        number of report steps for which the step size should be adapted
    factor : double
        factor by which the step size is adapted at each iteration durint the report interval
    min_acc_ratio : double
        minimum of target acceptance range
    max_acc_ratio: double
        maximum of target acceptance range
    single : bool
        True for single particle moves, False for global moves 
    nparticles : int
        number of particles, typically len(coords)/bdim
    bdim : int
        dimensionality of the space (box dimensionality)
    """

#===============================================================================
# GaussianCoordsDisplacement
#===============================================================================

cdef class _Cdef_GaussianCoordsDisplacement(_Cdef_TakeStep):
    cdef cppGaussianTakeStep* newptr
    def __cinit__(self, rseed, stepsize, ndim):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppGaussianCoordsDisplacement(rseed, stepsize, ndim))
        self.newptr = <cppGaussianTakeStep*> self.thisptr.get()
    
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
        
    def get_count(self):
        """get the total count of the number of steps taken
        
        Returns
        -------
        int
            total count of steps taken
        """
        return self.newptr.get_count()
    
    def get_stepsize(self):
        """get the step size
        
        Returns
        -------
        double
            stepsize
        """
        return self.newptr.get_stepsize()
    
class GaussianCoordsDisplacement(_Cdef_GaussianCoordsDisplacement):
    """Take a uniform random step in a ``bdim`` dimensional hypersphere of radius ``stepsize``
    
    this class is the Python interface for the c++ GaussianCoordsDisplacement implementation.
    Takes a step by sampling uniformly a ``bdim`` dimensional hypersphere or radius ``stepsize``
    
    Parameters
    ----------
    rseed : pos int
        seed for the random number generator (std:library 64 bits Merseene Twister)
    stepsize : double
        size of step in each dimension
    ndim : int
        dimensionality of coordinates array
    """

cdef class _Cdef_SampleGaussian(_Cdef_TakeStep):
    cdef cppGaussianTakeStep* newptr
    def __cinit__(self, rseed, stepsize, origin):
        cdef _pele.Array[double] origin_ = array_wrap_np(origin)
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppSampleGaussian(rseed, stepsize, origin_))
        self.newptr = <cppGaussianTakeStep*> self.thisptr.get()
    
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
        
    def get_count(self):
        """get the total count of the number of steps taken
        
        Returns
        -------
        int
            total count of steps taken
        """
        return self.newptr.get_count()
    
    def get_stepsize(self):
        """get the step size
        
        Returns
        -------
        double
            stepsize
        """
        return self.newptr.get_stepsize()

class SampleGaussian(_Cdef_SampleGaussian):
    """Sample directly a Gaussian centered at ``origin`` with standard deviation ``stepsize``
    
    this class is the Python interface for the c++ SampleGaussian implementation.
    Sample directly a Gaussian centered at ``origin`` with standard deviation ``stepsize`` 
    
    Parameters
    ----------
    rseed : pos int
        seed for the random number generator (std:library 64 bits Merseene Twister)
    stepsize : double
        standard deviation for direct sampling
    origin : numpy.array
        coordinates where the gaussian should be centered
    """
    
#
# ParticlePairSwap
#

cdef class _Cdef_ParticlePairSwap(_Cdef_TakeStep):
    cdef cppParticlePairSwap* newptr
    def __cinit__(self, seed, nr_particles):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppParticlePairSwap(seed, nr_particles))
        self.newptr = <cppParticlePairSwap*> self.thisptr.get()

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

class ParticlePairSwap(_Cdef_ParticlePairSwap):
    """Swap a pair of particles
    
    Python interface for c++ ParticlePairSwap
    
    Parameters
    ----------
    seed : pos integer
        Seed for random number generator.
    nr_particles : pos integer
        Number of particles.
    swap_every : pos integer
        Spacing for swapping attempts: particle pair swap is attempted every
        ``swap_every`` move.
    """
    
#
# TakeStepPattern
#

cdef class _Cdef_TakeStepPattern(_Cdef_TakeStep):
    cdef cppTakeStepPattern* newptr
    def __cinit__(self):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppTakeStepPattern())
        self.newptr = <cppTakeStepPattern*> self.thisptr.get()
    
    def add_step(self, _Cdef_TakeStep step, nr_repetitions):
        """add a step to a pattern
        
        Parameters
        ----------
        step : :class:`TakeStep`
            object of class :class:`TakeStep` constructed beforehand
        nr_repetitions: int
            number of Monte Carlo iterations in a row for which this
            move is performed
        """
        self.newptr.add_step(step.thisptr, nr_repetitions)

class TakeStepPattern(_Cdef_TakeStepPattern):
    """Takes multiple steps in a repeated pattern
    
    Python interface for c++ TakeStepPattern. This move
    takes multiple steps in a repeated deterministic pattern.
    
    .. warning:: breaks detailed balance locally
    """
    
#
# TakeStepProbabilities
#

cdef class _Cdef_TakeStepProbabilities(_Cdef_TakeStep):
    cdef cppTakeStepProbabilities* newptr
    def __cinit__(self, seed):
        self.thisptr = shared_ptr[cppTakeStep](<cppTakeStep*> new cppTakeStepProbabilities(seed))
        self.newptr = <cppTakeStepProbabilities*> self.thisptr.get()
    def add_step(self, _Cdef_TakeStep step, weight):
        """add a step to a pattern
        
        all the weights are combined in a normalised discrete distribution
        
        Parameters
        ----------
        step : :class:`TakeStep`
            object of class :class:`TakeStep` constructed beforehand
        weight: double or int
            weight to assign to each move
        """
        self.newptr.add_step(step.thisptr, weight)
        
class TakeStepProbabilities(_Cdef_TakeStepProbabilities):
    """Takes multiple steps in a repeated pattern
    
    Python interface for c++ TakeStepProbabilities. This move
    takes multiple steps, each with some probability thus
    not affecting the detailed balance condition.
    
    .. note:: it does NOT break detailed balance
              hence it is the recommended choice
    """
