# distutils: language = c++
import numpy as np
cimport numpy as np
cimport pele.potentials._pele as _pele
cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport *
import abc
from pele.optimize import Result

cdef class _Cdef_MC(_Cdef_BaseMC):
    
    # these are stored so that the memory is not freed
    cdef public potential
    cdef public niter
    cdef public temperature
    cdef public start_coords
    cdef public stepsize

    def __cinit__(self, pot, coords, temp, pstepsize, pniter, *args, **kwargs):
        cdef np.ndarray[double, ndim=1] start_coords = np.array(coords, dtype=float)        
        cdef _pele.BasePotential potential = pot
        cdef temperature = temp
        cdef niter = pniter
        cdef stepsize = pstepsize
        
        self.thisptr = <cppMC*>new cppMC(potential.thisptr, _pele.Array[double](<double*> start_coords.data, start_coords.size), 
                                         temperature, stepsize)
        
        # these are stored so that the memory is not freed
        self.potential = potential 
        self.start_coords = start_coords
        self.temperature = temperature
        self.niter = niter
        self.stepsize = stepsize
        
        
    def add_action(self, _Cdef_Action action):
        self.thisptr.add_action(action.thisptr)
    
    def add_accept_test(self, _Cdef_AcceptTest test):
        self.thisptr.add_accept_test(test.thisptr)
    
    def add_conf_test(self, _Cdef_ConfTest test):
        self.thisptr.add_conf_test(test.thisptr)
    
    def add_late_conf_test(self, _Cdef_ConfTest test):
        self.thisptr.add_late_conf_test(test.thisptr)
    
    def set_takestep(self, _Cdef_TakeStep takestep):
        self.thisptr.set_takestep(takestep.thisptr)
        
    def set_coordinates(self, coords, energy):
        cdef np.ndarray[double, ndim=1] ccoords = np.array(coords, dtype=float)
        self.thisptr.set_coordinates(_pele.Array[double](<double*> ccoords.data, ccoords.size), energy)
    
    def set_temperature(self, T):
        self.temperature = T
        self.thisptr.set_temperature(T)
    
    def reset_energy(self):
        self.thisptr.reset_energy()
    
    def get_energy(self):
        energy = self.thisptr.get_energy()
        return energy
    
    @cython.boundscheck(False)
    def get_coords(self):
        """return a histogram array"""
        cdef _pele.Array[double] xi = self.thisptr.get_coords()
        cdef double *xdata = xi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros(xi.size())
        cdef size_t i
        for i in xrange(xi.size()):
            x[i] = xdata[i]
              
        return x
    
    def get_norm_coords(self):
        return self.thisptr.get_norm_coords()
    
    def get_accepted_fraction(self):
        accepted_frac = self.thisptr.get_accepted_fraction()
        return accepted_frac
    
    def get_iterations_count(self):
        n = self.thisptr.get_iterations_count()
        return n
    
    def get_conf_rejection_fraction(self):
        f = self.thisptr.get_conf_rejection_fraction()
        return f
    
    def get_E_rejection_fraction(self):
        f = self.thisptr.get_E_rejection_fraction()
        return f
    
    def get_neval(self):
        neval = self.thisptr.get_neval()
        return neval
    
    def get_stepsize(self):
        s = self.thisptr.get_stepsize()
        return s
    
    def one_iteration(self):
        self.thisptr.one_iteration()
    
    def run(self):
        self.thisptr.run(self.niter)
    
    def __reduce__(self):
        return (_Cdef_MC,(self.potential, self.start_coords, self.temperature, self.stepsize, self.niter))

class _BaseMCRunner(_Cdef_MC):
    """
    Abstract method for MC runners, all MC runners should derive from this base class.
    The design of this class relies on a number of implementation choices made for the
    pele::MC cpp class. This is not limiting by any means, developers can easily modify
    this class to write a base class that uses a different basic MC class. 
    Using the pele::MC class is however recommended. 
    *potential should be constructed outside of this class and passed
    *coords are the initial coordinates
    *niter is the total number of MC iterations
    *set_control *MUST* be overwritten in any derived class
    """
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, potential, coords, temperature, stepsize, niter):
        super(_BaseMCRunner,self).__init__(potential, coords, temperature, stepsize, niter)
        
        self.ndim = len(coords)
        self.result = Result()
        self.result.message = []
    
    @abc.abstractmethod
    def set_control(self, c):
        """set control parameter, this could be temperature or some other control parameter like stiffness of the harmonic potential"""
    
    def get_config(self):
        """Return the coordinates of the current configuration and its associated energy"""
        coords = self.get_coords()
        energy = self.get_energy()
        return coords, energy
    
    def set_config(self, coords, energy):
        """set current configuration and its energy"""
        self.set_coordinates(coords, energy)
    
    def get_results(self):
        """Must return a result object, generally must contain at least final configuration and energy"""
        res = self.result
        res.coords = self.get_coords()
        res.energy = self.get_energy()
        return res
    
    def get_status(self):
        status = Result()
        status.iteration = self.get_iterations_count()
        status.acc_frac = self.get_accepted_fraction()
        status.conf_reject_frac = self.get_conf_rejection_fraction()
        status.E_reject_frac = self.get_E_rejection_fraction()
        status.step_size = self.get_stepsize()
        status.energy = self.get_energy()
        status.neval = self.get_neval()
        return status