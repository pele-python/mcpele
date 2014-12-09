# distutils: language = c++
# distutils: sources = mc.cpp

cimport cython
cimport numpy as np

import numpy as np
import sys
import abc
from pele.optimize import Result

cdef class _Cdef_MC(_Cdef_BaseMC):
    
    # these are stored so that the memory is not freed
    cdef public _pele.BasePotential potential
    cdef public start_coords
    cdef public double temperature
    cdef public size_t niter
    def __init__(self, _pele.BasePotential pot, coords, double temp, size_t pniter):
        cdef np.ndarray[double, ndim=1] cstart_coords = np.array(coords, dtype=float)        
        self.potential = pot 
        self.start_coords = cstart_coords
        self.temperature = temp
        self.niter = pniter
        self.thisptr = shared_ptr[cppMC]( <cppMC*>new cppMC(self.potential.thisptr,
                                                            _pele.Array[double](<double*> cstart_coords.data, cstart_coords.size),
                                                            self.temperature) )
    
    def add_action(self, _Cdef_Action action):
        self.thisptr.get().add_action(action.thisptr)
    
    def add_accept_test(self, _Cdef_AcceptTest test):
        self.thisptr.get().add_accept_test(test.thisptr)
    
    def add_conf_test(self, _Cdef_ConfTest test):
        self.thisptr.get().add_conf_test(test.thisptr)
    
    def add_late_conf_test(self, _Cdef_ConfTest test):
        self.thisptr.get().add_late_conf_test(test.thisptr)
    
    def set_takestep(self, _Cdef_TakeStep takestep):
        self.thisptr.get().set_takestep(takestep.thisptr)
        
    def set_report_steps(self, size_t steps):
        self.thisptr.get().set_report_steps(steps)
        
    def set_coordinates(self, coords, energy):
        cdef np.ndarray[double, ndim=1] ccoords = np.array(coords, dtype=float)
        self.thisptr.get().set_coordinates(_pele.Array[double](<double*> ccoords.data, ccoords.size), energy)
    
    def set_temperature(self, T):
        self.temperature = T
        self.thisptr.get().set_temperature(T)
    
    def get_temperature(self):
        T = self.thisptr.get().get_temperature()
        assert(self.temperature == T) #asserts that python and cpp temperatures match
        return T
        
    def reset_energy(self):
        self.thisptr.get().reset_energy()
    
    def get_energy(self):
        energy = self.thisptr.get().get_energy()
        return energy
    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def get_coords(self):
        """return a histogram array"""
        cdef _pele.Array[double] xi = self.thisptr.get().get_coords()
        cdef double *xdata = xi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] x = np.zeros(xi.size())
        cdef size_t i
        for i in xrange(xi.size()):
            x[i] = xdata[i]
              
        return x
    
    def get_norm_coords(self):
        return self.thisptr.get().get_norm_coords()
    
    def get_accepted_fraction(self):
        accepted_frac = self.thisptr.get().get_accepted_fraction()
        return accepted_frac
    
    def get_iterations_count(self):
        n = self.thisptr.get().get_iterations_count()
        return n
    
    def get_conf_rejection_fraction(self):
        f = self.thisptr.get().get_conf_rejection_fraction()
        return f
    
    def get_E_rejection_fraction(self):
        f = self.thisptr.get().get_E_rejection_fraction()
        return f
    
    def get_neval(self):
        neval = self.thisptr.get().get_neval()
        return neval
    
    def one_iteration(self):
        self.thisptr.get().one_iteration()
    
    def run(self):
        self.thisptr.get().run(self.niter)
    
    def abort(self):
        print "terminating MC, abort called: mc._niter->->infty"
        self.thisptr.get().abort()
    
    def __reduce__(self):
        return (_Cdef_MC,(self.potential, self.start_coords, self.temperature, self.niter))

class _BaseMCRunner(_Cdef_MC):
    """Abstract base class for MC runners, all MC runners should derive from this base class.
    
    Parameters
    ----------
    potential : pele::BasePotential 
        the potential (or cost function) return energy, gradient and hessian
        information given a set of coordinates
    coords : numpy.array
        these are the initial coordinates for the system
    niter : int
        total number of MC iterations to perform
       
    Attributes
    ----------
    
    
    
    .. note:: The design of this class relies on a number of implementation choices made for the
              pele::MC cpp class. This is not limiting by any means, developers can easily modify
              this class to write a base class that uses a different basic MC class. 
              Using the pele::MC class is however recommended.
              
    .. warning:: the cython wrapper to the pele::MC class **demands** that the first 4 parameters
                 of all inheriting classes constructors be positional as in the parent class.
                 Hence pay particular attention to the first 4 positional arguments when
                 constructing a MCrunner class! A workaround could be a pure python class
                 that has a member of type :class:`_BaseMCrunner`, then what you do with the
                 constructor will not matter as long as the member is constructed correctly
    """
    __metaclass__ = abc.ABCMeta
    #super(Metropolis_MCrunner, self).__init__(potential, coords, temperature, niter)
    def __init__(self, potential, coords, temperature, niter):
        super(_BaseMCRunner, self).__init__(potential, coords, temperature, niter)
        
        self.ndim = len(coords)
        self.result = Result()
        self.result.message = []
    
    @abc.abstractmethod
    def set_control(self, c):
        """
        sets the temperature or whichever control parameter that takes
        the role of the temperature, such as the stifness of an harmonic
        sprint to which the system is coupled. This abstract method **must** 
        be overwritten in any derived class.
        """
    
    def get_config(self):
        """
        Return the coordinates of the current configuration and its associated energy
        """
        coords = self.get_coords()
        energy = self.get_energy()
        return coords, energy
    
    def set_config(self, coords, energy):
        """sets current configuration and its energy
        
        Parameters
        ----------
        coords : numpy.array
            these are the initial coordinates for the system
        energy : double
            energy of coords
        """
        self.set_coordinates(coords, energy)
    
    def get_results(self):
        """Must return a result object, generally must contain at least final configuration and energy
        
        Returns
        -------
        res : pele::Result container
            tuple containing typically coords and energy
            accessible via: 
            
            * res.coords
            * res.energy
        """
        res = self.result
        res.coords = self.get_coords()
        res.energy = self.get_energy()
        return res
    
    def get_status(self):
        """Returns typical information about the status of the Monte Carlo walker"""
        status = Result()
        status.iteration = self.get_iterations_count()
        status.acc_frac = self.get_accepted_fraction()
        status.conf_reject_frac = self.get_conf_rejection_fraction()
        status.E_reject_frac = self.get_E_rejection_fraction()
        status.energy = self.get_energy()
        status.neval = self.get_neval()
        return status
    
