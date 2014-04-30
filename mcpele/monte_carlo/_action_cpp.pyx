# distutils: language = c++
import numpy as np
cimport numpy as np
cimport pele.potentials._pele as _pele
cimport cython
import sys
from libcpp cimport bool as cbool
from _pele_mc cimport cppAction,_Cdef_Action

#===============================================================================
# Adjust step size
#===============================================================================

cdef extern from "mcpele/actions.h" namespace "mcpele": 
    cdef cppclass cppAdjustStep "mcpele::AdjustStep":
        cppAdjustStep(double, double, size_t, size_t) except +
        
cdef class _Cdef_AdjustStep(_Cdef_Action):
    """This class is the python interface for the c++ pele::AdjustStep action class implementation
    """
    def __cinit__(self, target, factor, niter, navg):
        self.thisptr = <cppAction*>new cppAdjustStep(target, factor, niter, navg)
        
class AdjustStep(_Cdef_AdjustStep):
    """This class is the python interface for the c++ AdjustStep implementation.
    """

#===============================================================================
# Record Energy Histogram
#===============================================================================        

cdef extern from "mcpele/actions.h" namespace "mcpele":
    cdef cppclass cppRecordEnergyHistogram "mcpele::RecordEnergyHistogram":
        cppRecordEnergyHistogram(double, double, double, size_t) except +
        _pele.Array[double] get_histogram() except +
        void print_terminal(size_t) except +
        double get_max() except +
        double get_min() except +
    cdef cppclass cppRecordEnergyTimeseries "mcpele::RecordEnergyTimeseries":
        cppRecordEnergyTimeseries(const size_t) except +
        _pele.Array[double] get_time_series() except +
        
cdef class _Cdef_RecordEnergyHistogram(_Cdef_Action):
    """This class is the python interface for the c++ pele::RecordEnergyHistogram acceptance test class implementation
    """
    def __cinit__(self, min, max, bin, eqsteps):
        self.thisptr = <cppAction*>new cppRecordEnergyHistogram(min, max, bin, eqsteps)
    
    @cython.boundscheck(False)
    def get_histogram(self):
        """return a histogram array"""
        cdef cppRecordEnergyHistogram* newptr = <cppRecordEnergyHistogram*> self.thisptr
        cdef _pele.Array[double] histi = newptr.get_histogram()
        cdef double *histdata = histi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] hist = np.zeros(histi.size())
        cdef size_t i
        for i in xrange(histi.size()):
            hist[i] = histdata[i]
              
        return hist
        
    def print_terminal(self, ntot):
        cdef cppRecordEnergyHistogram* newptr2 = <cppRecordEnergyHistogram*> self.thisptr
        newptr2.print_terminal(ntot)
    
    def get_bounds_val(self):
        cdef cppRecordEnergyHistogram* newptr3 = <cppRecordEnergyHistogram*> self.thisptr
        Emin = newptr3.get_min()
        Emax = newptr3.get_max()
        return Emin, Emax
        
class RecordEnergyHistogram(_Cdef_RecordEnergyHistogram):
    """This class is the python interface for the c++ RecordEnergyHistogram implementation.
    """

#===============================================================================
# RecordEnergyTimeseries
#===============================================================================
        
cdef class _Cdef_RecordEnergyTimeseries(_Cdef_Action):
    """This class is the python interface for the c++ bv::RecordEnergyTimeseries action class implementation
    """
    def __cinit__(self, record_every):
        self.thisptr = <cppAction*>new cppRecordEnergyTimeseries(record_every)
    
    @cython.boundscheck(False)
    def get_time_series(self):
        """return a energy time series array"""
        cdef cppRecordEnergyTimeseries* newptr = <cppRecordEnergyTimeseries*> self.thisptr
        cdef _pele.Array[double] seriesi = newptr.get_time_series()
        cdef double *seriesdata = seriesi.data()
        cdef np.ndarray[double, ndim=1, mode="c"] series = np.zeros(seriesi.size())
        cdef size_t i
        for i in xrange(seriesi.size()):
            series[i] = seriesdata[i]
              
        return series
    
class RecordEnergyTimeseries(_Cdef_RecordEnergyTimeseries):
    """This class is the python interface for the c++ RecordEnergyTimeseries implementation.
    """